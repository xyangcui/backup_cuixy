load "/home/sunming/data5/cuixy/script/ncl_func/monte_caro.ncl"
;--------------------------------------------------------------------------
;function: 计算PCC.
;计算每个事件，每个位相和每一天 具体的PCC。
;--------------------------------------------------------------------------
function calc_pc(x,w)
local x,w,dims,dim0,pct,pc1
begin
  dims = dimsizes(x)
  dim0 = dims(0)
  pct  = new((/dim0,dim0/),float)
  pct  = 0.
  do i = 0, dim0-2
    do j = i+1, dim0-1
      pct(i,j) = pattern_cor(x(i,:,:), x(j,:,:), w, 0)
    end do
  end do
  pc1 = ndtooned(pct)
  return dimsizes(ind(pc1.ge.0.5))
end
;--------------------------------------------------------------------------
;function: 计算每种事件并返回矩阵的程序。
;--------------------------------------------------------------------------
function calc_pro(timt:integer,is,ie,isS,ieS,zt:float,pht:integer,clatt:float,nran)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,igood,pt,ip,out,pct1,nran
begin
  dims = dimsizes(zt)
  dim2 = dims(1)
  dim3 = dims(2)
  nrand = nran
  pct1 = new((/nrand,8,31/),float) 
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
;;IOfP
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ndaf = ndtooned(nda)
  igood = ind(.not.ismissing(ndaf))
  ndaf := ndaf(igood)  ;index for z.
  ptf = pht(ndaf)
  delete([/igood,nda,its,ite,ide/])
;;IOsP
  nda = new((/dimsizes(isS),100/),integer)
  do i = 0,dimsizes(isS)-1
    its = ind(timt.eq.isS(i))
    ite = ind(timt.eq.ieS(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ndas = ndtooned(nda)
  igood = ind(.not.ismissing(ndas))
  ndas := ndas(igood)  ;index for z.
  pts = pht(ndas)
  delete([/igood,nda,its,ite,ide/])

  do i = 0,nrand-1

  do m = 0, 7
    np = m+1
    if (np.eq.8) then
      ipf = ind(ptf.eq.8.or.ptf.eq.1)
      ips = ind(pts.eq.8.or.pts.eq.1)
    else
      ipf = ind(ptf.eq.np.or.ptf.eq.np+1)
      ips = ind(pts.eq.np.or.pts.eq.np+1)
    end if

    outf = new((/dimsizes(ipf),31,dim2,dim3/),float)
    outs = new((/dimsizes(ips),31,dim2,dim3/),float)
    do j = 0,nn
      outf(:,j,:,:) = zt(ndaf(ip)+j,:,:)
      outs(:,j,:,:) = zt(ndas(ip)+j,:,:)
    end do

    pctf = new((/31/),float) 
    pcts = new((/31/),float) 

    do l = 0, nn
      ;random seed.
      rseed1 = toint(systemfunc(" date +%s"))
      rseed2 = toint((12345l*rseed1)%2147483398l)
      random_setallseed(rseed1, rseed2)
    
      iwf = generate_sample_indices(dimsizes(outf(:,0,0)), 1)
      iws = generate_sample_indices(dimsizes(outs(:,0,0)), 1)

      pctf(l) = calc_pc(outf(iwf,l,:,:),clatt)
      pcts(l) = calc_pc(outs(iws,l,:,:),clatt)

      delete([/iwf,iws/])
    end do

    dimf = dimsizes(ipf)
    dims = dimsizes(ips)
    pct1(i,m,:) = pcts/(dims*(dims-1)/2.)-pctf/(dimf*(dimf-1)/2.)
    delete([/np,ipf,ips,outs,outf,pctf,pcts,dimf,dims/])
  end do

  return pct1

end
;----------------------------------------------------------------------
  lats = 20
  latn = 70
  lonw = 150
  lone = 240

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  zdir = ddir+"/"+"ERA5_z500_daily_anom_1979-2022_r2.5.nc"
;;jet events
  ndir = ddir+"/"+"IOnP_date_omi.csv"
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"
  idir = ddir+"/"+"IOiP_date_omi.csv"

  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
;-----------------------------------------------------------
; main code.
;-----------------------------------------------------------
;;times
  f = addfile(zdir,"r")
  z = f->z(:,{lats:latn},{lonw:lone})
  time = z&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.10.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.5.and.day.eq.31)

  mon := mon(n1:n2)
  z := z(n1:n2,:,:)
  nt = ind(mon.eq.10.or.mon.eq.11.or.mon.eq.12.or.mon.eq.1\
           .or.mon.eq.2.or.mon.eq.3.or.mon.eq.4.or.mon.eq.5)
  z := z(nt,:,:)
  delete([/TIME,time,year,mon,day,nt/])
;; ONDJFMAM
  time = z&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  timi = cd_calendar(time,-2)
;;OMI
  delim = ","
  rmmf  = asciiread(omidir,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
  phaset = str_get_field(rmmf,4,delim)
  ampt = str_get_field(rmmf,7,delim)

  iyear := toint(iyear(1::))
  imon  := toint(imon(1::))
  iday  := toint(iday(1::))
  phaset := tointeger(phaset(1::))
  ampt := tofloat(ampt(1::))

  phase = new((/dimsizes(time)/),integer)
  amp   = new((/dimsizes(time)/),float)
  phase = 0
  amp   = 0.

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i)\
             .and.iday.eq.day(i))
    phase(i) = phaset(nt)
    amp(i) = ampt(nt)
    delete(nt)
  end do
  delete([/phaset,ampt,iyear,imon,iday/])
  print("complete OMI selling "+ systemfunc("date"))
;----------------------------------------------------------------------
;计算PCC
  ind_time = ispan(0,dimsizes(time)-1,1) 
  nrand = 5000
  clat = tofloat(cos(0.01745329*z&lat))
  pc = new((/nrand,8,30+1/),float) ;(events,phases,days)

  nf = asciiread(fdir,-1,"string")
  is  = str_get_field(nf,1,delim)
  ie  = str_get_field(nf,2,delim)
  isf = toint(is(1::))
  ief = toint(ie(1::))
  delete([/is,ie/])

  ns  = asciiread(sdir,-1,"string")
  is  = str_get_field(ns,1,delim)
  ie  = str_get_field(ns,2,delim)
  iss = toint(is(1::))
  ies = toint(ie(1::))

  print("random calc begin. "+ systemfunc("date"))
  pc = calc_pro(timi,isf,ief,iss,ies,z,phase,clat,nrand) ;IOfP
  print("random calc done. "+ systemfunc("date"))

;;store.
  print("store begin. "+ systemfunc("date"))

  odir = ddir+"/"+"boostrap_pc_omi_days_2phase.nc"
  system("/bin/rm -f "+odir)   
  fout = addfile(odir,"c")  

  fAtt = True           
  fAtt@creation_date = systemfunc ("date")

  fileattdef( fout, fAtt )            
  fout->pc = pc   

  print("store done. "+ systemfunc("date"))  