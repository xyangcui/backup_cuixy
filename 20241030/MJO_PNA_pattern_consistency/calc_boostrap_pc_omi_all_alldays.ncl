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
function calc_pro(timt:integer,is,ie,zt:float,pht:integer,amt:float,clatt:float)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,igood,pt,at,ip,out,pct1
begin
  dims = dimsizes(zt)
  dim2 = dims(1)
  dim3 = dims(2)
  pct1 = new((/8/),float) 
  pct1 = 0.
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ;event all days.
  nda := ndtooned(nda)
  igood = ind(.not.ismissing(nda))
  nda := nda(igood)  ;index for z.
  pt = pht(nda)
  at = amt(nda)
  delete(igood)
  ;select amp >1. and phase.
  rseed1 = toint(systemfunc(" date +%s"))
  rseed2 = toint((12345l*rseed1)%2147483398l)+toint(random_normal(100, 10, 1))
  random_setallseed(rseed1, rseed2)
  do m = 0, 7
    ;random seed.
    rseed1 = toint(systemfunc(" date +%s"))
    rseed2 = toint((12345l*rseed1)%2147483398l)+toint(random_normal(100, 10, 1))
    random_setallseed(rseed1, rseed2)

    np = m+1
    ip = ind(pt.eq.m+1)
    ip@_FillValue = -999
    if (any(ip.eq.ip@_FillValue).eq."Missing")
      continue
    else
    out = zt(nda(ip),:,:)
    iw = generate_sample_indices(dimsizes(out(:,0,0)), 1)
    pct1(m) = calc_pc(out(iw,:,:),clatt)
    dim1 = dimsizes(ip)
    pct1(m) = pct1(m)/(dim1*(dim1-1)/2.)
    delete([/ip,out/])
    end if
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
  nrand = 1000
  clat = tofloat(cos(0.01745329*z&lat))
  pc = new((/nrand,8/),float) ;(events,phases,days)
  pc = 0.

;;jet events (year,mon,day) (events,days,lat,lon)
  nf  = asciiread(idir,-1,"string")
  is  = str_get_field(nf,1,delim)
  ie  = str_get_field(nf,2,delim)
  ist = toint(is(1::))
  iet = toint(ie(1::))

  nf := asciiread(ndir,-1,"string")
  is := str_get_field(nf,1,delim)
  ie := str_get_field(nf,2,delim)

  ist := array_append_record(ist, toint(is(1::)), 0)
  iet := array_append_record(iet, toint(ie(1::)), 0)

  nf := asciiread(fdir,-1,"string")
  is := str_get_field(nf,1,delim)
  ie := str_get_field(nf,2,delim)

  ist := array_append_record(ist, toint(is(1::)), 0)
  iet := array_append_record(iet, toint(ie(1::)), 0)

  nf := asciiread(sdir,-1,"string")
  is := str_get_field(nf,1,delim)
  ie := str_get_field(nf,2,delim) 

  ist := array_append_record(ist, toint(is(1::)), 0)
  iet := array_append_record(iet, toint(ie(1::)), 0)

  print("random calc begin. "+ systemfunc("date"))
  do k = 0, nrand - 1
    pc(k,:) = calc_pro(timi,ist,iet,z,phase,amp,clat) ;IOfP
  end do
  print("random calc done. "+ systemfunc("date"))

;;store.
  print("store begin. "+ systemfunc("date"))

  odir = ddir+"/"+"boostrap_pc_omi_all_alldays.nc"
  system("/bin/rm -f "+odir)   
  fout = addfile(odir,"c")  

  fAtt = True           
  fAtt@creation_date = systemfunc ("date")

  fileattdef( fout, fAtt )            
  fout->pc = pc   

  print("store done. "+ systemfunc("date"))  