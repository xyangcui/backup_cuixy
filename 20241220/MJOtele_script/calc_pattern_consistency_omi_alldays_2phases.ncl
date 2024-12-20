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
function calc_pro(timt:integer,xdir,zt:float,pht:integer,amt:float,nn:integer,clatt:float)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,igood,pt,at,ip,out,pct1
begin
  dims = dimsizes(zt)
  dim2 = dims(1)
  dim3 = dims(2)
  pct1 = new((/8,nn+1/),float) 
  pct1 = 0.

  delim = ","
  nint  = ispan(0,dimsizes(timt),1) ; time index.
  nf  = asciiread(xdir,-1,"string")
  ist = str_get_field(nf,1,delim)
  iet = str_get_field(nf,2,delim)
  is  = toint(ist(1::))
  ie  = toint(iet(1::))

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
;select amp >1. and phase. 2 phases
  do m = 0, 7
    np = m+1
    print("phase "+ np +" begin.")
    if (np.eq.8) then
      ip = ind(pt.eq.8.or.pt.eq.1)
    else
      ip = ind(pt.eq.np.or.pt.eq.np+1)
    end if
    ip@_FillValue = -999
    if (any(ip.eq.ip@_FillValue).eq."Missing")
      continue
    else
      out = new((/dimsizes(ip),31,dim2,dim3/),float)
      do j = 0,nn
        out(:,j,:,:) = zt(nda(ip)+j,:,:)
      end do
      do l = 0, nn
        pct1(m,l) = calc_pc(out(:,l,:,:),clatt)
      end do
      dim1 = dimsizes(ip)
      pct1(m,:) = pct1(m,:)/(dim1*(dim1-1)/2.)
      delete([/ip,out/])
      print("phase "+ np +" done.")
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
  delim  = ","
  rmmf   = asciiread(omidir,-1,"string")
  iyear  = str_get_field(rmmf,1,delim)
  imon   = str_get_field(rmmf,2,delim)
  iday   = str_get_field(rmmf,3,delim)
  phaset = str_get_field(rmmf,4,delim)
  ampt   = str_get_field(rmmf,7,delim)

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
  n = 30
  clat = tofloat(cos(0.01745329*z&lat))
  pc = new((/2,8,n+1/),float) ;(events,phases,days)
  pc = 0.
;;jet events (year,mon,day) (events,days,lat,lon)
  print("IOfP begin. "+ systemfunc("date"))
  pc(0,:,:) = calc_pro(timi,fdir,z,phase,amp,n,clat) ;IOfP
  print("IOfP done. "+ systemfunc("date"))

  print("IOsP begin. "+ systemfunc("date"))
  pc(1,:,:) = calc_pro(timi,sdir,z,phase,amp,n,clat) ;IOsP
  print("IOsP done. "+ systemfunc("date"))

  pc@dim1 = "IOfP,IOsP"
  pc@dim2 = "P12,P23,P34,P45,P56,P67,P78,P81"
  pc@dim3 = "lead days"

;;store.
  print("store begin. "+ systemfunc("date"))

  odir = ddir+"/"+"pattern_consistency_omi_alldays_2phases.nc"
  system("/bin/rm -f "+odir)   
  fout = addfile(odir,"c")  

  fAtt = True           
  fAtt@creation_date = systemfunc ("date")

  fileattdef( fout, fAtt )            
  fout->pc = pc   

  print("store done. "+ systemfunc("date"))  