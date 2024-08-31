load "/home/sunming/data5/cuixy/script/ncl_func/monte_caro.ncl"
;---------------------------------------------------------------
; draw temperature anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour.
; time: 2023-12-17
;---------------------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  idir = ddir + "/ERA5_daily_anom_ulevel_1979-2022_r2.5.nc"
  indir = ddir+"/mhf_events.csv"

  outdir = "/home/sunming/data5/cuixy/global_jets/data"+"/"+"monte_carlo/st_events_sig.nc"
  lats =  10
  latn =  80
  lonw = -60
  lone =  30
;-------------------------------------------------------------
; function calc_lagvalue.
; note that this function should change follow your intention.
;-------------------------------------------------------------
function calc(va,py,pm,pd,yeart,mont,dayt)
local va,py,pm,pd
begin
  DIMS = dimsizes(va)
  dim0 = dimsizes(py)
  dim1 = DIMS(1)
  dim2 = DIMS(2)
  llva = new((/dim0,dim1,dim2/),float)
  llva = 0.

  do i = 0, dim0-1
    nt = ind(yeart.eq.py(i).and.mont.eq.pm(i).and.dayt.eq.pd(i))
    llva(i,:,:) = dim_avg_n(va(nt+1:nt+5,:,:),0)
    delete(nt)
  end do

  llv = dim_avg_n(llva,0)

  return llv
end
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  f = addfile(idir,"r")
  x = f->u(:,{250},{lats:latn},{lonw:lone})

  time = x&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  print("complete data selling "+ systemfunc("date"))
;-------------------------------------------------------
;; selling events.
  delim  = ","
  rmmf   = asciiread(indir,-1,"string")
  iyear  = str_get_field(rmmf,1,delim)
  imon   = str_get_field(rmmf,2,delim)
  iday   = str_get_field(rmmf,3,delim)
  pct = str_get_field(rmmf,4,delim)

  iyear  := toint(iyear(1::))
  imon   := toint(imon(1::))
  iday   := toint(iday(1::))
  pct := toint(pct(1::))

  pev   = ind(pct.gt.0)
  pyear = iyear(pev)
  pmon  = imon(pev)
  pday  = iday(pev)

  nev   = ind(pct.lt.0)
  nyear = iyear(nev)
  nmon  = imon(nev)
  nday  = iday(nev)
;------------------------------------------------------------
; sort NDJFM.
  n1 = ind(year.eq.1979.and.mon.eq.10.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.4.and.day.eq.1)
  mont := mon(n1:n2)
  x   := x(n1:n2,:,:)
  n_djf = ind(mont.eq.12.or.mont.eq.1.or.mont.eq.2)
  x    := x(n_djf,:,:)  ; DJF can replace by variable.
  delete([/time,TIME,year,mon,day/])

  time = x&time
  TIME = cd_calendar(time, 0)
  year = toint( TIME(:,0) )
  mon  = toint( TIME(:,1) )
  day  = toint( TIME(:,2) )
  print("DJF done. "+ systemfunc("date"))
;-----------------------------------------------------------------------------
;;calculate phase-lag data. lag time 31 days. NDJFM. for shuttle.
  ind_time = ispan(0,dimsizes(time)-1,1) 
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
;; (nrand, pe/ne, lat ,lon)
  nrandom = 10000
  m_value = new((/nrandom,2,dimsizes(x&lat),dimsizes(x&lon)/),float)
  m_value = 0.

  print("random calc begin. "+ systemfunc("date"))
  do k = 0, nrandom - 1
    rd_time = monte_carlo1(ind_time)
    xt = x(rd_time,:,:)
    m_value(k,0,:,:) = calc(xt,pyear,pmon,pday,year,mon,day)
    m_value(k,1,:,:) = calc(xt,nyear,nmon,nday,year,mon,day)
    delete([/xt,rd_time/])
  end do
  copy_VarCoords(x(0,:,:), m_value(0,0,:,:))
  print("random calc done. "+ systemfunc("date"))
; store values.
; store values
  m_value_2 = dim_pqsort_n(m_value, 2, 0)

  mvalue  = new((/4,2,dimsizes(x&lat),dimsizes(x&lon)/),float) 
  mvalue  = 0.
  copy_VarCoords(x(0,:,:), mvalue(0,0,:,:))
  mvalue@long_name = "positive/negative"
  mvalue@store = "(u),0.05,0.95,0.1,0.9"

  sigl = 0.05
  sigh = 0.95
  nl2  = toint(dimsizes(m_value(:,0,0,0))*sigl) - 1
  nh2  = toint(dimsizes(m_value(:,0,0,0))*sigh) - 1
  mvalue(0,:,:,:) = m_value(nl2,:,:,:)
  mvalue(1,:,:,:) = m_value(nh2,:,:,:)

  sigl = 0.1
  sigh = 0.9
  nl2  = toint(dimsizes(m_value(:,0,0,0))*sigl) - 1
  nh2  = toint(dimsizes(m_value(:,0,0,0))*sigh) - 1
  mvalue(2,:,:,:) = m_value(nl2,:,:,:)
  mvalue(3,:,:,:) = m_value(nh2,:,:,:)

  system("/bin/rm -f "+outdir)
  ncdf = addfile(outdir,"c")
  ncdf->mvalue = mvalue
  print("store done. "+ systemfunc("date"))