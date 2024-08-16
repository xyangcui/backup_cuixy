load "/home/sunming/data5/cuixy/script/ncl_func/monte_caro.ncl"
;---------------------------------------------------------------
; draw temperature anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour.
; time: 2023-12-17
;---------------------------------------------------------------
  lats = 0
  latn = 60
  lonw = 0
  lone = 360

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  udir = ddir+"/u250_daily_anom_1979-2022_10-90_r2.5.nc"
  outdir = ddir+"/"+"monte_carlo/MLR_sig.nc"
;-------------------------------------------------------------
; function calc_lagvalue.
; note that this function should change follow your intention.
;-------------------------------------------------------------
function calc_lagvalue(UC,X)
local UC,X,R2
begin
  DIMS = dimsizes(UC)
  dim1 = DIMS(1)
  dim2 = DIMS(2)
  R2 = new((/dim1,dim2/),float)

  do i = 0,dim1-1
    do j = 0,dim2-1
      b = reg_multlin_stats(UC(:,i,j), X, False)
      R2(i,j) = b@SSR/b@SST
    end do
  end do
  return R2

end
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------
  f = addfile(udir,"r")

  u = f->u(:,{lats:latn},{lonw:lone})
  lat = u&lat
  lon = u&lon 

  time = u&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  print("complete data selling "+ systemfunc("date"))

  time = u&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
;-------------------------------------------------------
;; MJO phase and amp.
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  delim  = ","
  rmmf   = asciiread(rmmdir,-1,"string")
  iyear  = str_get_field(rmmf,1,delim)
  imon   = str_get_field(rmmf,2,delim)
  iday   = str_get_field(rmmf,3,delim)
  rmm1t  = str_get_field(rmmf,4,delim)
  rmm2t  = str_get_field(rmmf,5,delim)

  iyear  := toint(iyear(1::))
  imon   := toint(imon(1::))
  iday   := toint(iday(1::))
  rmm1t  := tofloat(rmm1t(1::))
  rmm2t  := tofloat(rmm2t(1::))

  rmm1 = new((/dimsizes(time)/),float)
  rmm2 = new((/dimsizes(time)/),float)
  rmm1 = 0.
  rmm2 = 0.

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
    rmm1(i) = rmm1t(nt)
    rmm2(i) = rmm2t(nt)
    delete(nt)
  end do

  nn = ispan(0, dimsizes(time)-1, 1)
  rmm1!0 = "timeind"
  rmm2!0 = "timeind"
  rmm1&timeind = nn
  rmm2&timeind = nn
  delete([/nn,rmm1t,rmm2t,iyear,imon,iday/])

  print("complete RMM selling "+ systemfunc("date"))
;----------------------------------------------------------
; prepare data for regression. (select DJF).
  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)
  mon  := mon(n1:n2)
  rmm1 := rmm1(n1:n2)
  rmm2 := rmm2(n1:n2)
  u    := u(n1:n2,:,:)

  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  rmm1 := rmm1(nt)
  rmm2 := rmm2(nt)
  u    := u(nt,:,:)
  mon  := mon(nt)

  amp = sqrt(rmm1*rmm1+rmm2*rmm2)
  nt := ind(amp.gt.1.)
  rmm1 := rmm1(nt)
  rmm2 := rmm2(nt)

  x = new((/dimsizes(nt),2/),float)
  x(:,0) = rmm1
  x(:,1) = rmm2 
  print("ind_time done. "+ systemfunc("date"))
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
;; (in22,phaset,ampt,value)
  ind_time = ispan(0,dimsizes(mon)-1,1)
  nrandom = 500
  m_value = new((/nrandom,dimsizes(u&lat),dimsizes(u&lon)/),float)
  m_value = 0.
  print("random calc begin. "+ systemfunc("date"))

  do k = 0, nrandom - 1
    rd_time = monte_carlo1(ind_time)
    xt = u(rd_time,:,:)
    uc := xt(nt,:,:)
    m_value(k,:,:) = calc_lagvalue(uc,x)
    delete([/xt,rd_time,uc/])
  end do

  print("random calc done. "+ systemfunc("date"))

; store values.
  copy_VarCoords(u(0,:,:), m_value(0,:,:))

  rdtimes = ispan(1, 500, 1)
  m_value!0 = "time"
  m_value&time = rdtimes
  m_valuet = dim_pqsort_n(m_value, 2, 0)

  system("/bin/rm -f "+outdir)
  fout = addfile(outdir,"c")
  fout->R2 = m_value
  print("store done. "+ systemfunc("date"))