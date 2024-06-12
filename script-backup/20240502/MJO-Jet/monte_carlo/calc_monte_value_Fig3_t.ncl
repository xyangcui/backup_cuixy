load "/home/sunming/data5/cuixy/script/ncl_func/monte_caro.ncl"
;--------------------------------------------------------
; draw temperature anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour.
; time: 2023-12-17
;--------------------------------------------------------
;-----------------------------------------------------------
  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  phidir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022_2.5x2.5.nc"
  tdir   = ddir+"/"+"ERA5_daily_anom_tvert_1979-2022_2.5x2.5.nc"

  outdir = ddir+"/"+"monte_carlo/Fig3_t_sig.nc"

  lats =  10
  latn =  55
  lonw =   0
  lone = 150
;--------------------------------------------------------
; function calc_lagvalue.
;--------------------------------------------------------
function calc_lagvalue(in22,phaset,ampt,value)
local in22,phaset,ampt,value,DIMS,dim1,dim2,llvalue,\
ip,nt,ntt,nt1,nindd
begin
  DIMS = dimsizes(value)
  dim1 = DIMS(1)
  dim2 = DIMS(2)
  llvalue = new((/dim1,dim2/),float)
  llvalue = 0.

  nt = ind(phaset.eq.3.and.ampt.gt.1.)
  nind = in22(nt)
  ntt   = 5
  nindd = nind + ntt
  llvalue = dim_avg_n(value(nindd,:,:), 0)
  delete([/nt,ntt,nindd,nind/])

  return llvalue

end

;--------------------------------------------------------
; function calc_lagvalue.
;--------------------------------------------------------
function calc_lagvalue2(in22,phaset,ampt,value)
local in22,phaset,ampt,value,DIMS,dim1,dim2,llvalue,\
ip,nt,ntt,nt1,nindd
begin
  DIMS = dimsizes(value)
  dim1 = DIMS(1)
  dim2 = DIMS(2)
  llvalue = new((/dim1,dim2/),float)
  llvalue = 0.

  nt = ind(phaset.eq.7.and.ampt.gt.1.)
  nind = in22(nt)
  ntt   = 5
  nindd = nind + ntt
  llvalue = dim_avg_n(value(nindd,:,:), 0)
  delete([/nt,ntt,nindd,nind/])

  return llvalue

end

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------
  ft   = addfile(tdir,"r")
  t    = ft->t(:,{lats:latn},{lonw:lone})

  print("complete data selling "+ systemfunc("date"))

  time = t&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

;-------------------------------------------------------
;; MJO phase and amp.
  RMMf  = rmmdir
  delim = ","
  rmmf  = asciiread(RMMf,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
  phase = str_get_field(rmmf,6,delim)
  amp   = str_get_field(rmmf,7,delim)

  iyear := toint(iyear(1::))
  imon  := toint(imon(1::))
  iday  := toint(iday(1::))
  phase := toint(phase(1::))
  amp   := tofloat(amp(1::))

  reamp   = new((/dimsizes(time)/),float)
  rephase = new((/dimsizes(time)/),integer)
  reamp   = 0.
  rephase = 0

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
    reamp(i)   = amp(nt)
    rephase(i) = phase(nt)
    delete(nt)
  end do

  delete([/amp,phase/])
  print("sort MJO done. "+systemfunc("date"))

;------------------------------------------------------------
; sort NDJFM.
  n1 = ind(year.eq.1979.and.mon.eq.10.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.4.and.day.eq.1)
  mont = mon(n1:n2)
  t := t(n1:n2,:,:)
  reamp   := reamp(n1:n2)
  rephase := rephase(n1:n2)
  n_ndjfm = ind(mont.eq.11.or.mont.eq.12.or.mont.eq.1.or.mont.eq.2.or.mont.eq.3)
  t  := t(n_ndjfm,:,:)  ; NDJFM can replace by variable.
  reamp   := reamp(n_ndjfm)
  rephase := rephase(n_ndjfm)

  delete([/time,TIME,year,mon,day/])

  time = t&time
  TIME = cd_calendar(time, 0)
  year = toint( TIME(:,0) )
  mon  = toint( TIME(:,1) )
  day  = toint( TIME(:,2) )

  print("NDJFM done. "+ systemfunc("date"))
;-----------------------------------------------------------------------------
;;calculate phase-lag data. lag time 31 days. NDJFM. for shuttle.
;;first sel data DJF.  

  n_djf = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  in1   = ispan(0,dimsizes(time)-1,1)   ; NDJFM.
  in2   = in1(n_djf)  ; DJF indices
  amp   = reamp(n_djf)
  phase = rephase(n_djf)

  ind_time = ispan(0,dimsizes(time)-1,1) 
  print("ind_time done. "+ systemfunc("date"))
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
;; (in22,phaset,ampt,value)
  nrandom = 10000

  m_value_2 = new((/nrandom,dimsizes(t&lat),dimsizes(t&lon)/),float)
  m_value_2 = 0.
  m_value_1 = m_value_2
  print("random calc t begin. "+ systemfunc("date"))
  do k = 0, nrandom - 1
    rd_time = monte_carlo1(ind_time)
    tt   = t(rd_time,:,:)
    m_value_2(k,:,:) = calc_lagvalue(in2,phase,amp,tt)
    m_value_1(k,:,:) = calc_lagvalue2(in2,phase,amp,tt)
    delete([/tt,rd_time/])
  end do
  print("random calc done. "+ systemfunc("date"))
  copy_VarCoords(t(0,:,:), m_value_2(0,:,:))
  copy_VarCoords(t(0,:,:), m_value_1(0,:,:))

; sort values
  m_value2 = dim_pqsort_n(m_value_2, 2, 0)
  m_value1 = dim_pqsort_n(m_value_1, 2, 0)

  mvalue1  = new((/4,dimsizes(t&lat),dimsizes(t&lon)/),float) 
  mvalue1  = 0.
  mvalue2  = mvalue1

  copy_VarCoords(t(0,:,:), mvalue2(0,:,:))
  copy_VarCoords(t(0,:,:), mvalue1(0,:,:))
  mvalue1@long_name = "phase 7"
  mvalue2@long_name = "phase 3"
  mvalue1@store = "(t),0.05,0.95,0.1,0.9"
  mvalue2@store = "(t),0.05,0.95,0.1,0.9"

  sigl = 0.05
  sigh = 0.95

  nl2  = toint(dimsizes(m_value_2(:,0,0))*sigl) - 1
  nh2  = toint(dimsizes(m_value_2(:,0,0))*sigh) - 1
  nl1  = toint(dimsizes(m_value_1(:,0,0))*sigl) - 1
  nh1  = toint(dimsizes(m_value_1(:,0,0))*sigh) - 1
  mvalue2(0,:,:) = m_value_2(nl2,:,:)
  mvalue2(1,:,:) = m_value_2(nh2,:,:)
  mvalue1(0,:,:) = m_value_1(nl1,:,:)
  mvalue1(1,:,:) = m_value_1(nh1,:,:)

  sigl = 0.1
  sigh = 0.9

  nl2  = toint(dimsizes(m_value_2(:,0,0))*sigl) - 1
  nh2  = toint(dimsizes(m_value_2(:,0,0))*sigh) - 1
  nl1  = toint(dimsizes(m_value_1(:,0,0))*sigl) - 1
  nh1  = toint(dimsizes(m_value_1(:,0,0))*sigh) - 1
  mvalue2(2,:,:) = m_value_2(nl2,:,:)
  mvalue2(3,:,:) = m_value_2(nh2,:,:)
  mvalue1(2,:,:) = m_value_1(nl1,:,:)
  mvalue1(3,:,:) = m_value_1(nh1,:,:)

  system("/bin/rm -f "+outdir)
  ncdf = addfile(outdir,"c")
  ncdf->mvalue1 = mvalue1
  ncdf->mvalue2 = mvalue2

  print("store done. "+ systemfunc("date"))