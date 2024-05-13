load "/home/sunming/data5/cuixy/script/ncl_func/monte_caro.ncl"
;--------------------------------------------------------
; time: 2024-12-17
;--------------------------------------------------------
;-----------------------------------------------------------
  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
  outdir = ddir+"/"+"monte_carlo/Fig13_uv_p5_sig_day-5.nc"
  out2dir = ddir+"/"+"monte_carlo/Fig13_olr_p5_sig_day-5.nc"

  olrdir = ddir+"/"+"olr_daily_anom_1979-2022.nc"

  udir = "/home/sunming/data5/cuixy/global_jets/data/ERA5_daily_iso_uvert_1979-2022_r2.5.nc"
  vdir = "/home/sunming/data5/cuixy/global_jets/data/ERA5_daily_iso_vvert_1979-2022_r2.5.nc"

  lats = -10
  latn =  50
  lonw =  90
  lone =  210
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
  llvalue = new((/5,dim1,dim2/),float)
  llvalue = 0.

  nt = ind(phaset.eq.5.and.ampt.gt.1.)
  nt1 = in22(nt)
  do j = 0, 4
    ntt   = -7 + j*1
    nindd = nt1 + ntt
    llvalue(j,:,:) = dim_avg_n(value(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do

  llv = dim_avg_n(llvalue,0)

  return llv

end

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  folr = addfile(olrdir,"r")

  olr = folr->olr(:,{lats:latn},{lonw:lone})

  fu = addfile(udir,"r")
  u  = fu->u(:,{lats:latn},{lonw:lone})

  fv = addfile(vdir,"r")
  v  = fv->v(:,{lats:latn},{lonw:lone})

  time = olr&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  print("complete data selling "+ systemfunc("date"))
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
  olr := olr(n1:n2,:,:)
  u   := u(n1:n2,:,:)
  v   := v(n1:n2,:,:)
  reamp   := reamp(n1:n2)
  rephase := rephase(n1:n2)
  n_ndjfm = ind(mont.eq.11.or.mont.eq.12.or.mont.eq.1.or.mont.eq.2.or.mont.eq.3)
  olr     := olr(n_ndjfm,:,:)  ; NDJFM can replace by variable.
  u       := u(n_ndjfm,:,:)    ; NDJFM can replace by variable.
  v       := v(n_ndjfm,:,:)    ; NDJFM can replace by variable.
  reamp   := reamp(n_ndjfm)
  rephase := rephase(n_ndjfm)

  delete([/time,TIME,year,mon,day/])

  time = olr&time
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
  nrandom = 10000

  m_value_2 = new((/nrandom,2,dimsizes(u&lat),dimsizes(u&lon)/),float)
  m_value_2 = 0.

  m_value_3 = new((/nrandom,dimsizes(olr&lat),dimsizes(olr&lon)/),float)
  m_value_3 = 0.

  print("random calc u begin. "+ systemfunc("date"))
  do k = 0, nrandom - 1
    rd_time = monte_carlo1(ind_time)
    ut   = u(rd_time,:,:)
    vt   = v(rd_time,:,:)
    olrt = olr(rd_time,:,:)
    m_value_3(k,:,:) = calc_lagvalue(in2,phase,amp,olrt)
    m_value_2(k,0,:,:) = calc_lagvalue(in2,phase,amp,ut)
    m_value_2(k,1,:,:) = calc_lagvalue(in2,phase,amp,vt)
    delete([/ut,vt,olrt,rd_time/])
  end do

  copy_VarCoords(u(0,:,:), m_value_2(0,0,:,:))
  copy_VarCoords(olr(0,:,:), m_value_3(0,:,:))

  print("random calc done. "+ systemfunc("date"))

; sort values for 0.05, 0.1, 0.9,0.95
  sigl0  = 0.05
  sigl1  = 0.1
  sigh0  = 0.9
  sigh1  = 0.95

  m_value2 = dim_pqsort_n(m_value_2, 2, 0)

  m_value  = new((/4,2,dimsizes(u&lat),dimsizes(u&lon)/),float) 
  m_value  = 0.

  nl0  = toint(dimsizes(m_value_2(:,0,0,0))*sigl0) - 1
  nh0  = toint(dimsizes(m_value_2(:,0,0,0))*sigh0) - 1
  nl1  = toint(dimsizes(m_value_2(:,0,0,0))*sigl1) - 1
  nh1  = toint(dimsizes(m_value_2(:,0,0,0))*sigh1) - 1

  m_value(0,:,:,:) = m_value_2(nl0,:,:,:)
  m_value(1,:,:,:) = m_value_2(nl1,:,:,:)
  m_value(2,:,:,:) = m_value_2(nh0,:,:,:)
  m_value(3,:,:,:) = m_value_2(nh1,:,:,:)

  delete([/nl0,nl1,nh0,nh1/])
;;
  m_value3 = dim_pqsort_n(m_value_3, 2, 0)

  m_olr  = new((/4,dimsizes(olr&lat),dimsizes(olr&lon)/),float) 
  m_olr  = 0.

  nl0  = toint(dimsizes(m_value_3(:,0,0))*sigl0) - 1
  nh0  = toint(dimsizes(m_value_3(:,0,0))*sigh0) - 1
  nl1  = toint(dimsizes(m_value_3(:,0,0))*sigl1) - 1
  nh1  = toint(dimsizes(m_value_3(:,0,0))*sigh1) - 1

  m_olr(0,:,:) = m_value_3(nl0,:,:)
  m_olr(1,:,:) = m_value_3(nl1,:,:)
  m_olr(2,:,:) = m_value_3(nh0,:,:)
  m_olr(3,:,:) = m_value_3(nh1,:,:)

  system("/bin/rm -f "+outdir)
  ncdf = addfile(outdir,"c")
  ncdf->m_value = m_value

  system("/bin/rm -f "+out2dir)
  ncdf = addfile(out2dir,"c")
  ncdf->m_olr = m_olr

  print("sort done. "+ systemfunc("date"))