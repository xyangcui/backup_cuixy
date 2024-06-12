load "/home/sunming/data5/cuixy/script/ncl_func/monte_caro.ncl"
;--------------------------------------------------------
; draw temperature anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour.
; time: 2023-12-17
;--------------------------------------------------------
;-----------------------------------------------------------
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  ddir = "/home/sunming/data5/cuixy/global_jets/data/meri_circu"
  wdir = ddir + "/anom_"+"w"+"level_"+"EAJS"+"_r2.5.nc"
  vdir = ddir + "/anom_"+"v"+"level_"+"EAJS"+"_r2.5.nc"

  olrdir = "/home/sunming/data5/cuixy/global_jets/data/olr_daily_anom_1979-2022.nc"

  outdir = "/home/sunming/data5/cuixy/global_jets/data"+"/"+"monte_carlo/Fig7_meri_sig.nc"

  latn =  40
  lats = -20
  lonw =  105
  lone =  140
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
  llv = new((/4,dim1,dim2/),float)
  llv = 0.
  nt = ind(phaset.eq.2.and.ampt.gt.1.)
  nt1 = in22(nt)
  do j = 0, 3
    ntt   = -10 + j*3
    nindd = nt1 + ntt
    llv(j,:,:) = dim_avg_n(value(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do
  return llv
end
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------
  fw  = addfile(wdir,"r")
  fv  = addfile(vdir,"r")

  w = fw->w(:,:,{lats:latn},{lonw:lone})
  v = fv->v(:,:,{lats:latn},{lonw:lone})

  w := dim_avg_n_Wrap(w, 3)
  v := dim_avg_n_Wrap(v, 3)
;---------------------------------------------------------------=
;; MJO phases.
  time = fw->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

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
  w := w(n1:n2,:,:)
  v := v(n1:n2,:,:)
  reamp   := reamp(n1:n2)
  rephase := rephase(n1:n2)

  n_ndjfm  = ind(mont.eq.11.or.mont.eq.12.or.mont.eq.1.or.mont.eq.2.or.mont.eq.3)
  w := w(n_ndjfm,:,:)    ; NDJFM can replace by variable.
  v := v(n_ndjfm,:,:)
  reamp   := reamp(n_ndjfm)
  rephase := rephase(n_ndjfm)

  delete([/time,TIME,year,mon,day/])

  time = w&time
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

  m_value_1 = new((/2,4,nrandom,dimsizes(w&level),dimsizes(w&lat)/),float)
  m_value_1 = 0.
  print("random calc begin. "+ systemfunc("date"))
  do k = 0, nrandom - 1
    rd_time = monte_carlo1(ind_time)
    wt   = w(rd_time,:,:)
    vt   = v(rd_time,:,:)
    m_value_1(0,:,k,:,:) = calc_lagvalue(in2,phase,amp,wt)
    m_value_1(1,:,k,:,:) = calc_lagvalue(in2,phase,amp,vt)
    delete([/wt,vt,rd_time/])
  end do
  print("random calc done. "+ systemfunc("date"))

; sort values
  m_value1 = dim_pqsort_n(m_value_1, 2, 2)

  sigl = 0.05
  sigh = 0.95

  sigl1 = 0.1
  sigh1 = 0.9

  m_value = new((/2,4,4,dimsizes(w&level),dimsizes(w&lat)/),float) 
  m_value = 0.
  copy_VarCoords(w(0,:,:), m_value(0,0,0,:,:))

  nl  = toint(dimsizes(m_value_1(0,0,:,0,0))*sigl) - 1
  nh  = toint(dimsizes(m_value_1(0,0,:,0,0))*sigh) - 1
  m_value(:,:,0,:,:) = m_value_1(:,:,nl,:,:)
  m_value(:,:,1,:,:) = m_value_1(:,:,nh,:,:)

  nl1  = toint(dimsizes(m_value_1(:,0,0,0,0))*sigl1) - 1
  nh1  = toint(dimsizes(m_value_1(:,0,0,0,0))*sigh1) - 1
  m_value(:,:,2,:,:) = m_value_1(:,:,nl1,:,:)
  m_value(:,:,3,:,:) = m_value_1(:,:,nh1,:,:)

  system("/bin/rm -f "+outdir)
  ncdf = addfile(outdir,"c")
  ncdf->mvalue = m_value
  print("store done. "+ systemfunc("date"))