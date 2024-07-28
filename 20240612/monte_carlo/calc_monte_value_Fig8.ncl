load "/home/sunming/data5/cuixy/script/ncl_func/monte_caro.ncl"
;--------------------------------------------------------
; calc u states. MEJS and EAJS.
;--------------------------------------------------------
  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  phidir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022_2.5x2.5.nc"
  udir   = ddir+"/"+"ERA5_daily_anom_u250_1979-2022_2.5x2.5.nc"
  vdir   = ddir+"/"+"ERA5_daily_anom_v250_1979-2022_2.5x2.5.nc"
  tdir   = ddir+"/"+"ERA5_daily_anom_tvert_1979-2022_2.5x2.5.nc"

  outdir = ddir+"/"+"monte_carlo/Fig8_usft_sig.nc"

  lats = 0
  latn = 60
  lonw = 200
  lone = 350
;--------------------------------------------------------
; function calc_lagvalue.
;--------------------------------------------------------
function calc_lagvalue(in22,phaset,ampt,value)
local in22,phaset,ampt,value,DIMS,dim1,dim2,llvalue,\
nt,ntt,nt1,nindd
begin
  DIMS = dimsizes(value)
  dim1 = DIMS(1)
  dim2 = DIMS(2)
  llvalue = new((/dim1,dim2/),float)
  llvalue = 0.

  nt = ind(phaset.eq.4.and.ampt.gt.1.)
  nind = in22(nt)
  ntt  = 1
  nindd = nind + ntt
  do ntt = 2,10
    nit = nind + ntt
    nindd := array_append_record(nindd, nit, 0)
    delete(nit)
  end do
  llvalue = dim_avg_n(value(nindd,:,:), 0)
  delete([/nt,ntt,nindd,nind/])

  return llvalue

end
;--------------------------------------------------------
; function calc_lagvalue Fig8.
;--------------------------------------------------------
function calc_lagvaluet(in22,phaset,ampt,value)
local in22,phaset,ampt,value,DIMS,dim1,dim2,llvalue,\
nt,ntt,nt1,nindd
begin
  DIMS = dimsizes(value)
  dim1 = DIMS(1)
  dim2 = DIMS(2)
  llvalue = new((/dim1,dim2/),float)
  llvalue = 0.

  nt = ind(phaset.eq.8.and.ampt.gt.1.)
  nind = in22(nt)
  ntt  = 6
  nindd = nind + ntt
  do ntt = 7,15
    nit = nind + ntt
    nindd := array_append_record(nindd, nit, 0)
    delete(nit)
  end do
  llvalue = dim_avg_n(value(nindd,:,:), 0)
  delete([/nt,ntt,nindd,nind/])

  return llvalue

end
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  fphi = addfile(phidir,"r")
  fu   = addfile(udir,"r")
  fv   = addfile(vdir,"r")
  ft = addfile(tdir,"r")

  t    = ft->t(:,{lats:latn},{lonw:lone})
  phi  = fphi->phi(:,{lats:latn},{lonw:lone})
  u    = fu->u(:,{lats:latn},{lonw:lone})
  v    = fv->v(:,{lats:latn},{lonw:lone})

  time = phi&time
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
  u   := u(n1:n2,:,:)
  v   := v(n1:n2,:,:)
  t   := t(n1:n2,:,:)
  phi := phi(n1:n2,:,:)
  reamp   := reamp(n1:n2)
  rephase := rephase(n1:n2)
  n_ndjfm = ind(mont.eq.11.or.mont.eq.12.or.mont.eq.1.or.mont.eq.2.or.mont.eq.3)
  u   := u(n_ndjfm,:,:)    ; NDJFM can replace by variable.
  v   := v(n_ndjfm,:,:)
  t   := t(n_ndjfm,:,:)
  phi := phi(n_ndjfm,:,:)
  reamp   := reamp(n_ndjfm)
  rephase := rephase(n_ndjfm)

  delete([/time,TIME,year,mon,day/])

  time = u&time
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

  m_value_2 = new((/4,2,nrandom,dimsizes(u&lat),dimsizes(u&lon)/),float)
  m_value_2 = 0.
  print("random calc begin. "+ systemfunc("date"))

  do k = 0, nrandom - 1
    rd_time = monte_carlo1(ind_time)
    ut = u(rd_time,:,:)
    vt = v(rd_time,:,:)
    tt = t(rd_time,:,:)
    phit = phi(rd_time,:,:)
    m_value_2(0,0,k,:,:) = calc_lagvalue(in2,phase,amp,ut)
    m_value_2(0,1,k,:,:) = calc_lagvaluet(in2,phase,amp,ut)

    m_value_2(1,0,k,:,:) = calc_lagvalue(in2,phase,amp,vt)
    m_value_2(1,1,k,:,:) = calc_lagvaluet(in2,phase,amp,vt)

    m_value_2(2,0,k,:,:) = calc_lagvalue(in2,phase,amp,phit)
    m_value_2(2,1,k,:,:) = calc_lagvaluet(in2,phase,amp,phit)

    m_value_2(3,0,k,:,:) = calc_lagvalue(in2,phase,amp,tt)
    m_value_2(3,1,k,:,:) = calc_lagvaluet(in2,phase,amp,tt)
    delete([/ut,vt,tt,phit,rd_time/])
  end do
  print("random calc done. "+ systemfunc("date"))
; store values
  m_value2 = dim_pqsort_n(m_value_2, 2, 2)

  mvalue2  = new((/4,2,4,dimsizes(u&lat),dimsizes(u&lon)/),float) 
  mvalue2  = 0.
  copy_VarCoords(u(0,:,:), mvalue2(0,0,0,:,:))
  mvalue2@long_name = "p4/P8"
  mvalue2@store = "(u),0.05,0.95,0.1,0.9"

  sigl = 0.05
  sigh = 0.95
  nl2  = toint(dimsizes(m_value_2(0,0,:,0,0))*sigl) - 1
  nh2  = toint(dimsizes(m_value_2(0,0,:,0,0))*sigh) - 1
  mvalue2(:,:,0,:,:) = m_value_2(:,:,nl2,:,:)
  mvalue2(:,:,1,:,:) = m_value_2(:,:,nh2,:,:)

  sigl = 0.1
  sigh = 0.9
  nl2  = toint(dimsizes(m_value_2(0,0,:,0,0))*sigl) - 1
  nh2  = toint(dimsizes(m_value_2(0,0,:,0,0))*sigh) - 1
  mvalue2(:,:,2,:,:) = m_value_2(:,:,nl2,:,:)
  mvalue2(:,:,3,:,:) = m_value_2(:,:,nh2,:,:)

  system("/bin/rm -f "+outdir)
  ncdf = addfile(outdir,"c")
  ncdf->mvalue = mvalue2
  print("store done. "+ systemfunc("date"))