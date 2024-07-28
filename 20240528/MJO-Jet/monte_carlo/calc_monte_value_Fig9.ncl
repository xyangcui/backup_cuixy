load "/home/sunming/data5/cuixy/script/ncl_func/monte_caro.ncl"
;--------------------------------------------------------
; calc sf and u and v states. MEJS and EAJS.
;--------------------------------------------------------
  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
  
  phidir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022_2.5x2.5.nc"
  outdir = ddir+"/"+"monte_carlo/Fig9_sf_sig.nc"

  lats = -20
  latn =  90
  lonw = 100
  lone = 360
;--------------------------------------------------------
; function calc_lagvalue.
;--------------------------------------------------------
function calc_lagvalue(in22,phaset,ampt,value,ip)
local in22,phaset,ampt,value,DIMS,dim1,dim2,llvalue,\
ip,nt,ntt,nt1,nindd
begin
  DIMS = dimsizes(value)
  dim1 = DIMS(1)
  dim2 = DIMS(2)
  llv = new((/4,dim1,dim2/),float)
  llv = 0.

  nt = ind(phaset.eq.ip.and.ampt.gt.1.)
  nind = in22(nt)
  do j = 0,3
    ntt = -10 + j*5
    nit = nind + ntt
    llv(j,:,:) = dim_avg_n(value(nit,:,:), 0)
    delete([/ntt,nit/])
  end do

  return llv
end
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------
  fphi = addfile(phidir,"r")

  phi  = fphi->phi(:,{lats:latn},{lonw:lone})

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
  phi := phi(n1:n2,:,:)
  reamp   := reamp(n1:n2)
  rephase := rephase(n1:n2)
  n_ndjfm = ind(mont.eq.11.or.mont.eq.12.or.mont.eq.1.or.mont.eq.2.or.mont.eq.3)
  phi := phi(n_ndjfm,:,:)  
  reamp   := reamp(n_ndjfm)
  rephase := rephase(n_ndjfm)

  delete([/time,TIME,year,mon,day/])

  time = phi&time
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

  m_value_2 = new((/2,nrandom,4,dimsizes(phi&lat),dimsizes(phi&lon)/),float)
  m_value_2 = 0.
  print("random calc u begin. "+ systemfunc("date"))

  do k = 0, nrandom - 1
    rd_time = monte_carlo1(ind_time)
    phit = phi(rd_time,:,:)

    m_value_2(0,k,:,:,:) = calc_lagvalue(in2,phase,amp,phit,4)
    m_value_2(1,k,:,:,:) = calc_lagvalue(in2,phase,amp,phit,8)
    delete([/phit,rd_time/])
  end do
  print("random calc done. "+ systemfunc("date"))
; store values
  m_value2 = dim_pqsort_n(m_value_2, 2, 1)

  mvalue2  = new((/2,4,4,dimsizes(phi&lat),dimsizes(phi&lon)/),float) 
  mvalue2  = 0.
  copy_VarCoords(phi(0,:,:), mvalue2(0,0,0,:,:))
  mvalue2@long_name = "p2/P6"
  mvalue2@store = "(sf,u,v),0.05,0.95,0.1,0.9"

  sigl = 0.05
  sigh = 0.95
  nl2  = toint(dimsizes(m_value_2(0,:,0,0,0))*sigl) - 1
  nh2  = toint(dimsizes(m_value_2(0,:,0,0,0))*sigh) - 1
  mvalue2(:,0,:,:,:) = m_value_2(:,nl2,:,:,:)
  mvalue2(:,1,:,:,:) = m_value_2(:,nh2,:,:,:)

  sigl = 0.1
  sigh = 0.9
  nl2  = toint(dimsizes(m_value_2(0,:,0,0,0))*sigl) - 1
  nh2  = toint(dimsizes(m_value_2(0,:,0,0,0))*sigh) - 1
  mvalue2(:,2,:,:,:) = m_value_2(:,nl2,:,:,:)
  mvalue2(:,3,:,:,:) = m_value_2(:,nh2,:,:,:)

  system("/bin/rm -f "+outdir)
  ncdf = addfile(outdir,"c")
  ncdf->mvalue = mvalue2
  print("store done. "+ systemfunc("date"))