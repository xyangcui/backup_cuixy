load "/home/sunming/data5/cuixy/script/ncl_func/monte_caro.ncl"
;---------------------------------------------------------------
; draw temperature anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour.
; time: 2023-12-17
;---------------------------------------------------------------
  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  vdir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022_2.5x2.5.nc"
  outdir = ddir+"/"+"monte_carlo/phiu_sig_phase15.nc"
  latmin = -30
  latmax = 90
;-------------------------------------------------------------
; function calc_lagvalue.
; note that this function should change follow your intention.
;-------------------------------------------------------------
function calc_lagvalue(in22,phaset,ampt,value)
local in22,phaset,ampt,value,DIMS,dim1,dim2,llvalue,\
ip,nt,ntt,nt1,nindd
begin
  DIMS = dimsizes(value)
  dim1 = DIMS(1)
  dim2 = DIMS(2)
  llvalue = new((/2,8,dim1,dim2/),float)
  llvalue = 0.
  iphase  = (/1,5/)
  do i = 0, 1
    ip = iphase(i)
    nt = ind(phaset.eq.ip.and.ampt.gt.1.)
    nt1 = in22(nt)
    do j = 0, 7
      ntt   = -15 + j*5
      nindd = nt1 + ntt
      llvalue(i,j,:,:) = dim_avg_n(value(nindd,:,:), 0)
      delete([/ntt,nindd/])
    end do
    delete([/ip,nt,nt1/])
  end do

  return llvalue
end
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  f  = addfile(vdir,"r")

  x  = f->$var$(:,{latmin:latmax},:)
  lat = x&lat
  lon = x&lon
  print("complete data selling "+ systemfunc("date"))

  time = x&time
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
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i)\
             .and.iday.eq.day(i))
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
  x   := x(n1:n2,:,:)
  reamp   := reamp(n1:n2)
  rephase := rephase(n1:n2)
  n_ndjfm = ind(mont.eq.11.or.mont.eq.12.or.mont.eq.1\
                .or.mont.eq.2.or.mont.eq.3)
  x       := x(n_ndjfm,:,:)  ; NDJFM can replace by variable.
  reamp   := reamp(n_ndjfm)
  rephase := rephase(n_ndjfm)

  delete([/time,TIME,year,mon,day/])

  time = x&time
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
  m_value = new((/nrandom,1,8,dimsizes(x&lat),dimsizes(x&lon)/),float)
  m_value = 0.

  print("random calc phi begin. "+ systemfunc("date"))
  do k = 0, nrandom - 1
    rd_time = monte_carlo1(ind_time)
    xt = x(rd_time,:,:)
    m_value(k,:,:,:,:) = calc_lagvalue(in2,phase,amp,xt)
    delete([/xt,rd_time/])
  end do
  copy_VarCoords(x(0,:,:), m_value(0,0,0,:,:))

  rdtimes = ispan(1, 10000, 1)
  m_value!0 = "time"
  m_value&time = rdtimes
  print("random calc done. "+ systemfunc("date"))

; store values.
  m_valuet = dim_pqsort_n(m_value, 2, 0)

  system("/bin/rm -f "+outdir)
  fout = addfile(outdir,"c")

  dimNames = (/"time", "lat", "lon"/)  
  dimSizes = (/ -1   ,  dimsizes(lat),  dimsizes(lon)/) 
  dimUnlim = (/ True , False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout, "time" ,typeof(rdtimes),getvardims(rdtimes))                           
  filevardef(fout, "lat"  ,typeof(lat),getvardims(lat))                          
  filevardef(fout, "lon"  ,typeof(lon),getvardims(lon))                          
  filevardef(fout, VAR    ,typeof(x)  ,getvardims(x))    
                                                                                             
  filevarattdef(fout,"time" ,rdtimes)                    
  filevarattdef(fout,"lat"  ,lat)                   
  filevarattdef(fout,"lon"  ,lon)                     
  filevarattdef(fout,VAR   ,x)   

  setfileoption(fout,"DefineMode",False)

  fout->time   = (/rdtimes/)     
  fout->lat    = (/lat/)
  fout->lon    = (/lon/) 
  fout->$VAR$  = (/x/)

  print("store done. "+ systemfunc("date"))