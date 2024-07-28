;------------------------------------------
;; calculate correlation between Jet index
;;   and filtered daily NAO index.
;; daily NAO index contains 10-30 days.
;; Keely et al. 2009 Bendict et al. 2004.
;; date: 20240712
;------------------------------------------

;------------------------------------------
; Function to filter.
;------------------------------------------
undef("ISO_Filter")
function ISO_Filter(x,time,npt,nwa,nwb)
local time
  begin
    TIME1 = cd_calendar(time, 0)
    year1 = toint(TIME1(:,0))
    mon1  = toint(TIME1(:,1))
    day1  = toint(TIME1(:,2))
    ddd   = day_of_year(year1, mon1, day1)
    yyyyddd = year1*1000 + ddd
    wts     = filwgts_lanczos (npt, 2, 1./nwa,1./nwb , 1.0)
    xtmp  = new((/dimsizes(x),2,2/),float)
    xtmp  = 0.0
    xtmp(:,0,0) = x
    xClim = clmDayTLL(xtmp, yyyyddd)
    xAnom = calcDayAnomTLL(xtmp, yyyyddd, xClim)
    xr    = xAnom(:,0,0)
    xr := wgt_runave_Wrap(xr, wts, 0)
  
    return(xr)
    delete([/TIME1,year1,mon1,day1,xtmp,xAnom,xr/])
  
  end

ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
jdir = ddir+"/Ulevel_PC.csv";"/meri_heatflux.csv";
mdir = ddir+"/meri_heatflux.csv";"/eke_index.csv";
ndir = "/home/sunming/data5/cuixy/data/clim_indices/NAO-daily-index.csv"
;-------------------------------------------
; Main code.
;-------------------------------------------
; read in jet index. DJF.
  delim  = ","
  jf  = asciiread(jdir,-1,"string")
  ji  = str_get_field(jf,1,delim)
  ji := tofloat(ji)

  mf  = asciiread(mdir,-1,"string")
  mi  = str_get_field(mf,1,delim)
  mi := tofloat(mi)*(-1.)

; time table.
  idir = ddir + "/ERA5_daily_anom_ulevel_1979-2022_r2.5.nc"
  f = addfile(idir,"r")
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

; read in ndir.
  nf    = asciiread(ndir,-1,"string")
  iyear = str_get_field(nf,1,delim)
  imon  = str_get_field(nf,2,delim)
  iday  = str_get_field(nf,3,delim)
  njt   = str_get_field(nf,4,delim)

  iyear := toint(iyear(1::))
  imon  := toint(imon(1::))
  iday  := toint(iday(1::))
  njt   := tofloat(njt(1::))

  nj = new((/dimsizes(time)/),float)
  nj = 0.

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
    nj(i) = njt(nt)
    delete(nt)
  end do

  nn = ispan(0, dimsizes(time)-1, 1)
  delete([/nn,njt,iyear,imon,iday/])

;; calc corr.
  mxlag    = 10
  x_Lead_y = esccr(mi,ji,mxlag)
  y_Lead_x = esccr(ji,mi,mxlag)    ; switch the order of the series

  ccr = new ( 2*mxlag+1, float)    
  ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
  ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"

  print(ccr)