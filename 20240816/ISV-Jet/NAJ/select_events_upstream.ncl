;-----------------------------------------------------------------
;; this script is just used to filter values by removing seasonal
;; cycle. Lanczos filter.
;; time: 2023-12-18.
;-----------------------------------------------------------------
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_3D.ncl"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
  lats =  28
  latn =  45
  lonw = -90
  lone = -60

  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  fu   = addfile(ddir+"/ERA5_daily_anom_ulevel_1979-2022_r2.5.nc","r")

  u   = fu->u(:,{250},{lats:latn},{lonw:lone})

  time = u&time
  lat  = u&lat
  lon  = u&lon

  rad = 4.0*atan(1.0)/180.0
  clat = cos(lat * rad)
  u := wgt_areaave_Wrap(u, clat, 1.0, 1)
;; calc mht index.
; time table.
  idir = ddir + "/EKE250_daily_anom_1979-2022_r2.5.nc"
  f = addfile(idir,"r")
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.4.and.day.eq.1)

  mon:= mon(n1:n2)
  u  := u(n1:n2)
  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  ut = u(nt)

  ji = dim_standardize_n_Wrap(ut, 0, 0)
;;
  mdir = ddir+"/meri_heatflux.csv"
  delim  = ","
  mf  = asciiread(mdir,-1,"string")
  mi  = str_get_field(mf,1,delim)
  mi := tofloat(mi)
;; calc corr.
  mxlag    = 10
  x_Lead_y = esccr(ji,mi,mxlag)
  y_Lead_x = esccr(mi,ji,mxlag)    

  ccr = new ( 2*mxlag+1, float)    
  ccr(0:mxlag-1) = x_Lead_y(1:mxlag:-1)  
  ccr(mxlag:)    = y_Lead_x(0:mxlag) 

  print(ccr)  

  mxlag = 2
  mac = esacr(mi,mxlag)
  vm = -1*dimsizes(mi)*0.5*log(mac(1))

  jac = esacr(ji,mxlag)
  vj = -1*dimsizes(ji)*0.5*log(jac(1))
  vv = (vj+vm)/2

  print(vv)  