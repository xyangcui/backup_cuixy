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
  lats =  10
  latn =  80
  lonw = -120
  lone =  30

  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  fu   = addfile(ddir+"/ERA5_daily_high_"+"u"+"250_1979-2022.nc","r")
  fv   = addfile(ddir+"/ERA5_daily_high_"+"v"+"250_1979-2022.nc","r")

  u = fu->u(:,{lats:latn},{lonw:lone})
  v = fv->v(:,{lats:latn},{lonw:lone})

  time = u&time
  lat  = u&lat
  lon  = u&lon

  EKE = (u*u + v*v)*0.5
  copy_VarCoords(u,EKE)

  nwa  = 30.
  nwb  = 10.
  xa   = iso_filter(EKE,201,nwa,nwb,0)

  print("calc done. "+ systemfunc("date"))
  copy_VarCoords(EKE,xa)

  create_nc_file(xa,time,lat,lon,ddir+\
                 "/EKE250_daily_anom_1979-2022_r2.5.nc","EKE")

  create_nc_file(EKE,time,lat,lon,ddir+\
                 "/EKE250_daily_mean_1979-2022_r2.5.nc","EKE")