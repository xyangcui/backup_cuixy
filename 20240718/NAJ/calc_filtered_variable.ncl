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

  var = "v"

  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  fin  = addfile(ddir+"/ERA5_daily_mean_"+var+"level_1979-2022_r2.5.nc","r")

  x   = fin->$var$(:,{250},{lats:latn},{lonw:lone})

  time = x&time
  lat  = x&lat
  lon  = x&lon

  wgth = filwgts_lanczos(101, 1, 1./8., -999., 1.)
  xh   = wgt_runave_n_Wrap(x, wgth, 0, 0)
  copy_VarCoords(x, xh)
  create_nc_file(xh,time,lat,lon,\
                   ddir+"/ERA5_daily_"+"high"+"_"+var+"250_1979-2022.nc",var)
  ;system("/bin/rm -f "+ddir+"/ERA5_daily_"+"high"+"_"+var+"250_1979-2022.nc")
  ;fout = addfile(ddir+"/ERA5_daily_"+"high"+"_"+var+"250_1979-2022.nc","c")
  ;fout->$var$ = xh
  ;delete([/fout,xh/])
  ;delete(xh)
  print("calc high done. "+ systemfunc("date"))