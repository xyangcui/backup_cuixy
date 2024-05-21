;-----------------------------------------------------------------
;; this script is just used to filter values by removing seasonal
;; cycle. Lanczos filter.
;; time: 2023-11-22
;-----------------------------------------------------------------
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_3D.ncl"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
begin
  
  var  = "uwnd"

  nwa  = 90.
  nwb  = 10.

  ddir = "/home/sunming/data5/cuixy/ENSO-MJO-Jet/data"

  fin  = addfile(ddir+"/ERA5_u200_daily_mean_1979-2022_r1.nc","r")

  xv   = fin->$var$(:,{200},:,:)

  time = xv&time
  lat  = xv&lat
  lon  = xv&lon
  
  xa   = iso_filter(xv,201,nwa,nwb,1)

  print("calc done. "+ systemfunc("date"))
  copy_VarCoords(xv,xa)

  create_nc_file(xa,time,lat,lon,ddir+"/ERA5_"+var+"200_daily_anom_1979-2022_r1.nc",var)

end