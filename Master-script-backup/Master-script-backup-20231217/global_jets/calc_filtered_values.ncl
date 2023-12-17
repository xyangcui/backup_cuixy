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
  var  = "z"
  nwa  = 100.
  nwb  = 20.

  ddir = "/home/sunming/data5/cuixy/global_jets/data"

  fin  = addfile(ddir+"/ERA5_daily_"+var+"250_1979-2022.nc","r")

  time = fin->time
  lat  = fin->latitude
  lon  = fin->longitude

  xv   = fin->$var$(:,0,:,:)
  
  xa   = iso_filter(xv,201,nwa,nwb,0)

  print("calc done. "+ systemfunc("date"))
  copy_VarCoords(xv,xa)

  create_nc_file(xa,time,lat,lon,ddir+"/ERA5_daily_anom_"+var+"250_1979-2022.nc",var)

end

