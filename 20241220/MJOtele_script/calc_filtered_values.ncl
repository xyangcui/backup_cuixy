;-----------------------------------------------------------------
;; this script is just used to filter values by removing seasonal
;; cycle. Lanczos filter.
;; time: 2024-05-13
;-----------------------------------------------------------------
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_3D.ncl"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
begin
  
  var  = "v"

  nwa  = 100.
  nwb  = 20.

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  fin  = addfile(ddir+"/ERA5_v200_daily_mean_1979-2022_r2.5.nc","r")

  xv = fin->$var$(:,{200},:,:)
  x  = xv;dtrend_n(xv, False, 0)
  time = xv&time
  lat  = xv&lat
  lon  = xv&lon
  copy_VarCoords(xv,x)
  xa   = iso_filter(x,201,nwa,nwb,1)
  print("calc done. "+ systemfunc("date"))
  copy_VarCoords(xv,xa)

  create_nc_file(xa,time,lat,lon,ddir+\
                "/"+var+"200_daily_anom_1979-2022_r2.5.nc",var)
end