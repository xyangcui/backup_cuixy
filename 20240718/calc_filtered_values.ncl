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
  
  var  = "u"

  nwa  = 30.
  nwb  = 10.

  ddir = "/home/sunming/data5/cuixy/global_jets/data"

  fin  = addfile(ddir+"/ERA5_daily_u250_1979-2022_r2.5.nc","r")

  xv   = fin->$var$(:,{250},:,:)

  time = xv&time
  lon  = xv&lon
  lat  = xv&lat

  printVarSummary(xv)

  xa   = iso_filter(xv,201,nwa,nwb,0)

  print("calc done. "+ systemfunc("date"))
  copy_VarCoords(xv,xa)

  create_nc_file(xa,time,lat,lon,ddir+\
                 "/"+var+"250_daily_anom_1979-2022_10-100_r2.5.nc",var)
end