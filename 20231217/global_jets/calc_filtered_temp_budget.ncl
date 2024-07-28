;-----------------------------------------------------------------
;; this script is  used to filter values by removing seasonal
;; cycle. Lanczos filter.
;; only get budget term anomalies.
;; time: 2023-12-13
;-----------------------------------------------------------------
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_3D.ncl"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
begin
  nwa  = 100.
  nwb  = 20.

  ddir = "/home/sunming/data5/cuixy/global_jets/data/temp_budget"

  fin  = addfile(ddir+"/mean_"+var+".nc","r")

  time = fin->time
  lat  = fin->latitude
  lon  = fin->longitude

  xv   = fin->$var$
  
  xa   = iso_filter(xv,201,nwa,nwb,0)

  print("calc "+var+" done. "+ systemfunc("date"))
  copy_VarCoords(xv,xa)

  create_nc_file(xa,time,lat,lon,ddir+"/anom_"+var+".nc",var)

end

