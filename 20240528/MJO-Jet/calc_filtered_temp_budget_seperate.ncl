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

  fultl  = addfile(ddir+"/mean_"+"ultl"+".nc","r")
  fulti  = addfile(ddir+"/mean_"+"ulti"+".nc","r")
  fulth  = addfile(ddir+"/mean_"+"ulth"+".nc","r")
  fuitl  = addfile(ddir+"/mean_"+"uitl"+".nc","r")
  fuiti  = addfile(ddir+"/mean_"+"uiti"+".nc","r")
  fuith  = addfile(ddir+"/mean_"+"uith"+".nc","r")
  fuhtl  = addfile(ddir+"/mean_"+"uhtl"+".nc","r")
  fuhti  = addfile(ddir+"/mean_"+"uhti"+".nc","r")
  fuhth  = addfile(ddir+"/mean_"+"uhth"+".nc","r")

  time = fultl->time
  lat  = fultl->latitude
  lon  = fultl->longitude

;ultl
  var = "ultl"
  x  = fultl->udtdx
  xa = iso_filter(x,201,nwa,nwb,0)
  copy_VarCoords(x,xa)
  create_nc_file(xa,time,lat,lon,ddir+"/anom_"+var+".nc",var)
  delete([/x,xa/])
  
;ulti
  var = "ulti"
  x  = fulti->udtdx
  xa = iso_filter(x,201,nwa,nwb,0)
  copy_VarCoords(x,xa)
  create_nc_file(xa,time,lat,lon,ddir+"/anom_"+var+".nc",var)
  delete([/x,xa/])

;ulth
  var = "ulth"
  x  = fulth->udtdx
  xa = iso_filter(x,201,nwa,nwb,0)
  copy_VarCoords(x,xa)
  create_nc_file(xa,time,lat,lon,ddir+"/anom_"+var+".nc",var)
  delete([/x,xa/])
  
;uitl
  var = "uitl"
  x  = fuitl->udtdx
  xa = iso_filter(x,201,nwa,nwb,0)
  copy_VarCoords(x,xa)
  create_nc_file(xa,time,lat,lon,ddir+"/anom_"+var+".nc",var)
  delete([/x,xa/])
;uiti
  var = "uiti"
  x  = fuiti->udtdx
  xa = iso_filter(x,201,nwa,nwb,0)
  copy_VarCoords(x,xa)
  create_nc_file(xa,time,lat,lon,ddir+"/anom_"+var+".nc",var)
  delete([/x,xa/])
;uith
  var = "uith"
  x  = fuith->udtdx
  xa = iso_filter(x,201,nwa,nwb,0)
  copy_VarCoords(x,xa)
  create_nc_file(xa,time,lat,lon,ddir+"/anom_"+var+".nc",var)
  delete([/x,xa/])
  
;uhtl
  var = "uhtl"
  x  = fuhtl->udtdx
  xa = iso_filter(x,201,nwa,nwb,0)
  copy_VarCoords(x,xa)
  create_nc_file(xa,time,lat,lon,ddir+"/anom_"+var+".nc",var)
  delete([/x,xa/])

;uhti
  var = "uhti"
  x  = fuhti->udtdx
  xa = iso_filter(x,201,nwa,nwb,0)
  copy_VarCoords(x,xa)
  create_nc_file(xa,time,lat,lon,ddir+"/anom_"+var+".nc",var)
  delete([/x,xa/])

;uhth
  var = "uhth"
  x  = fuhth->udtdx
  xa = iso_filter(x,201,nwa,nwb,0)
  copy_VarCoords(x,xa)
  create_nc_file(xa,time,lat,lon,ddir+"/anom_"+var+".nc",var)
  delete([/x,xa/])

  print("calc "+var+" done. "+ systemfunc("date"))

end

