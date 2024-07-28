load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_4D.ncl"
;-----------------------------------------------------------------
;; this script is just used to filter values by removing seasonal
;; cycle. Lanczos filter.
;; Logicical: read in single level values.
;; calculate filtered values.
;; selected one year value into a martirx.
;; calc all level and form the martrix and save as a file.
;; time: 2023-11-29
;-----------------------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  idir = ddir + "/ERA5_daily_mean_vlevel_1979-2022_r2.5.nc"
  odir = ddir + "/ERA5_daily_anom_vlevel_1979-2022_r2.5.nc"

  lats =  0
  latn =  90
  lonw = -90
  lone =  40

  nwa  = 30.
  nwb  = 10.

  var  = "v"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
  fin = addfile(idir,"r")

  x = fin->$var$(:,:,{lats:latn},{lonw:lone})
  time  = x&time
  level = x&level
  lat   = x&lat
  lon   = x&lon

  xp = new((/dimsizes(time),dimsizes(level),\
          dimsizes(lat),dimsizes(lon)/),float)
  xp = 0.0

  do i = 0, dimsizes(level)-1
    print("level " + i + " begin. "+systemfunc("date"))
    xt = x(:,i,:,:)
    x&time = time
    xa = iso_filter(xt,201,nwa,nwb,0)
    xp(:,i,:,:) = xa
    delete([/xa,xt/])
    print("level " + i + " done. "+systemfunc("date"))
  end do

  copy_VarCoords(x, xp)
    
  create_nc_file(xp,time,level,lat,lon,odir,var)