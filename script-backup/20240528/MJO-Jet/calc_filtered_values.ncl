;-----------------------------------------------------------------
;; this script is just used to filter values by removing seasonal
;; cycle. Lanczos filter.
;; time: 2023-11-22
;-----------------------------------------------------------------
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_4D.ncl"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
begin
  
  var  = "w"

  nwa  = 100.
  nwb  = 20.

  ;lats =  0
  ;latn =  90
  ;lonw =  90
  ;lone =  200

  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJS_temp_budget"

  fin  = addfile(ddir+"/ERA5_daily_wlevel_1979-2022_r2.5.nc","r")

  xv   = fin->$var$

  printVarSummary(xv)

  time  = xv&time
  level = xv&level
  lat   = xv&lat
  lon   = xv&lon

  ymdhms  = cd_calendar(time, 0)
  iyear   = floattoint( ymdhms(:,0) )
  imon    = floattoint( ymdhms(:,1) )   
  iday    = floattoint( ymdhms(:,2) )
  ddd     = day_of_year(iyear, imon, iday) 
  yyyyddd = iyear*1000 + ddd 
  wts     = filwgts_lanczos (201, 2, 1./nwa,1./nwb , 1.0) 

  print("calc begin. "+ systemfunc("date"))

  xClim = clmDayTLLL(xv, yyyyddd)

  xAnom = xv
  xAnom = 0.
  do i = 0,8
    xAnom(:,i,:,:) = calcDayAnomTLL(xv(:,i,:,:), yyyyddd, xClim(:,i,:,:))
  end do

  xa = wgt_runave_n_Wrap(xAnom, wts, 0, 0)

  print("calc done. "+ systemfunc("date"))
  copy_VarCoords(xv,xa)

  create_nc_file(xa,time,level,lat,lon,ddir+"/ERA5_daily_anom_"+var+"level_1979-2022_r2.5.nc",var)

end