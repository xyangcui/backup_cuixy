;-----------------------------------------------------------------
;; this script is just used to filter values by removing seasonal
;; cycle. Lanczos filter.
;; time: 2023-12-18.
;-----------------------------------------------------------------
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_4D.ncl"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
begin
  latmin = 20-1
  latmax = 40+1
  lonmin = 10-1
  lonmax = 180+1
  ddir = "/home/sunming/data5/cuixy/global_jets/data"
  ;tdir = ddir + "/ERA5_daily_mean_tlevel"
  ;files  = systemfunc("ls " + tdir + "/day_t_*.nc")
  ;fin    = addfiles(files,"r")
  fin   = addfile(ddir+"/ERA5_daily_mean_"+var+"level_1979-2022.nc","r")
  ;time  = fin1->time
  x    = fin->$var$(:,{250:500},{latmin:latmax},{lonmin:lonmax})
  ;x&time = time
  ;delete(time)
  time   = x&time
  lat    = x&latitude
  lon    = x&longitude
  level  = x&level

  wgtl = filwgts_lanczos(101, 0, 1./100., -999., 1.)
  xl   = wgt_runave_n_Wrap(x, wgtl, 0, 0)
  copy_VarCoords(x, xl)
  create_nc_file(xl,time,level,lat,lon,\
                 ddir+"/temp_budget/ERA5_daily_"+"low"+"_"+var+"250_1979-2022.nc",var)
  ;system("/bin/rm -f "+ddir+"/temp_budget/ERA5_daily_"+"low"+"_"+var+"250_1979-2022.nc")
  ;fout = addfile(ddir+"/temp_budget/ERA5_daily_"+"low"+"_"+var+"250_1979-2022.nc","c")
  ;fout->$var$ = xl
  ;delete([/fout,xl/])
  delete(xl)
  print("calc low done. "+ systemfunc("date"))

  wgti = filwgts_lanczos(101, 2, 1./100., 1./20., 1.)
  xi   = wgt_runave_n_Wrap(x, wgti, 0, 0)
  copy_VarCoords(x, xi)
  create_nc_file(xi,time,level,lat,lon,\
                 ddir+"/temp_budget/ERA5_daily_"+"iso"+"_"+var+"250_1979-2022.nc",var)
  ;system("/bin/rm -f "+ddir+"/temp_budget/ERA5_daily_"+"iso"+"_"+var+"250_1979-2022.nc")
  ;fout = addfile(ddir+"/temp_budget/ERA5_daily_"+"iso"+"_"+var+"250_1979-2022.nc","c")
  ;fout->$var$ = xi
  ;delete([/fout,xi/])
  delete(xi)
  print("calc iso done. "+ systemfunc("date"))

  wgth = filwgts_lanczos(101, 1, 1./20., -999., 1.)
  xh   = wgt_runave_n_Wrap(x, wgth, 0, 0)
  copy_VarCoords(x, xh)
  create_nc_file(xh,time,level,lat,lon,\
                   ddir+"/temp_budget/ERA5_daily_"+"high"+"_"+var+"250_1979-2022.nc",var)
  ;system("/bin/rm -f "+ddir+"/temp_budget/ERA5_daily_"+"high"+"_"+var+"250_1979-2022.nc")
  ;fout = addfile(ddir+"/temp_budget/ERA5_daily_"+"high"+"_"+var+"250_1979-2022.nc","c")
  ;fout->$var$ = xh
  ;delete([/fout,xh/])
  delete(xh)
  print("calc high done. "+ systemfunc("date"))

end

