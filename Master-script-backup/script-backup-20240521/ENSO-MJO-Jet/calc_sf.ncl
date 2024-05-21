load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_3D.ncl"
;--------------------------------------------------------
; calculate streamfunction by sperical function.
; Note that u and v must global data.
; author: cui xiangyang time: 2023-11-21
;--------------------------------------------------------

  ddir = "/home/sunming/data5/cuixy/ENSO-MJO-Jet/data"
  udir = ddir+"/"+"u250_daily_mean_1979-2022_r2.5.nc"
  vdir = ddir+"/"+"v250_daily_mean_1979-2022_r2.5.nc"
  sdir = ddir+"/"+"sf250_daily_mean_1979-2022_r2.5.nc"

begin

  fu = addfile(udir,"r")
  fv = addfile(vdir,"r")

  lt = fu->lat

  if(lt(0).lt.lt(1)) then
    u  = fu->u(:,{250},:,:)
    v  = fv->v(:,{250},:,:)
  else
    u  = fu->u(:,{250},::-1,:)
    v  = fv->v(:,{250},::-1,:)
  end if

  printVarSummary(u)
  delete(lt)

  time = u&time
  lat  = u&lat
  lon  = u&lon

  print("calc begin. "+ systemfunc("date"))

  sfvp = uv2sfvpF(u,v)
  delete(v)

  sf   = sfvp(0,:,:,:)
  copy_VarCoords(u, sf)

  print("calc done. "+ systemfunc("date"))

  sf@unit = "m2/s"
  sf@long_name = "streamfunction at 250hPa"
  sf@standard_name = "streamfunction"
  delete([/sfvp,u/])
  
  create_nc_file(sf,time,lat,lon,sdir,"sf")
   
end
