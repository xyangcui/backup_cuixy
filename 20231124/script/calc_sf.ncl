;--------------------------------------------------------
; calculate streamfunction by sperical function.
; Note that u and v must global data.
; author: cui xiangyang time: 2023-11-21
;--------------------------------------------------------
begin
  ddir = "/home/sunming/data5/cuixy/global_jets/data"

  fu = addfile(ddir+"/"+"ERA5_daily_u250_1979-2022_all.nc","r")
  fv = addfile(ddir+"/"+"ERA5_daily_v250_1979-2022_all.nc","r")

  u  = fu->u(:,0,::-1,:)
  v  = fv->v(:,0,::-1,:)
  time = u&time
  lat  = u&latitude
  lon  = u&longitude
  sfvp = uv2sfvpF(u,v)
  delete([/v/])
  sf   = sfvp(0,:,:,:)
  copy_VarCoords(u, sf)
  sf@unit = "m2/s"
  sf@long_name = "streamfunction at 250hPa"
  sf@standard_name = "Streamfunction"
  delete([/sfvp,u/])
  
  system("/bin/rm -f "+ddir+"/"+"ERA5_daily_sf250_1979-2022_all1.nc")

  fout = addfile(ddir+"/"+"ERA5_daily_sf250_1979-2022_all1.nc","c")
   
  ;setfileoption(fout,"DefineMode",True)

  ;dimNames = (/"time","latitude","longitude"/)
  ;dimSizes = (/-1, dimsizes(lat), dimsizes(lon)/)
  ;dimUlim  = (/True, False, False/)
  ;filedimdef(fout, dimNames, dimSizes, dimUlim)

  ;filevardef(fout, "time", typeof(time), getvardims(time))
  ;filevardef(fout, "latitude", typeof(lat), getvardims(lat))
  ;filevardef(fout, "longitude", typeof(lon), getvardims(lon))
  ;filevardef(fout, "phi", typeof(sf), getvardims(sf))

  ;filevarattdef(fout,"time",time)
  ;filevarattdef(fout,"latitude",lat)
  ;filevarattdef(fout,"longitude",lon)
  ;filevarattdef(fout,"phi",sf)
 
  ;setfileoption(fout,"DefineMode",False)
  
  ;fout->time      = (/time/)
  ;fout->latitude  = (/lat/)
  ;fout->longitude = (/lon/)
  ;fout->phi       = (/sf/)
  fout->phi = sf
end
