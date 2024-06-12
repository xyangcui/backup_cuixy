
;----------------------------------------------------
;----------------------------------------------------

  ddir = "/home/sunming/data5/cuixy/global_jets/data"
  udir = ddir+"/ERA5_daily_u250_1979-2022_r2.5.nc"
  vdir = ddir+"/ERA5_daily_v250_1979-2022_r2.5.nc"
  odir = ddir+"/ERA5_daily_vor_1979-2022_r2.5.nc"
  dvudir = ddir+"/ERA5_daily_divu_1979-2022_r2.5.nc"
  dvvdir = ddir+"/ERA5_daily_divv_1979-2022_r2.5.nc"

;----------------------------------------------------
procedure create_nc_file(x:float,TIME,LAT,LON,DIR:string,VAR:string)
local x,LON,LAT,TIME,DIR,VAR,fout
begin
    system("/bin/rm -f "+DIR)

  fout = addfile(DIR,"c")

  dimNames = (/"time", "lat", "lon"/)  
  dimSizes = (/ -1   ,  dimsizes(LAT),  dimsizes(LON)/) 
  dimUnlim = (/ True , False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout, "time" ,typeof(TIME),getvardims(TIME))                           
  filevardef(fout, "lat"  ,typeof(LAT),getvardims(LAT))                          
  filevardef(fout, "lon"  ,typeof(LON),getvardims(LON))                          
  filevardef(fout, VAR    ,typeof(x)  ,getvardims(x))    
                                                                                             
  filevarattdef(fout,"time" ,TIME)                    
  filevarattdef(fout,"lat"  ,LAT)                   
  filevarattdef(fout,"lon"  ,LON)                     
  filevarattdef(fout,VAR   ,x)   

  setfileoption(fout,"DefineMode",False)

  fout->time    = (/TIME/)     
  fout->lat     = (/LAT/)
  fout->lon     = (/LON/) 
  fout->$VAR$   = (/x/)
  
end

;----------------------------------------------------
; Main code.
;----------------------------------------------------
begin

  fu = addfile(udir,"r")
  fv = addfile(vdir,"r")

  u  = fu->u(:,0,:,:)
  v  = fv->v(:,0,:,:)

  print("calc divwind begin. "+systemfunc("date"))
  div = uv2dvF_Wrap(u, v)
  dvw = dv2uvF_Wrap(div)
  delete(div)
  print("calc divwind done. "+systemfunc("date"))

  print("store divwind begin. "+systemfunc("date"))
  dvu = dvw(0,:,:,:)
  dvv = dvw(1,:,:,:)
  delete(dvw)
  copy_VarCoords(u, dvu)
  copy_VarCoords(v, dvv)
  dvu@units = "m/s"
  dvu@long_name = "divergent zonal wind"
  dvv@units = "m/s"
  dvv@long_name = "divergent meridional wind"
  create_nc_file(dvu,u&time,u&lat,u&lon,dvudir,"u")
  create_nc_file(dvv,v&time,v&lat,v&lon,dvvdir,"v")
  delete([/dvu,dvv/])
  print("store divwind done. "+systemfunc("date"))

  print("calc vorticity begin. "+systemfunc("date"))
  vor = uv2vrF_Wrap(u, v)
  create_nc_file(vor,u&time,u&lat,u&lon,odir,"vor")
  print("store vorticity done. "+systemfunc("date"))

end