;--------------------------------------------------------
; calculate streamfunction by sperical function.
; Note that u and v must global data.
; author: cui xiangyang time: 2023-11-21
;--------------------------------------------------------
begin
  ddir = "/home/sunming/data5/cuixy/global_jets/data"

  fu = addfile(ddir+"/"+"ERA5_daily_u250_1979-2022_all.nc","r")
  fv = addfile(ddir+"/"+"ERA5_daily_sf250_1979-2022_all1.nc","r")

  u  = fu->u(:,0,::-1,:)
  phi  = fv->phi
  time = u&time
  lat  = u&latitude
  lon  = u&longitude
  delete(u)
  phi!1 = "latitude"
  phi!2 = "longitude"

  system("/bin/rm -f "+ddir+"/"+"ERA5_daily_sf250_1979-2022_all.nc")

  fout = addfile(ddir+"/"+"ERA5_daily_sf250_1979-2022_all.nc","c")
 ;===================================================================
    ; predefine the coordinate variables and their dimensionality
    ; Note: to get an UNLIMITED record dimension, we set the dimensionality
    ; to -1 (or the actual size) and set the dimension name to True.
    ;===================================================================
  dimNames = (/"time", "latitude", "longitude"/)  
  dimSizes = (/ -1   ,  dimsizes(lat),  dimsizes(lon)/) 
  dimUnlim = (/ True , False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

    ;===================================================================
    ; predefine the the dimensionality of the variables to be written out
    ;===================================================================
  filevardef(fout, "time" ,typeof(time),getvardims(time))                           
  filevardef(fout, "latitude"  ,typeof(lat),getvardims(lat))                          
  filevardef(fout, "longitude"  ,typeof(lon),getvardims(lon))                          
  filevardef(fout, "phi"    ,typeof(phi)  ,getvardims(phi))    
                                                                                             
  filevarattdef(fout,"time" ,time)                    ; copy time attributes
  filevarattdef(fout,"latitude"   ,lat)                     ; copy lat attributes
  filevarattdef(fout,"longitude"  ,lon)                     ; copy lon attributes
  filevarattdef(fout,"phi"   ,phi)                      ; copy PS attributes

  setfileoption(fout,"DefineMode",False)

  fout->time         = (/time/)     
  fout->latitude     = (/lat/)
  fout->longitude    = (/lon/) 
  fout->phi          = (/phi/)
end