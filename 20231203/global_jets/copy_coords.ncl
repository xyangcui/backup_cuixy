;--------------------------------------------------------
; calculate streamfunction by sperical function.
; Note that u and v must global data.
; author: cui xiangyang time: 2023-11-21
;--------------------------------------------------------
begin
  var  = "phi"
  ddir = "/home/sunming/data5/cuixy/global_jets/data"

  fin1  = addfile(ddir+"/"+"ERA5_daily_"+var+"250_1979-2022_all.nc","r")
  fin   = addfile(ddir+"/"+"ERA5_daily_anom_"+var+"250_1979-2022.nc","r")

  x    = fin->$var$
  time = fin1->time
  lat  = fin1->latitude
  lon  = fin1->longitude

  x!1 = "latitude"
  x!2 = "longitude"

  system("/bin/rm -f "+ddir+"/"+"ERA5_daily_anom_"+var+"250_1979-2022.nc")

  fout = addfile(ddir+"/"+"ERA5_daily_anom_"+var+"250_1979-2022.nc","c")

  dimNames = (/"time", "latitude", "longitude"/)  
  dimSizes = (/ -1   ,  dimsizes(lat),  dimsizes(lon)/) 
  dimUnlim = (/ True , False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout, "time" ,typeof(time),getvardims(time))                           
  filevardef(fout, "latitude"  ,typeof(lat),getvardims(lat))                          
  filevardef(fout, "longitude"  ,typeof(lon),getvardims(lon))                          
  filevardef(fout, var    ,typeof(x)  ,getvardims(x))    
                                                                                             
  filevarattdef(fout,"time" ,time)                    
  filevarattdef(fout,"latitude"   ,lat)                   
  filevarattdef(fout,"longitude"  ,lon)                     
  filevarattdef(fout,var   ,x)   

  setfileoption(fout,"DefineMode",False)

  fout->time         = (/time/)     
  fout->latitude     = (/lat/)
  fout->longitude    = (/lon/) 
  fout->$var$          = (/x/)
end