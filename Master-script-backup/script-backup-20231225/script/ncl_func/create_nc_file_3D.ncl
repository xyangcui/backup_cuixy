;---------------------------------------------------------
;; use effective method to create netcdf file.
; variables info:
; x: required data (time,lat,lon)
; LON: longitude
; LAT: latitude
; TIME: time.
; DIR: outdir.
; VAR: variables name.
; default: x(time|:, lat|:, lon|:)
;---------------------------------------------------------

undef("create_nc_file")

procedure create_nc_file(x:float,TIME,LAT,LON,DIR:string,VAR:string)
local x,LON,LAT,TIME,DIR,VAR,fout
begin
    system("/bin/rm -f "+DIR)

  fout = addfile(DIR,"c")

  dimNames = (/"time", "latitude", "longitude"/)  
  dimSizes = (/ -1   ,  dimsizes(LAT),  dimsizes(LON)/) 
  dimUnlim = (/ True , False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout, "time" ,typeof(TIME),getvardims(TIME))                           
  filevardef(fout, "latitude"  ,typeof(LAT),getvardims(LAT))                          
  filevardef(fout, "longitude"  ,typeof(LON),getvardims(LON))                          
  filevardef(fout, VAR    ,typeof(x)  ,getvardims(x))    
                                                                                             
  filevarattdef(fout,"time" ,TIME)                    
  filevarattdef(fout,"latitude"  ,LAT)                   
  filevarattdef(fout,"longitude"  ,LON)                     
  filevarattdef(fout,VAR   ,x)   

  setfileoption(fout,"DefineMode",False)

  fout->time         = (/TIME/)     
  fout->latitude     = (/LAT/)
  fout->longitude    = (/LON/) 
  fout->$VAR$        = (/x/)
  
end
