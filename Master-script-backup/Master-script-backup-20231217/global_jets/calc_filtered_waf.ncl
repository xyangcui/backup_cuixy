;-----------------------------------------------------------------
;; this script is just used to filter values by removing seasonal
;; cycle. Lanczos filter.
;; time: 2023-11-22
;-----------------------------------------------------------------
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_3D.ncl"


undef("create_nc_file")

procedure create_nc_file(x:float,y:float,TIME,LAT,LON,DIR:string,VAR1:string,VAR2:string)
local x,LON,LAT,TIME,DIR,VAR1,VAR2,fout,y
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
  filevardef(fout, VAR1    ,typeof(x)  ,getvardims(x))   
  filevardef(fout, VAR2    ,typeof(y)  ,getvardims(y))   
                                                                                             
  filevarattdef(fout,"time" ,TIME)                    
  filevarattdef(fout,"latitude"  ,LAT)                   
  filevarattdef(fout,"longitude"  ,LON)                     
  filevarattdef(fout,VAR1   ,x) 
  filevarattdef(fout,VAR2   ,y)   

  setfileoption(fout,"DefineMode",False)

  fout->time         = (/TIME/)     
  fout->latitude     = (/LAT/)
  fout->longitude    = (/LON/) 
  fout->$VAR1$       = (/x/)
  fout->$VAR2$       = (/y/)
  
end

;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
begin
  var  = "waf"
  nwa  = 100.
  nwb  = 20.

  ddir = "/home/sunming/data5/cuixy/global_jets/data"

  fin  = addfile(ddir+"/ERA5_daily_"+var+"_1979-2022.nc","r")

  outdir = ddir+"/"+"ERA5_anom_waf_1979-2022.nc"

  time = fin->time
  lat  = fin->latitude
  lon  = fin->longitude

  xu   = fin->u 
  u   = iso_filter(xu,201,nwa,nwb,0)

  xv   = fin->v
  v   = iso_filter(xv,201,nwa,nwb,0)

  print("calc done. "+ systemfunc("date"))
  copy_VarCoords(xv,v)
  copy_VarCoords(xu,u)

  create_nc_file(u,v,time,lat,lon,outdir,"u","v")

end

