;--------------------------------------------
; ncep: 1999-2010 years
; grid: (10-80N 140-270E)
; previous 120 days

lats = 10
latn = 80
lonw = 140
lone = 270
ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/"
outdir = ddir+"/hindcast/ecmwf/2016/data/ecmwf_z500_daily_anom_1979-2022_r2.5.nc"

  f = addfile(ddir+"ERA5_z500_daily_mean_1979-2022_r2.5_with0229.nc","r")
  time = f->time
  TIME = cd_calendar(time,0)
  year = toint(TIME(:,0))
  nt = ind(year.ge.1996.and.year.le.2015)

  zc = f->z(nt,{500},{lats:latn},{lonw:lone})
  delete([/time,TIME,year,nt/])
  lat  = zc&lat
  lon  = zc&lon
  time = zc&time
  TIME = cd_calendar(time,0)
  year = toint(TIME(:,0))
  mont = toint(TIME(:,1))
  dayt = toint(TIME(:,2))
  ddd  = day_of_year(year, mont, dayt)
  yyyyddd = year*1000 + ddd
  hClmDay = clmDayTLL(zc, yyyyddd) 
  delete([/zc,time,year,TIME,mont,dayt,ddd,yyyyddd/])

  z = f->z(:,{500},{lats:latn},{lonw:lone}) 
  time = z&time
  TIME = cd_calendar(time,0)
  year = toint(TIME(:,0))
  mont = toint(TIME(:,1))
  dayt = toint(TIME(:,2))
  ddd  = day_of_year(year, mont, dayt)
  yyyyddd = year*1000 + ddd                                                                                                                               
  hAnom = calcDayAnomTLL(z, yyyyddd, hClmDay)     

  system("/bin/rm -f "+outdir)

  setfileoption("nc", "Format", "NetCDF4Classic")
  setfileoption("nc","CompressionLevel",5)

  fout = addfile(outdir,"c")
  setfileoption(fout,"DefineMode",True)
  
  fAtt               = True           
  fAtt@creation_date = systemfunc ("date")
  fileattdef( fout, fAtt )

  dimNames = (/"time", "lat", "lon"/)  
  dimSizes = (/ -1   ,  dimsizes(lat),  dimsizes(lon)/) 
  dimUnlim = (/ True , False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout, "time" ,typeof(time),getvardims(time))                           
  filevardef(fout, "lat"  ,typeof(lat),getvardims(lat))                          
  filevardef(fout, "lon"  ,typeof(lon),getvardims(lon))                          
  filevardef(fout, "z"    ,typeof(hAnom)  ,getvardims(hAnom))    
                                                                                             
  filevarattdef(fout,"time" ,time)                    
  filevarattdef(fout,"lat"  ,lat)                   
  filevarattdef(fout,"lon"  ,lon)                     
  filevarattdef(fout,"z" ,z)   

  setfileoption(fout,"DefineMode",False)

  fout->time  = (/time/)     
  fout->lat   = (/lat/)
  fout->lon   = (/lon/) 
  fout->z     = (/hAnom/)