load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_3D.ncl"

var ="u"
ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"

  ;f = addfile(ddir+"/NCEP_daily_mean_"+"z"+lev+"_1979-2022_r2.5.nc","r")
  ;f = addfile(ddir+"/NCEP_daily_mean_mslp_1979-2022_r2.5.nc","r")
  ;f = addfile(ddir+"/NOAA_daily_mean_olr_1979-2022_r1.nc","r")
  f = addfile(ddir+"/ERA5_daily_mean_u200_1979-2024.nc","r")
  ;f = addfile(ddir+"/OISST_daily_mean_sst_1982-2022.nc","r")
  TIME = cd_calendar(f->time,-2)
  nt = ind(TIME.ge.19821201.and.TIME.le.20220331)
  x = f->u(nt,:,:)
  delete([/TIME,nt/])   
  time    = x&time                       
  TIME    = cd_calendar(time, 0)            
  year    = toint( TIME(:,0) )               
  month   = toint( TIME(:,1) )
  day     = toint( TIME(:,2) ) 
  ddd     = day_of_year(year, month, day) 
  yyyyddd = year*1000 + ddd  

  hClmDay = clmDayTLL(x, yyyyddd) 
  hClmsmth= smthClmDayTLL(hClmDay, 3)

  x := f->u
  time    := x&time                       
  TIME    := cd_calendar(time, 0)            
  year    := toint( TIME(:,0) )               
  month   := toint( TIME(:,1) )
  day     := toint( TIME(:,2) ) 
  ddd     := day_of_year(year, month, day) 
  yyyyddd := year*1000 + ddd 

  xa = calcDayAnomTLL(x, yyyyddd, hClmDay)
  xa:= wgt_runave_n_Wrap(xa, 7., 0, 0)
  create_nc_file(xa,xa&time,xa&lat,xa&lon,ddir+\
                "/ERA5_daily_anom_"+"u"+"200"+"_1979-2022.nc",var)
  ;create_nc_file(xa,xa&time,xa&lat,xa&lon,ddir+\
  ;              "/NOAA_daily_anom_"+"olr"+"_1982-2022.nc",var)
  ;create_nc_file(xa,xa&time,xa&lat,xa&lon,ddir+\
  ;              "/NCEP_daily_anom_"+"mslp"+"_1979-2022.nc",var)
  ;create_nc_file(xa,xa&time,xa&lat,xa&lon,ddir+\
  ;              "/OISST_daily_anom_"+"sst"+"_1982-2022.nc",var)