;-----------------------------------------------------------------
;; zonal average.
;-----------------------------------------------------------------

;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/vs/"

  fw = addfile(ddir+"ERA5_"+"w"+"level_daily_anom_1979-2022_r2.5.nc","r")
  x  = fw->w(:,{1000:100},{-10:5},:)
  w = dim_avg_n_Wrap(x, 2)
  delete(x)

  fu = addfile(ddir+"ERA5_"+"u"+"level_daily_anom_1979-2022_r2.5.nc","r")
  x  = fu->u(:,{1000:100},{-10:5},:)
  u = dim_avg_n_Wrap(x, 2)
  delete(x)

  fe = addfile(ddir+"ERA5_"+"ept"+"level_daily_anom_1979-2022_r2.5.nc","r")
  x  = fe->ept(:,{1000:100},{-10:5},:)
  e = dim_avg_n_Wrap(x, 2)
  delete(x)
  e@long_name = "EPT"
  e@units = "K"

  fh = addfile(ddir+"ERA5_"+"h"+"level_daily_anom_1979-2022_r2.5.nc","r")
  x  = fh->h(:,{1000:100},{-10:5},:)
  h = dim_avg_n_Wrap(x, 2)
  delete(x)
  h@long_name = "precipitation heating"
  h@units = "J/(kg s)"

  time  = u&time
  level = u&level
  lon   = u&lon

  system("/bin/rm -f "+ddir+"ERA5_uwlevel_zonalavg_r2.5.nc")
  fout = addfile(ddir+"ERA5_uwlevel_zonalavg_r2.5.nc","c")

  fout->w = w 
  fout->u = u 
  fout->e = e
  fout->h = h
  fout->time  = (/time/)
  fout->level = (/level/)
  fout->lon   = (/lon/)