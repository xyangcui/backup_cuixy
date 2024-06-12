load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_4D.ncl"
;-----------------------------------------------------------------
;; this script is just used to filter values by removing seasonal
;; cycle. Lanczos filter.
;; Logicical: read in single level values.
;; calculate filtered values.
;; selected one year value into a martirx.
;; calc all level and form the martrix and save as a file.
;; time: 2023-11-29
;-----------------------------------------------------------------

  ddir = "/home/sunming/data5/cuixy/global_jets/data"
  idir = ddir + "/ERA5_daily_mean_tlevel"
  odir = ddir + "/ERA5_daily_anom_tlevel"
  tdir = ddir + "/ERA5_daily_t250_1979-2022_all.nc"

  lats = 0
  latn = 60
  lonw = 0
  lone = 200
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------

  files = systemfunc("ls " + idir + "/day_t_flt_*.nc")
  fin   = addfiles(files,"r")
  ListSetType (fin, "cat")
  ft    = addfile(tdir,"r")
  time  = ft->time
  level = fin[0]->level;({200,250,300,350,400,450,500,600,700,850,925,1000})
  lat   = fin[0]->latitude({latmin:latmax})
  lon   = fin[0]->longitude

  TIME   = cd_calendar(time, 0)
  year   = toint(TIME(:,0))
  nt     = ind(year.eq.j)

  time_t = time(nt)
  xp = new((/dimsizes(time_t),dimsizes(level),\
              dimsizes(lat),dimsizes(lon)/),float)
  xp!0 = "time"
  xp!1 = "level"
  xp!2 = "lat"
  xp!3 = "lon"
  xp&time  = time_t
  xp&level = level
  xp&latitude  = lat
  xp&longitude = lon
  xp = 0.0
  do i = 0, dimsizes(level)-1
      print("level " + i + " begin. "+systemfunc("date"))
      x  = fin[:]->$var$(:,i,{latmin:latmax},:)
      x&time = time
      xa = iso_filter(x,201,nwa,nwb,0)

      xp(:,i,:,:) = xa(nt,:,:)
      delete([/xa,x/])
      print("level " + i + " done. "+systemfunc("date"))
      
  end do
    
  create_nc_file(xp,time_t,level,lat,lon,odir+\
                   "/day_t_"+j+".nc",var)

  delete([/files,fin,time,level,lat,lon,TIME,year,nt,xp,time_t/])

