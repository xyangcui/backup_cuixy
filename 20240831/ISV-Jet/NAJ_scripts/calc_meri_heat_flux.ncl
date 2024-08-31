;-----------------------------------------------------------------
;; this script is just used to filter values by removing seasonal
;; cycle. Lanczos filter.
;; time: 2023-12-18.
;-----------------------------------------------------------------
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_3D.ncl"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
  lats =  10
  latn =  80
  lonw = -120
  lone =  0

  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  fv  = addfile(ddir+"/ERA5_daily_mean_v850_1979-2022_r2.5.nc","r")
  ft  = addfile(ddir+"/ERA5_daily_mean_t850_1979-2022_r2.5.nc","r")

  v   = fv->v(:,{850},{lats:latn},{lonw:lone})
  t   = ft->t(:,{850},{lats:latn},{lonw:lone})

  time = v&time
  lat  = v&lat
  lon  = v&lon

  wgth = filwgts_lanczos(101, 1, 1./8., -999., 1.)
  vh   = wgt_runave_n_Wrap(v, wgth, 0, 0)
  copy_VarCoords(v, vh)
  th   = wgt_runave_n_Wrap(t, wgth, 0, 0)
  copy_VarCoords(t, th)

;; 
  x = vh*th
  copy_VarCoords(t, x)
  create_nc_file(x,time,lat,lon,\
                   ddir+"/ERA5_daily_mean_mhf850_1979-2022_r2.5.nc","mhf")

  xa   = iso_filter(x,201,30.,10.,0)
  copy_VarCoords(x, xa)

  create_nc_file(xa,time,lat,lon,\
                   ddir+"/ERA5_daily_anom_mhf850_1979-2022_r2.5.nc","mhf")

  print("calc high done. "+ systemfunc("date"))

;; calc mht index.
; time table.
  idir = ddir + "/EKE250_daily_anom_1979-2022_r2.5.nc"
  f = addfile(idir,"r")
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.4.and.day.eq.1)

  mon := mon(n1:n2)
  xa  := xa(n1:n2,:,:)
  nt   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  mht  = xa(nt,{35:50},{-70:-40})

  rad = 4.0*atan(1.0)/180.0
  clat = cos(mht&lat * rad)
  mind = wgt_areaave_Wrap(mht, clat, 1.0, 1)

  mind := dim_standardize_n_Wrap(mind, 0, 0)

  filename = ddir+"/meri_heatflux.csv"
  system("rm -rf " + filename)
  write_table(filename,"w",[/mind/],"%5.2f")
