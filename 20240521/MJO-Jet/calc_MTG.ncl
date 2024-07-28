load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_3D.ncl"
;----------------------------------------------
;; calculate meridinoal temperature gradient.
; author: cui xiangyang time: 2023-12-04
;----------------------------------------------

ddir   = "/home/sunming/data5/cuixy/global_jets/data"
pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

phidir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022.nc"
idir   = ddir + "/ERA5_daily_anom_tlevel"
files  = systemfunc("ls " + idir + "/day_t_*.nc")

var = "mtg"

;---------------------------------------------------
; Main code.
;---------------------------------------------------
begin

  fphi = addfile(phidir,"r")
  ft   = addfiles(files,"r")
  ListSetType (ft, "cat")

  t    = ft[:]->t(:,{250:500},:,:)

  print("complete data selling "+ systemfunc("date"))

  time = fphi->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  t&time = time
  lat    = t&latitude
  lon    = t&longitude
  level  = t&level

  ;calc MTG and vertical average.
                                
  gradLatLon  = grad_latlon_cfd (t, lat, lon, True, False)                               
  ori_MTG     = gradLatLon[0]    
  delete(gradLatLon)

  print("all layers complete. "+ systemfunc("date"))

  ori_MTG@long_name = "latitudinal gradient (derivative)"
  ori_MTG@units     = "K/m"

  copy_VarCoords(t,ori_MTG)

  pres = t(:,0,:,:)
  pres = 1000.

  opt  = (/ 1, 250, 500 /)
  MTG  = wgt_vert_avg_beta(level, ori_MTG, pres, 0, opt)

  print("average complete. "+ systemfunc("date"))

  copy_VarCoords(t(:,0,:,:), MTG)

  create_nc_file(MTG,time,lat,lon,ddir+"/ERA5_daily_anom_"\
                 +var+"_1979-2022.nc",var)

  print("sort complete. "+ systemfunc("date"))

end




