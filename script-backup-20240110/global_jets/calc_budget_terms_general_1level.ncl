load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_3D.ncl"
;---------------------------------------------------
;; calc temperature budget.
;; temp, uwind, vwind, omega.
;; 5 terms: dTdt, udTdx, vdTdy, wdTdp, Q/cp.
;---------------------------------------------------
; constant.
pi  = atan(1.0)*4.
a   = 6371000.
td  = 24.*60.*60.
cp  = 1004.
latmin = 0
latmax = 60
lonmin = 0
lonmax = 200

; dir
ddir = "/home/sunming/data5/cuixy/global_jets/data"
tdir = ddir + "/ERA5_daily_mean_tlevel"
udir = ddir+"/"+"ERA5_daily_mean_ulevel_1979-2022.nc"
vdir = ddir+"/"+"ERA5_daily_mean_vlevel_1979-2022.nc"
wdir = ddir+"/"+"ERA5_daily_mean_wlevel_1979-2022.nc"
files  = systemfunc("ls " + tdir + "/day_t_*.nc")
;---------------------------------------------------
; function of vertical average.
;---------------------------------------------------
function vertavg(var)
local pres,opt,out
begin
  pres = var(:,0,:,:)
  pres = 1000.

  opt  = (/ 1, 250, 500 /)
  out  = wgt_vert_avg_beta(var&level, var, pres, 0, opt)
  copy_VarCoords(var(:,0,:,:), out)
  return out
end
;---------------------------------------------------

;---------------------------------------------------
; main code.
;---------------------------------------------------
  fu = addfile(udir,"r")
  fv = addfile(vdir,"r")
  fw = addfile(wdir,"r")
  ft = addfiles(files,"r")
  ListSetType (ft, "cat")

  t    = ft[:]->t(:,{250},{latmin:latmax},{lonmin:lonmax})
  u    = fu->u(:,{250},{latmin:latmax},{lonmin:lonmax})
  time = u&time
  lat  = u&latitude
  lon  = u&longitude

  t&time = time
print("sel T done. "+ systemfunc("date"))
; term 1.
  
  dTdt = center_finite_diff_n(t, 1., False, 0, 0)
  copy_VarCoords(t,dTdt)
  dTdt@unit = "K/day"
  copy_VarCoords(t,dTdt)

print("calc term1 done. "+ systemfunc("date"))  
; term 2 and term 3.
  v    = fv->v(:,{250},{latmin:latmax},{lonmin:lonmax})   

  ;gradLatLon  = grad_latlon_cfd (t, t&latitude, t&longitude, True, False)
  ;T_grad_lat  = gradLatLon[0]    
  ;T_grad_lon  = gradLatLon[1]
  ;delete(gradLatLon)
  dlon = (lon(2)-lon(1))*0.0174533 ; convert to radia

  dTdX = new ( dimsizes(t), typeof(t), t@_FillValue)
  dTdX = 0.
  do nl=0,dimsizes(lat)-1                      ; loop over each latitude
      dX = 6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
      dTdX(:,nl:nl,:) = center_finite_diff_n (t(:,nl:nl,:), dX , True,0,2)
  end do

  dTdY = new ( dimsizes(t), typeof(t), t@_FillValue)
  dTdY = 0.
  do nl=0,dimsizes(lon)-1                      ; loop over each latitude
      dY = 6378388.*dlon  ; constant at this latitude
      dTdY(:,:,nl:nl) = center_finite_diff_n (t(:,:,nl:nl), dY , False,0,1)
  end do

  udTdx = u*dTdX
  vdTdy = v*dTdY

  delete([/u,v,dTdX,dTdY/])

  copy_VarCoords(t,udTdx)
  copy_VarCoords(t,vdTdy)

  udTdx@long_name = "zonal advection"
  vdTdy@long_name = "meridinoal advection"
  udTdx@unit = "K/s"
  vdTdy@unit = "K/s"

print("calc term2 and term3 done. "+ systemfunc("date"))
; term 3.
  t1   = ft[:]->t(:,{200:300},{latmin:latmax},{lonmin:lonmax})
  t1&time = time
  w    = fw->w(:,{250},{latmin:latmax},{lonmin:lonmax})
;; dTdp.
  R  = 287.1
  cp = 1004.0
  dTdp  = center_finite_diff_n(t1,t1&level*100,False,0,1)
  copy_VarCoords(t1,dTdp)
  delete(t1)
  dTdp := dTdp(:,{250},:,:)

  ttmp = t*R/(cp*25000.)
  ss = dTdp-ttmp
  wdtdp = w*ss
  copy_VarCoords(t,wdtdp)
  delete([/w/])
  wdtdp@long_name = "vertical advection"
  wdtdp@unit = "K/s"
  wdtdp1 = tofloat(wdtdp)
  copy_VarCoords(t,wdtdp1)
print("calc term4 done. "+ systemfunc("date"))

; output.
  t1dir = ddir+"/temp_budget/"+"mean_dtdt.nc"
  create_nc_file(dTdt,time,lat,lon,t1dir,"dtdt")
print("store term1 done. "+ systemfunc("date"))

  t2dir = ddir+"/temp_budget/"+"mean_udtdx.nc"
  create_nc_file(udTdx,time,lat,lon,t2dir,"udtdx")
print("store term2 done. "+ systemfunc("date"))

  t3dir = ddir+"/temp_budget/"+"mean_vdtdy.nc"
  create_nc_file(vdTdy,time,lat,lon,t3dir,"vdtdy")
print("store term3 done. "+ systemfunc("date"))

  t4dir = ddir+"/temp_budget/"+"mean_wdtdp.nc"
  ;system("/bin/rm -f "+t4dir)

  ;fout = addfile(t4dir,"c")
  ;fout->wdtdp = wdtdp1
  create_nc_file(wdtdp1,time,lat,lon,t4dir,"wdtdp")
print("store term4 done. "+ systemfunc("date"))