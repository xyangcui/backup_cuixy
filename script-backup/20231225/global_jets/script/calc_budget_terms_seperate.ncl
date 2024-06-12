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
tdir = ddir+"/temp_budget/"+"ERA5_daily_"+tin+"_t250_1979-2022.nc"
udir = ddir+"/temp_budget/"+"ERA5_daily_"+uin+"_u250_1979-2022.nc"
;vdir = ddir+"/"+"ERA5_daily_mean_vlevel_1979-2022.nc"
;wdir = ddir+"/"+"ERA5_daily_mean_wlevel_1979-2022.nc"
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
  ft = addfile(tdir,"r")

  t    = ft->t(:,{250:500},:,:)
  u    = fu->u(:,{250:500},:,:)
  time = u&time
  lat  = u&latitude
  lon  = u&longitude

print("sel T done. "+ systemfunc("date"))

;print("calc term1 done. "+ systemfunc("date"))  
; term 2 and term 3.  
  ;gradLatLon  = grad_latlon_cfd (t, t&latitude, t&longitude, True, False)
  ;dTdX  = gradLatLon[0]    
  ;delete(gradLatLon)
  dlon = (lon(2)-lon(1))*0.0174533 ; convert to radia

  dTdX = new ( dimsizes(t), typeof(t), t@_FillValue)
  dTdX = 0.
  do nl=0,dimsizes(lat)-1                      ; loop over each latitude
      dX = 6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
      dTdX(:,:,nl:nl,:) = center_finite_diff_n (t(:,:,nl:nl,:), dX , True,0,3)
  end do

  udTdxp = u*dTdX

  delete([/u,dTdX/])

  copy_VarCoords(t,udTdxp)

  udTdx  = vertavg(udTdxp)
  copy_VarCoords(t(:,0,:,:),udTdx)
  udTdx@long_name = "zonal advection"
  udTdx@unit = "K/s*Pa"
  delete(udTdxp)

  print("calc term2 and term3 done. "+ systemfunc("date"))

; output.
  t2dir = ddir+"/temp_budget/"+"mean_"+varname+".nc"
  create_nc_file(udTdx,time,lat,lon,t2dir,"udtdx")
  print("store term2 done. "+ systemfunc("date"))