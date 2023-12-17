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

  t    = ft[:]->t(:,{250:500},{latmin:latmax},{lonmin:lonmax})
  u    = fu->u(:,{250:500},{latmin:latmax},{lonmin:lonmax})
  time = u&time
  lat  = u&latitude
  lon  = u&longitude

  t&time = time
print("sel T done. "+ systemfunc("date"))
; term 1.
  
  dTdtp = center_finite_diff_n(t, 1., False, 0, 0)
  copy_VarCoords(t,dTdtp)
  dTdt = vertavg(dTdtp)
  dTdt@unit = "K/day"
  copy_VarCoords(t(:,0,:,:),dTdt)
  delete([/dTdtp/])

;print("calc term1 done. "+ systemfunc("date"))  
; term 2 and term 3.
  v    = fv->v(:,{250:500},{latmin:latmax},{lonmin:lonmax})   

  ;gradLatLon  = grad_latlon_cfd (t, t&latitude, t&longitude, True, False)
  ;dTdX  = gradLatLon[0]    
  ;dTdY  = gradLatLon[1]
  ;delete(gradLatLon)
  dlon = (lon(2)-lon(1))*0.0174533 ; convert to radia
  dlat = dlon

  dTdX = new ( dimsizes(t), typeof(t), t@_FillValue)
  dTdX = 0.
  do nl=0,dimsizes(lat)-1                      ; loop over each latitude
      dX = 6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
      dTdX(:,:,nl:nl,:) = center_finite_diff_n (t(:,:,nl:nl,:), dX , True,0,3)
  end do

  dY = 6378388.*dlat 
  dTdY = center_finite_diff_n (t, dY , False,0,2)

  udTdxp = u*dTdX
  vdTdyp = v*dTdY

  delete([/u,v,dTdX,dTdY/])

  copy_VarCoords(t,udTdxp)
  copy_VarCoords(t,vdTdyp)

  udTdx  = vertavg(udTdxp)
  vdTdy  = vertavg(vdTdyp)
  copy_VarCoords(t(:,0,:,:),udTdx)
  copy_VarCoords(t(:,0,:,:),vdTdy)
  udTdx@long_name = "zonal advection"
  vdTdy@long_name = "meridinoal advection"
  udTdx@unit = "K/s*Pa"
  vdTdy@unit = "K/s*Pa"
  delete([/udTdxp,vdTdyp/])

  print("calc term2 and term3 done. "+ systemfunc("date"))
; term 4.
  t1   = ft[:]->t(:,{200:600},{latmin:latmax},{lonmin:lonmax})
  t1&time = time
  w    = fw->w(:,{250:500},{latmin:latmax},{lonmin:lonmax})
;; dTdp.
  ;sstmp = static_stability(t1&level*100., t1, 1, 0)
  ;copy_VarCoords(t1,sstmp)
  ;ss = sstmp(:,{250:500},:,:)
  R  = 287.0
  cp = 1005.7
  dTdp  = center_finite_diff_n(t1,t1&level*100.,False,0,1)
  copy_VarCoords(t1,dTdp)
  dTdp := dTdp(:,{250:500},:,:)
  ss2 = dTdp
  ss2 = 0.
  ss2(:,{250},:,:) = t(:,{250},:,:)*R/(cp*25000.)
  ss2(:,{300},:,:) = t(:,{300},:,:)*R/(cp*30000.)
  ss2(:,{350},:,:) = t(:,{350},:,:)*R/(cp*35000.)
  ss2(:,{400},:,:) = t(:,{400},:,:)*R/(cp*40000.)
  ss2(:,{450},:,:) = t(:,{450},:,:)*R/(cp*45000.)
  ss2(:,{500},:,:) = t(:,{500},:,:)*R/(cp*50000.)

  print("calc ss done. "+ systemfunc("date"))
  wdtdpp1 = w*dTdp
  wdtdpp2 = w*ss2
  copy_VarCoords(t,wdtdpp1)
  copy_VarCoords(t,wdtdpp2)
  delete([/w/])

  wdtdp1 = vertavg(wdtdpp1)
  copy_VarCoords(t(:,0,:,:),wdtdp1)
  wdtdp1@long_name = "vertical advection"
  wdtdp1@unit = "K/s*Pa"
  wdtdp11 = tofloat(wdtdp1)
  copy_VarCoords(t(:,0,:,:),wdtdp11)

  wdtdp2 = vertavg(wdtdpp2)
  copy_VarCoords(t(:,0,:,:),wdtdp2)
  wdtdp2@long_name = "adiabatic heating"
  wdtdp2@unit = "K/s*Pa"
  wdtdp21 = tofloat(wdtdp2)
  copy_VarCoords(t(:,0,:,:),wdtdp21)

  delete([/wdtdp1,wdtdp2,wdtdpp1,wdtdpp2/])
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
  create_nc_file(wdtdp11,time,lat,lon,t4dir,"wdtdp")
  print("store term4 done. "+ systemfunc("date"))

  t5dir = ddir+"/temp_budget/"+"mean_wrtcp.nc"
  create_nc_file(wdtdp21,time,lat,lon,t5dir,"wrtcp")
  print("store term5 done. "+ systemfunc("date"))