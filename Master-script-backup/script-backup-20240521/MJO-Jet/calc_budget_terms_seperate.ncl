load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_3D.ncl"
;load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_4D.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
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

  lats = 10
  latn = 50
  lonw = -160
  lone =  20

; dir
  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJS_temp_budget"
  tdir = ddir+"/"+"ERA5_daily_mean_tlevel_1979-2022_r2.5.nc"
  udir = ddir+"/"+"ERA5_daily_mean_ulevel_1979-2022_r2.5.nc"
  ;vdir = ddir+"/"+"ERA5_daily_mean_vlevel_1979-2022_r2.5.nc"
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
  print("prepare data"+" begin. "+ systemfunc("date"))
  fu = addfile(udir,"r")
  ft = addfile(tdir,"r")
  ;fv = addfile(vdir,"r")

  t1    = ft->t(:,{250:500},{lats:latn},{lonw:lone})
  u1    = fu->u(:,{250:500},{lats:latn},{lonw:lone})
  ;u1    = fv->v(:,{250:500},{lats:latn},{lonw:lone})

  printVarSummary(u1)
  time = u1&time
  lat  = u1&lat
  lon  = u1&lon
  level = u1&level

  wgtl = filwgts_lancos(201, 0, 1/100., -999., 1.)
  wgti = filwgts_lancos(201, 2, 1/100., 1/20., 1.)
  wgth = filwgts_lancos(201, 1, 1/20., -999., 1.)
  print("prepare data"+" done. "+ systemfunc("date"))

  ;u2 = wgt_runave_n_Wrap(u1, wgtl, 0, 0)
  ;t2 = wgt_runave_n_Wrap(t1, wgti, 0, 0)

  ;u  = vertavg(u2)
  ;copy_VarCoords(u1(:,0,:,:),u)

  ;t  = vertavg(t2)
  ;copy_VarCoords(t1(:,0,:,:),t)

  ;uoutdir = ddir+"/"+"ERA5_daily_low_uvert_1979-2022_r2.5.nc" 
  ;toutdir = ddir+"/"+"ERA5_daily_iso_tvert_1979-2022_r2.5.nc" 

  ;create_nc_file(u,time,lat,lon,uoutdir,"u")

  ;create_nc_file(t,time,lat,lon,toutdir,"t")

  print("prepare data"+" done. "+ systemfunc("date"))

  if (varname.eq."ultl") then
    u = wgt_runave_n_Wrap(u1, wgtl, 0, 0)
    t = wgt_runave_n_Wrap(t1, wgtl, 0, 0)
    else if (varname.eq."ulti") then
      u = wgt_runave_n_Wrap(u1, wgtl, 0, 0)
      t = wgt_runave_n_Wrap(t1, wgti, 0, 0)
    else if (varname.eq."ulth") then
      u = wgt_runave_n_Wrap(u1, wgtl, 0, 0)
      t = wgt_runave_n_Wrap(t1, wgth, 0, 0)
    else if (varname.eq."uitl") then
      u = wgt_runave_n_Wrap(u1, wgti, 0, 0)
      t = wgt_runave_n_Wrap(t1, wgtl, 0, 0)
    else if (varname.eq."uiti") then
      u = wgt_runave_n_Wrap(u1, wgti, 0, 0)
      t = wgt_runave_n_Wrap(t1, wgti, 0, 0)
    else if (varname.eq."uith") then
      u = wgt_runave_n_Wrap(u1, wgti, 0, 0)
      t = wgt_runave_n_Wrap(t1, wgth, 0, 0)
    else if (varname.eq."uhtl") then
      u = wgt_runave_n_Wrap(u1, wgth, 0, 0)
      t = wgt_runave_n_Wrap(t1, wgtl, 0, 0)
    else if (varname.eq."uhti") then
      u = wgt_runave_n_Wrap(u1, wgth, 0, 0)
      t = wgt_runave_n_Wrap(t1, wgti, 0, 0)
  else
    u = wgt_runave_n_Wrap(u1, wgth, 0, 0)
    t = wgt_runave_n_Wrap(t1, wgth, 0, 0)
    end if
    end if
    end if
    end if
    end if
    end if
    end if
  end if

  print("prepare "+varname+" done. "+ systemfunc("date"))

  dlon = (lon(2)-lon(1))*0.0174533 ; convert to radia

  dTdX = new ( dimsizes(t1), typeof(t1), t@_FillValue)
  dTdX = 0.
  do nl=0,dimsizes(lat)-1                      ; loop over each latitude
    dX = tofloat(6378388.*cos(0.0174533*lat(nl))*dlon)  ; constant at this latitude
    dTdX(:,:,nl:nl,:) = center_finite_diff_n (t(:,:,nl:nl,:), dX , True,0,3)
  end do

  ;do nl=0,dimsizes(lon)-1                      ; loop over each latitude
  ;  dY = tofloat(6378388.*dlon)  ; constant at this latitude
  ;  dTdX(:,:,:,nl:nl) = center_finite_diff_n (t(:,:,:,nl:nl), dY , True,0,2)
  ;end do

  udTdxp = u*dTdX

  delete([/u,dTdX/])

  copy_VarCoords(t,udTdxp)

  udTdx  = vertavg(udTdxp)
  copy_VarCoords(t(:,0,:,:),udTdx)
  udTdx@long_name = "zonal advection"
  udTdx@unit = "K/s*Pa"
  delete(udTdxp)

  x = iso_filter(udTdx,201,100.,20.,0)

  print("calc "+varname+" done. "+ systemfunc("date"))

; output.
  t2dir = ddir+"/anom_"+varname+"_r2.5.nc"
  create_nc_file(x,time,lat,lon,t2dir,varname)

  print("store "+varname+" done. "+ systemfunc("date"))