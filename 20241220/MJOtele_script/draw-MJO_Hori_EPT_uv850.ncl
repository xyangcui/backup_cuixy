
;--------------------------------------------------------------------------
;function: 计算每种事件并返回矩阵的程序。
;--------------------------------------------------------------------------
function calc_pro(timt:integer,xdir,zt:float,pht:integer,ampt:float)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,igood,pt,at,ip,out,pct1
begin
  dims = dimsizes(zt)
  dim2 = dims(1)
  dim3 = dims(2)
  pct1 = new((/8,dim2,dim3/),float) 
  pct1 = 0.
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nf = asciiread(xdir,-1,"string")
  is = str_get_field(nf,1,delim)
  ie = str_get_field(nf,2,delim)
  is := toint(is(1::))
  ie := toint(ie(1::))
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ;event all days.
  nda := ndtooned(nda)
  igood = ind(.not.ismissing(nda))
  nda := nda(igood)  ;index for z.
  pt = pht(nda)
  at = ampt(nda)
  pt  := pt(ind(at.gt.1.))
  nda := nda(ind(at.gt.1.))
  delete(igood)
  ;select amp >1. and phase.
  do m = 0, 7
    np = m+1
    ip = ind(pt.eq.m+1)
    ip@_FillValue = -999
    if (any(ip.eq.ip@_FillValue).eq."Missing")
      continue
    else
      pct1(m,:,:) = dim_avg_n(zt(nda(ip),:,:),0)
      delete(ip)
    end if
  end do

  return pct1

end
;------------------------------------------------------------------

  lats = -20
  latn =  20
  lonw =  90
  lone =  220

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  odir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
  udir = ddir+"/vs/"+"ERA5_"+"u"+"level_daily_anom_1979-2022_r2.5.nc"
  vdir = ddir+"/vs/"+"ERA5_"+"v"+"level_daily_anom_1979-2022_r2.5.nc"
  edir = ddir+"/vs/"+"ERA5_"+"ept"+"level_daily_anom_1979-2022_r2.5.nc"
;;jet events
  ndir = ddir+"/"+"IOnP_date_omi.csv"
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"
  idir = ddir+"/"+"IOiP_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
;-----------------------------------------------------------------
; main code.
;-----------------------------------------------------------------
  f = addfile(odir,"r")
  fu = addfile(udir,"r")
  fv = addfile(vdir,"r")
  fe = addfile(edir,"r")
  olr = f->olr(:,{lats:latn},{lonw:lone})
  u = fu->u(:,{850},{lats:latn},{lonw:lone})
  v = fv->v(:,{850},{lats:latn},{lonw:lone})
  e = fe->ept(:,:,{lats:latn},{lonw:lone})
  cii = e(:,{850},:,:) - e(:,{400},:,:)
  copy_VarCoords(e(:,0,:,:),cii)

  time = olr&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  timi = cd_calendar(time, -2) 
;;OMI
  delim = ","
  rmmf  = asciiread(omidir,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
  phaset = str_get_field(rmmf,4,delim)
  ampt = str_get_field(rmmf,7,delim)

  iyear := toint(iyear(1::))
  imon  := toint(imon(1::))
  iday  := toint(iday(1::))
  phaset := tointeger(phaset(1::))
  ampt := tofloat(ampt(1::))

  phase = new((/dimsizes(time)/),integer)
  amp   = new((/dimsizes(time)/),float)
  phase = 0
  amp   = 0.

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i)\
             .and.iday.eq.day(i))
    phase(i) = phaset(nt)
    amp(i) = ampt(nt)
    delete(nt)
  end do
  delete([/phaset,ampt,iyear,imon,iday/])
  print("complete OMI selling "+ systemfunc("date"))

;;MJO events.
  IOnP = new((/4,8,dimsizes(u&lat),dimsizes(u&lon)/),float)
  IOnP = 0.
  IOfP = IOnP
  IOsP = IOfP
  IOiP = IOsP
;;jet events (year,mon,day) (events,days,lat,lon)
  print("jet events begin. "+ systemfunc("date"))
  IOnP(0,:,:,:) = calc_pro(timi,ndir,olr,phase,amp) ;IOnP
  IOiP(0,:,:,:) = calc_pro(timi,idir,olr,phase,amp) ;IOiP
  IOfP(0,:,:,:) = calc_pro(timi,fdir,olr,phase,amp) ;IOfP
  IOsP(0,:,:,:) = calc_pro(timi,sdir,olr,phase,amp) ;IOsP

  IOnP(1,:,:,:) = calc_pro(timi,ndir,u,phase,amp) ;IOnP
  IOiP(1,:,:,:) = calc_pro(timi,idir,u,phase,amp) ;IOiP
  IOfP(1,:,:,:) = calc_pro(timi,fdir,u,phase,amp) ;IOfP
  IOsP(1,:,:,:) = calc_pro(timi,sdir,u,phase,amp) ;IOsP

  IOnP(2,:,:,:) = calc_pro(timi,ndir,v,phase,amp) ;IOnP
  IOiP(2,:,:,:) = calc_pro(timi,idir,v,phase,amp) ;IOiP
  IOfP(2,:,:,:) = calc_pro(timi,fdir,v,phase,amp) ;IOfP
  IOsP(2,:,:,:) = calc_pro(timi,sdir,v,phase,amp) ;IOsP

  IOnP(3,:,:,:) = calc_pro(timi,ndir,cii,phase,amp) ;IOnP
  IOiP(3,:,:,:) = calc_pro(timi,idir,cii,phase,amp) ;IOiP
  IOfP(3,:,:,:) = calc_pro(timi,fdir,cii,phase,amp) ;IOfP
  IOsP(3,:,:,:) = calc_pro(timi,sdir,cii,phase,amp) ;IOsP
  print("jet events done. "+ systemfunc("date"))

  IOsP(3,:,:,:) = smth9(IOsP(3,:,:,:), 0.5, 0.25, False)
  IOfP(3,:,:,:) = smth9(IOfP(3,:,:,:), 0.5, 0.25, False)
  IOnP(3,:,:,:) = smth9(IOnP(3,:,:,:), 0.5, 0.25, False)
  IOiP(3,:,:,:) = smth9(IOiP(3,:,:,:), 0.5, 0.25, False)

  IOP = (IOnP+IOiP+IOfP+IOsP)/4.
  copy_VarCoords(olr(0,:,:),IOP(0,0,:,:))
  copy_VarCoords(olr(0,:,:),IOfP(0,0,:,:))
  copy_VarCoords(olr(0,:,:),IOsP(0,0,:,:))
  copy_VarCoords(olr(0,:,:),IOnP(0,0,:,:))
  copy_VarCoords(olr(0,:,:),IOiP(0,0,:,:))
;--------------------------------------------------------------------
  wks  = gsn_open_wks ("x11", "hov" )  

  cmap = read_colormap_file("BlueWhiteOrangeRed")         
;shading u plot
  res = True
  res@gsnFrame           = False
  res@gsnDraw            = False
  res@gsnMaximize        = False
  res@gsnAddCyclic       = False
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30
  ;res@gsnRightStringFontHeightF = 0.015 
  res@gsnLeftStringFontHeightF  = 0.015

  ;res@gsnRightStringOrthogonalPosF = 0.01
  res@gsnLeftStringOrthogonalPosF  = 0.01

  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBLabelFont = 21
  res@tmYLLabelFont = 21
  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01

  res@mpMinLatF = lats
  res@mpMaxLatF = latn
  res@mpMinLonF = lonw
  res@mpMaxLonF = lone
  res@mpCenterLonF = 150
  res@mpFillOn        = False
  res@mpLandFillColor = "white"
  res@mpGridAndLimbDrawOrder = "PostDraw"
  res@mpGridLineColor        = "Gray60"

  res@cnFillOn               = True
  res@cnFillPalette          = cmap
  res@cnLinesOn              = False
  res@cnLineLabelsOn         = False
  res@cnLevelSelectionMode   = "Manuallevels"
  res@cnMinLevelValF         = -1.2
  res@cnMaxLevelValF         =  1.2
  res@cnLevelSpacingF        =  0.2
  res@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  res@lbLabelBarOn           = True

;; vector res.
  resw = True
  resw@gsnFrame     = False
  resw@gsnDraw      = False
  resw@gsnMaximize  = False
  resw@gsnAddCyclic = False
  resw@vcGlyphStyle            = "LineArrow"
  resw@vcLineArrowThicknessF   = 2.5
  resw@vcRefMagnitudeF         = 5.
  resw@vcRefLengthF            = 0.08
  resw@vcLineArrowColor        = "mediumblue"
  resw@vcRefAnnoString2On      = False
  resw@vcRefAnnoString1        = "$VMG$" +" m s~S~-1~N~"
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoFontHeightF    = 0.02
  resw@vcRefAnnoBackgroundColor = -1


  label = (/"p1","p2","p3","p4","p5","p6","p7","p8"/)

  pv = IOsP
  copy_VarCoords(IOfP,pv)
  res@gsnLeftString = "~F21~"+"IOsP"
  plot  = gsn_csm_contour_map(wks, pv(3,5,:,:), res) 
  plotv = gsn_csm_vector(wks, pv(1,5,::2,::2), pv(2,5,::2,::2), resw)

  overlay(plot, plotv)
  draw(plot)

  i = 0
  do while(i.eq.0)
    i = 0
  end do