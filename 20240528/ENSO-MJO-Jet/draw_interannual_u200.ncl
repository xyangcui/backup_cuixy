; ==============================================================

begin
  
  yrStrt = 19790101
  yrLast = 20221201

  lats = 0
  latn = 90

  season = "DJF"    

  ddir = "/home/sunming/data5/cuixy/ENSO-MJO-Jet/data"

  f  = addfile (ddir+"/ERA5_u200_mon_mean_1979-2022_r1.nc", "r")
  ;f  = addfile ("/home/sunming/data5/cuixy/global_jets/data/uwnd.mon.mean.nc", "r")

  YYYY   = cd_calendar(f->time,-2) ;cd_calendar(f->time,-2)   
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  print(iYYYY)  
  u  = f->uwnd(iYYYY,{200},{lats:latn},:)
        
  u := month_to_season (u, season)
  printVarSummary(u)
  um = dim_avg_n_Wrap(u, 0)

;; FFT for remove the first three harmonics.
  cf = ezfftf_n(u, 0)
  cf(:,0:2,:,:) = 0.
  printVarSummary(cf)
  ui = ezfftb_n(cf, cf@xbar, 1)

  uv = dim_variance_n(ui, 0)
  copy_VarCoords(u(0,:,:), uv)
  printVarSummary(uv)

;===========================================================
; PLOTS
;============================================================
  wks = gsn_open_wks("x11","eof")         
  plot  = new(1,graphic)                
  plotu = plot                                       

;--plot
  res = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False
  res@gsnAddCyclic = True
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringFontHeightF = 0.015
  res@gsnLeftStringOrthogonalPosF  = 0.01
  res@gsnRightStringOrthogonalPosF = 0.01

  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01

  res@mpFillOn = True
  res@mpLandFillColor = "white"
  res@mpCenterLonF      = 180.
  res@mpMinLatF         = 10
  res@mpMaxLatF         = 70.

  res@cnFillOn = True
  ;res@cnLabelBarEndStyle   = "IncludeMinMaxLabels"
  res@cnFillPalette        = "WhiteBlueGreenYellowRed"
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "Manuallevels"
  res@cnMinLevelValF       =  10
  res@cnMaxLevelValF       =  100
  res@cnLevelSpacingF      =  10

  res@lbBoxLinesOn          = True
  res@lbBoxSeparatorLinesOn = False
  ;res@lbTitleString         = "~F21~(unit: m~/s)"
  ;res@lbTitlePosition       = "Bottom"
  ;res@lbTitleFontHeightF    = 0.01
  ;res@lbOrientation         = "Vertical"
  ;res@pmLabelBarHeightF     = 0.2
  ;res@pmLabelBarWidthF      = 0.1
  res@pmLabelBarOrthogonalPosF = 0.2

  res@tmXBLabelFont = 21
  res@tmYRLabelFont = 21

  res@gsnLeftString    = "~F21~(a) ~F21~200hPa U"
  res@gsnRightString   = "~F21~DJF"

  resc = True
  resc@gsnFrame     = False
  resc@gsnDraw      = False
  resc@gsnMaximize  = False
  resc@gsnAddCyclic = True
  resc@gsnContourNegLineDashPattern = 1
  resc@gsnLeftString      = ""
  resc@gsnRightString     = ""
  resc@cnMonoLineColor    = True
  resc@cnLinesOn          = True
  resc@cnLineLabelsOn     = False
  resc@cnLineColor        = "red"
  resc@cnLineThicknessF   = 2.0
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels      = (/30.,40.,50.,60.,70./)


  plot = gsn_csm_contour_map(wks, uv, res)

  plotu = gsn_csm_contour(wks, um, resc)

  overlay(plot, plotu)

  resp = True
  resp@gsLineColor = "black"
  resp@gsLineThicknessF = 3.0

  resp@gsLineDashPattern = 0
  dum1 = gsn_add_polyline(wks, plot, (/45,45,75,75,45/), (/25,35,35,25,25/), resp)

  resp@gsLineDashPattern = 1
  dum2 = gsn_add_polyline(wks, plot, (/170,170,210,210,170/), (/25,40,40,25,25/), resp)

  resp@gsLineDashPattern = 2
  dum3 = gsn_add_polyline(wks, plot, (/-60,-60,-25,-25,-60/), (/30,55,55,30,30/), resp)

  draw(plot)

  i = 0
  do while (i.eq.0)
    i = 0
  end do

end