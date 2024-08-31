load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
;-----------------------------------------------------------------------------
;
;-----------------------------------------------------------------------------
ddir = "/home/sunming/data5/cuixy/global_jets/data"
udir = ddir+"/ERA5_daily_u250_1979-2022_r2.5.nc"

rmmdir  = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
plotdir = "/home/sunming/data5/cuixy/global_jets/plot/Fig1_uclim"
;----------------------------------------------------------------------------
; draw u250.
;----------------------------------------------------------------------------
  f    = addfile(udir,"r")
  ut = f->u(:,0,:,:)
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)

  mon := mon(n1:n2)
  ut  := ut(n1:n2,:,:)
  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  u = dim_avg_n_Wrap(ut(nt,:,:),0)
  delete([/n1,n2,nt/])
; u250 resources.
;--plot
  res1 = True
  res1@gsnFrame     = False
  res1@gsnDraw      = False
  res1@gsnMaximize  = False
  res1@gsnAddCyclic = True
  res1@gsnMajorLatSpacing = 15
  res1@gsnMajorLonSpacing = 30

  res1@gsnLeftStringOrthogonalPosF  = 0.012
  res1@gsnRightStringOrthogonalPosF = 0.012
  res1@gsnLeftStringFontHeightF = 0.015
  res1@gsnRightStringFontHeightF = 0.015
  res1@tmXBLabelFontHeightF  = 0.012
  res1@tmYLLabelFontHeightF  = 0.012

  res1@vpXF     = 0.15
  res1@vpYF     = 0.9
  res1@vpWidthF = 0.7
  res1@vpHeightF= 0.4

  res1@mpFillOn = True
  res1@mpLandFillColor = "white"
  res1@mpCenterLonF      = 180.
  res1@mpMinLatF         = 0

  cmap1  = read_colormap_file("MPL_jet")
  cmap2 = cmap1(49::,:)
  res1@cnFillOn = True
  res1@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
  res1@cnFillPalette        = "WhiteBlueGreenYellowRed";cmap2
  res1@cnLinesOn            = False
  res1@cnLineLabelsOn       = False
  res1@cnLevelSelectionMode = "Manuallevels"

  res1@lbBoxLinesOn          = True
  res1@lbBoxSeparatorLinesOn = False
  res1@lbTitleString         = "~F21~[m s~S~-1~N~]"
  res1@lbTitlePosition       = "Bottom"
  res1@lbTitleFontHeightF    = 0.01
  res1@lbOrientation         = "Vertical"
  res1@pmLabelBarHeightF     = 0.2
  res1@pmLabelBarWidthF      = 0.1
  res1@pmLabelBarOrthogonalPosF = -0.01

  res1@tmXBLabelFont = 21
  res1@tmYRLabelFont = 21

  res1@gsnLeftString    = "~F21~(a) 250hPa U"
  res1@gsnRightString   = "~F21~DJF"

  wks = gsn_open_wks("eps",plotdir)
  plot  = new(1,graphic)
  plot  = gsn_csm_contour_map(wks, u, res1)

  respoly = True
  respoly@gsLineColor      = "black"
  respoly@gsLineThicknessF = 3.0

  txres               = True 
  txres@txFontHeightF = 0.012
  txres@txFont = 22
;; EAJS
  lonmin = 110
  lonmax = 180
  latmin = 15
  latmax = 50
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum1 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)
  tex1 = gsn_add_text(wks, plot, "EAJS", 125, 55, txres)
;; MEJS
  lonmin = 10
  lonmax = 70
  latmin = 10
  latmax = 45
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum2 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)
  txt2 = gsn_add_text(wks, plot, "MEJS", 35, 50, txres)
  ;; NAJS
  lonmin = -110
  lonmax = -40
  latmin =  15
  latmax =  60
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum3 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)
  txt3 = gsn_add_text(wks, plot, "NAJS", -45, 65, txres)

  draw(plot)

  ddir = "/home/sunming/data5/cuixy/global_jets/data"
  udir = ddir+"/"+"u250_daily_anom_1979-2022_10-90_r2.5.nc"

  mlats = 25
  mlatn = 35
  mlonw = 50
  mlone = 80

  elats = 27
  elatn = 38
  elonw = 170
  elone = 200

  nlats =  45
  nlatn =  55
  nlonw = -25
  nlone =  5
  undef("ut")
  f    = addfile(udir,"r")
  ut  = f->u
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  undef("mon")
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)

  mon := mon(n1:n2)
  ut  := ut(n1:n2,:,:)
  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  ut  := ut(nt,:,:)
  delete([/n1,n2,nt/])

  ustd = dim_stddev_n_Wrap(ut,0)
; u250 resources.
;--plot
  res1 = True
  res1@gsnFrame     = False
  res1@gsnDraw      = False
  res1@gsnMaximize  = False
  res1@gsnAddCyclic = True
  res1@gsnMajorLatSpacing = 15
  res1@gsnMajorLonSpacing = 30

  res1@gsnLeftStringOrthogonalPosF  = 0.015
  res1@gsnRightStringOrthogonalPosF = 0.015
  res1@gsnLeftStringFontHeightF  = 0.015
  res1@gsnRightStringFontHeightF = 0.015
  res1@tmXBLabelFontHeightF  = 0.012
  res1@tmYLLabelFontHeightF  = 0.012

  res1@vpXF     = 0.15
  res1@vpYF     = 0.63
  res1@vpWidthF = 0.7
  res1@vpHeightF= 0.4

  res1@mpFillOn = True
  res1@mpLandFillColor = "white"
  res1@mpCenterLonF      = -160
  res1@mpMinLatF         = 0

  cmap1  = read_colormap_file("MPL_jet")
  cmap2 = cmap1(49::,:)
  res1@cnFillOn = True
  res1@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
  res1@cnFillPalette        = "WhiteBlueGreenYellowRed";cmap2
  res1@cnLinesOn            = False
  res1@cnLineLabelsOn       = False
  res1@cnLevelSelectionMode = "Manuallevels"

  res1@lbBoxLinesOn          = True
  res1@lbBoxSeparatorLinesOn = False
  res1@lbTitleString         = "~F21~[m s~S~-1~N~]"
  res1@lbTitlePosition       = "Bottom"
  res1@lbTitleFontHeightF    = 0.01
  res1@lbOrientation         = "Vertical"
  res1@pmLabelBarHeightF     = 0.2
  res1@pmLabelBarWidthF      = 0.1
  res1@pmLabelBarOrthogonalPosF = -0.01

  res1@tmXBLabelFont = 21
  res1@tmYRLabelFont = 21

  res1@gsnLeftString    = "~F21~(b) ~F21~Filtered_250hPa_U_std"
  res1@gsnRightString   = "~F21~DJF"

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
  resc@cnLineColor        = "magenta"
  resc@cnLineThicknessF   = 3.0
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels      = (/30.,45.,60.,75./)

  plotc  = new(1,graphic)
  plotu = plot
  plotc = gsn_csm_contour_map(wks, ustd, res1)
  plotu = gsn_csm_contour(wks,u,resc)
  overlay(plotc,plotu)

  respoly = True
  respoly@gsLineColor      = "black"
  respoly@gsLineThicknessF = 3.0
;; EAJS
  lonmin = elonw
  lonmax = elone
  latmin = elats
  latmax = elatn
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum5 = gsn_add_polyline(wks, plotc, polyx, polyy, respoly)
;; MEJS
  lonmin = mlonw
  lonmax = mlone
  latmin = mlats
  latmax = mlatn
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum6 = gsn_add_polyline(wks, plotc, polyx, polyy, respoly)
;; NAJS
  lonmin = nlonw
  lonmax = nlone
  latmin = nlats
  latmax = nlatn
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum7 = gsn_add_polyline(wks, plotc, polyx, polyy, respoly)

  draw(plotc)

  ;i = 0
  ;do while(i.eq.0)
  ;  i = 0
  ;end do