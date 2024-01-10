load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
;-----------------------------------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/global_jets/data"
  udir = ddir+"/ERA5_daily_u250_1979-2022.nc"
  stddir = ddir+"/ERA5_u250_std.nc"
;; jet indices.
  eadir = ddir+"/jet_indices_EAJS.nc"
  medir = ddir+"/jet_indices_MEJS.nc"
  nadir = ddir+"/jet_indices_NAJS.nc"

  plotdir = "/home/sunming/data5/cuixy/global_jets/plot/std_jet_indices"

;-----------------------------------------------------------------------------
; Main code.
;-----------------------------------------------------------------------------
  f    = addfile(udir,"r")
  fs   = addfile(stddir,"r")
  ustd = fs->u(0,:,:)
  time = f->time
  TIME = cd_calendar(time, 0)
  mon  = toint(TIME(:,1))
  nt   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)

  u    = dim_avg_n_Wrap(f->u(nt,0,:,:),0)
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
  res1@gsnLeftStringFontHeightF = 0.012
  res1@gsnRightStringFontHeightF = 0.012
  res1@tmXBLabelFontHeightF  = 0.012
  res1@tmYLLabelFontHeightF  = 0.012

  res1@vpXF     = 0.15
  res1@vpYF     = 0.98
  res1@vpWidthF = 0.7
  res1@vpHeightF= 0.25

  res1@mpFillOn = True
  res1@mpLandFillColor = "white"
  res1@mpCenterLonF      = 60.
  res1@mpMinLatF         = 0

  cmap1  = read_colormap_file("MPL_jet")
  cmap2 = cmap1(49::,:)
  res1@cnFillOn = True
  res1@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
  res1@cnFillPalette        = "WhiteBlueGreenYellowRed";cmap2
  res1@cnLinesOn            = False
  res1@cnLineLabelsOn       = False
  res1@cnLevelSelectionMode = "Manuallevels"
  ;res1@cnMinLevelValF       =  5.
  ;res1@cnMaxLevelValF       =  70.
  ;res1@cnLevelSpacingF      =  1.

  res1@lbBoxLinesOn          = True
  res1@lbBoxSeparatorLinesOn = False
  ;res1@lbLabelStride         = 10
  res1@lbTitleString         = "~F21~(unit: m~S~2~N~/s~S~2~N~)"
  res1@lbTitlePosition       = "Bottom"
  res1@lbTitleFontHeightF    = 0.01
  res1@lbOrientation         = "Vertical"
  res1@pmLabelBarHeightF     = 0.2
  res1@pmLabelBarWidthF      = 0.1
  ;res1@pmLabelBarParallelPosF   = -0.003
  res1@pmLabelBarOrthogonalPosF = -0.01

  res1@tmXBLabelFont = 21
  res1@tmYRLabelFont = 21

  res1@gsnLeftString    = "~F21~(a) ~F21~250hPa U_std"
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
  resc@cnLevels      = (/25.,30.,40.,50.,60.,70./)

  wks = gsn_open_wks("x11",plotdir)
  plot  = new(1,graphic)
  plotu = plot
  plot  = gsn_csm_contour_map(wks, ustd, res1)
  plotu = gsn_csm_contour(wks,u,resc)
  overlay(plot,plotu)

  respoly = True
  respoly@gsLineColor      = "black"
  respoly@gsLineThicknessF = 3.0
;; EAJS
  lonmin = 110
  lonmax = 180
  latmin = 30
  latmax = 40
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum1 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)
;; MEJS
  lonmin = 10
  lonmax = 60
  latmin = 25
  latmax = 35
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum2 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)
;; NAJS
  lonmin = -60
  lonmax = 0
  latmin = 30
  latmax = 60
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum3 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)

  draw(plot)
  delete([/u,time,TIME,mon,nt,cmap1,cmap2/])

;;----------------------------------------------------------------------
;; draw spec.
;-----------------------------------------------------------------------

  i = 0
  do while(i.eq.0)
    i = 0
  end do
