load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
;-----------------------------------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/global_jets/data"
  ucdir = ddir+"/ERA5_daily_u250_1979-2022_r2.5.nc"
  udir = ddir+"/"+"ERA5_daily_anom_u250_1979-2022_2.5x2.5.nc"  ;;20-100
  u2dir = ddir+"/"+"u250_daily_anom_1979-2022_10-90_r2.5.nc"   ;;10-100
  ;u2dir = ddir+"/"+"u250_daily_anom_1979-2022_30-90_r2.5.nc"   ;;10-100
  u3dir = ddir+"/"+"Reconstruct_u250_anom_30-60.nc"            ;;MJO-related                 
  plotdir = "/home/sunming/data5/cuixy/global_jets/plot/Fig1_spec"

  pi = atan(1.0)*4.0

  mlats = 25
  mlatn = 35
  mlonw = 50
  mlone = 80

  elats = 25;25
  elatn = 37;45
  elonw = 170;165
  elone = 200;200

  ntlats =  23
  ntlatn =  35
  ntlonw = -100
  ntlone = -60

  nlats =  45
  nlatn =  55
  nlonw = -25
  nlone =  5
;-----------------------------------------------------------------------------
; Main code.
;-----------------------------------------------------------------------------
  fuc = addfile(ucdir,"r")
  uc  = fuc->u(:,0,:,:)
;; uwind
  fu  = addfile(udir,"r")
  fu2 = addfile(u2dir,"r")
  fu3 = addfile(u3dir,"r")

  f1 = addfile(ddir+"/MLRegression.nc","r")
  u4 = f1->coAmp

  u   = fu->u   ;20-100
  u2  = fu2->u  ;10-100   
  u3  = fu3->u  ;reconstruct

  time = fu->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)

  mon := mon(n1:n2)
  u   := u(n1:n2,:,:)
  u2  := u2(n1:n2,:,:)
  uc  := uc(n1:n2,:,:)

  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  u  := u(nt,:,:)
  u2 := u2(nt,:,:)
  uc := uc(nt,:,:)

;; MEJS
  ut   = u(:,{mlats:mlatn},{mlonw:mlone})
  ut1  = u2(:,{mlats:mlatn},{mlonw:mlone})
  ut3  = u3(:,{mlats:mlatn},{mlonw:mlone})
  us3  = dim_stddev_n_Wrap(ut3, 0)
  us   = dim_stddev_n_Wrap(ut, 0)
  us2  = dim_stddev_n_Wrap(ut1, 0)

  clat = cos(ut&lat*pi/180.)
  ume  = wgt_areaave_Wrap(us,  clat, 1.0, 1)
  ume2 = wgt_areaave_Wrap(us2, clat, 1.0, 1)
  ume3 = wgt_areaave_Wrap(us3, clat, 1.0, 1)

  us4  = u4({mlats:mlatn},{mlonw:mlone})
  clat1 = cos(us4&lat*pi/180.)
  ume4  = wgt_areaave_Wrap(us4, clat1, 1.0, 1)

  delete([/ut,ut1,ut3,us,us2,us3,us4,clat1,clat/])

;;EAJS
  ut   = u(:,{elats:elatn},{elonw:elone})
  ut1  = u2(:,{elats:elatn},{elonw:elone})
  ut3  = u3(:,{elats:elatn},{elonw:elone})
  us3  = dim_stddev_n_Wrap(ut3, 0)
  us   = dim_stddev_n_Wrap(ut, 0)
  us2  = dim_stddev_n_Wrap(ut1, 0)

  clat = cos(ut&lat*pi/180.)
  uea  = wgt_areaave_Wrap(us,  clat, 1.0, 1)
  uea2 = wgt_areaave_Wrap(us2, clat, 1.0, 1)
  uea3 = wgt_areaave_Wrap(us3, clat, 1.0, 1)

  us4  = u4({elats:elatn},{elonw:elone})
  clat1 = cos(us4&lat*pi/180.)
  uea4  = wgt_areaave_Wrap(us4, clat1, 1.0, 1)
  
  delete([/ut,ut1,ut3,us,us2,us3,us4,clat1,clat/])

; NAJS
  utt1 = lonFlip(u)
  utt2 = lonFlip(u2)
  utt3 = lonFlip(u3)

  ut   = utt1(:,{ntlats:ntlatn},{ntlonw:ntlone})
  ut1  = utt2(:,{ntlats:ntlatn},{ntlonw:ntlone})
  ut3  = utt3(:,{ntlats:ntlatn},{ntlonw:ntlone})
  us3  = dim_stddev_n_Wrap(ut3, 0)
  us   = dim_stddev_n_Wrap(ut, 0)
  us2  = dim_stddev_n_Wrap(ut1, 0)

  clat = cos(ut&lat*pi/180.)
  una  = wgt_areaave_Wrap(us,  clat, 1.0, 1)
  una2 = wgt_areaave_Wrap(us2, clat, 1.0, 1)
  una3 = wgt_areaave_Wrap(us3, clat, 1.0, 1)

  u4 = lonFlip(u4)
  us4  = u4({ntlats:ntlatn},{ntlonw:ntlone})
  clat1 = cos(us4&lat*pi/180.)
  una4  = wgt_areaave_Wrap(us4, clat1, 1.0, 1)
  
  delete([/ut,ut1,ut3,us,us2,us3,us4,clat1,clat/])

;; std. ratio.
  mej = ume3/ume2;*1.7
  print(mej)

  eaj = uea3/uea2;*1.7
  print(eaj)

  naj = una3/una2;*1.7
  print(naj)

  lats = 0
  latn = 90
  lonw = 0
  lone = 360
;;recalcute

  ust1 = dim_stddev_n_Wrap(u3(:,{lats:latn},{lonw:lone}),0)
  ust2 = dim_stddev_n_Wrap(u2(:,{lats:latn},{lonw:lone}),0)

  undef("ust")
  ust = ust1/ust2*100
  copy_VarCoords(ust1, ust)

  ust3 = ust1({20:45},{30:210})
  ust4 = ust2({20:45},{30:210})
  clat = cos(0.01745329*ust4&lat)
  pcc = tofloat( pattern_cor(ust3, ust4, clat, 0) )
  print(pcc)

  res1 = True
  res1@gsnFrame     = False
  res1@gsnDraw      = False
  res1@gsnMaximize  = False
  res1@gsnAddCyclic = True
  res1@gsnMajorLatSpacing = 15
  res1@gsnMajorLonSpacing = 30

  res1@gsnLeftStringOrthogonalPosF  = 0.0
  res1@gsnRightStringOrthogonalPosF = 0.0
  res1@gsnLeftStringFontHeightF  = 0.015
  res1@gsnRightStringFontHeightF = 0.015
  res1@tmXBLabelFontHeightF  = 0.012
  res1@tmYLLabelFontHeightF  = 0.012

  res1@vpXF     = 0.15
  res1@vpYF     = 0.9
  res1@vpWidthF = 0.7
  res1@vpHeightF= 0.4

  res1@mpFillOn = True
  res1@mpLandFillColor = "white"
  res1@mpCenterLonF      = -160
  res1@mpMinLatF         = lats
  res1@mpMaxLatF         = latn

  cmap1  = read_colormap_file("MPL_jet")
  cmap2 = cmap1(49::,:)
  res1@cnFillOn = True
  res1@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
  res1@cnFillPalette        = "WhiteBlueGreenYellowRed";cmap2
  res1@cnLinesOn            = False
  res1@cnLineLabelsOn       = False
  res1@cnLevelSelectionMode = "Manuallevels"
  ;res1@cnMinLevelValF       =  0.
  ;res1@cnMaxLevelValF       =  16.
  ;res1@cnLevelSpacingF      =  2.0

  res1@lbBoxLinesOn          = True
  res1@lbBoxSeparatorLinesOn = False
  res1@lbTitleString         = "~F21~[%]"
  res1@lbTitlePosition       = "Bottom"
  res1@lbTitleFontHeightF    = 0.01
  res1@lbOrientation         = "Vertical"
  res1@pmLabelBarHeightF     = 0.2
  res1@pmLabelBarWidthF      = 0.1
  res1@pmLabelBarOrthogonalPosF = -0.01

  res1@tmXBLabelFont = 21
  res1@tmYRLabelFont = 21

  res1@gsnLeftString    = "~F21~(a) ~F21~Ratio"
  res1@gsnRightString   = "~F21~PCC:"+sprintf("%5.2f", pcc)

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

  wks = gsn_open_wks("x11",plotdir)
  plot  = new(1,graphic)
  plotu = plot
  plot  = gsn_csm_contour_map(wks, ust, res1)

  undef("u")
  u = dim_avg_n_Wrap(uc,0)
  plotu = gsn_csm_contour(wks, u, resc)
  overlay(plot,plotu)

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
  dum1 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)
;; MEJS
  lonmin = mlonw
  lonmax = mlone
  latmin = mlats
  latmax = mlatn
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum2 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)
;; NAJS
  lonmin = nlonw
  lonmax = nlone
  latmin = nlats
  latmax = nlatn
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum3 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)

  respoly@gsLineColor  = "red"

  lonmin = ntlonw
  lonmax = ntlone
  latmin = ntlats
  latmax = ntlatn
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum4 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)

  draw(plot)

  i = 0
  do while (i.eq.0)
    i = 0
  end do