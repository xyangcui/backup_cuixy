;; draw ratio.
  ddir  = "/home/sunming/data5/cuixy/global_jets/data"
  ucdir = ddir+"/ERA5_daily_u250_1979-2022_r2.5.nc"
  u2dir = ddir+"/"+"u250_daily_anom_1979-2022_10-90_r2.5.nc"   ;;10-100
  u3dir = ddir+"/"+"Reconstruct_u250_anom_30-60.nc"            ;;MJO-related                 
  pdir  = "/home/sunming/data5/cuixy/global_jets/plot/Fig3_ratio"
  wks = gsn_open_wks("eps",pdir)
  pi = atan(1.0)*4.0

  mlats = 25
  mlatn = 35
  mlonw = 50
  mlone = 80

  elats = 25
  elatn = 37
  elonw = 170
  elone = 200

  ntlats =  23
  ntlatn =  35
  ntlonw = -100
  ntlone = -60

  nlats =  45
  nlatn =  55
  nlonw = -25
  nlone =  5

  ax = 0.15
  ay = 0.95
  bx = 0.15
  by = 0.64
  cx = 0.15
  cy = 0.41
;-----------------------------------------------------------------------------
; Main code.
;-----------------------------------------------------------------------------
  fuc = addfile(ucdir,"r")
  uc  = fuc->u(:,0,:,:)
;; uwind
  fu2 = addfile(u2dir,"r")
  fu3 = addfile(u3dir,"r")

  u2  = fu2->u    
  u3  = fu3->u  ;reconstruct

  time = fu2->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)

  mon := mon(n1:n2)
  u2  := u2(n1:n2,:,:)
  uc  := uc(n1:n2,:,:)

  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  u2 := u2(nt,:,:)
  uc := uc(nt,:,:)

;; MEJS
  ut1  = u2(:,{mlats:mlatn},{mlonw:mlone})
  ut3  = u3(:,{mlats:mlatn},{mlonw:mlone})
  us3  = dim_stddev_n_Wrap(ut3, 0)
  us2  = dim_stddev_n_Wrap(ut1, 0)

  clat = cos(ut1&lat*pi/180.)
  ume2 = wgt_areaave_Wrap(us2, clat, 1.0, 1)
  ume3 = wgt_areaave_Wrap(us3, clat, 1.0, 1)

  delete([/ut1,ut3,us2,us3,clat/])

;;EAJS
  ut1  = u2(:,{elats:elatn},{elonw:elone})
  ut3  = u3(:,{elats:elatn},{elonw:elone})
  us3  = dim_stddev_n_Wrap(ut3, 0)
  us2  = dim_stddev_n_Wrap(ut1, 0)

  clat = cos(ut1&lat*pi/180.)
  uea2 = wgt_areaave_Wrap(us2, clat, 1.0, 1)
  uea3 = wgt_areaave_Wrap(us3, clat, 1.0, 1)
  
  delete([/ut1,ut3,us2,us3,clat/])

; NAJS
  utt2 = lonFlip(u2)
  utt3 = lonFlip(u3)

  ut1  = utt2(:,{ntlats:ntlatn},{ntlonw:ntlone})
  ut3  = utt3(:,{ntlats:ntlatn},{ntlonw:ntlone})
  us3  = dim_stddev_n_Wrap(ut3, 0)
  us2  = dim_stddev_n_Wrap(ut1, 0)

  clat = cos(ut1&lat*pi/180.)
  una2 = wgt_areaave_Wrap(us2, clat, 1.0, 1)
  una3 = wgt_areaave_Wrap(us3, clat, 1.0, 1)
  
  delete([/ut1,ut3,us2,us3,clat/])

;; std. ratio.
  mej = ume3/ume2
  print(mej)

  eaj = uea3/uea2
  print(eaj)

  naj = una3/una2
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

  res = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False
  res@gsnAddCyclic = True
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30

  res@gsnLeftStringOrthogonalPosF  = 0.0
  res@gsnRightStringOrthogonalPosF = 0.0
  res@gsnLeftStringFontHeightF  = 0.016
  res@gsnRightStringFontHeightF = 0.016
  res@tmXBLabelFontHeightF  = 0.012
  res@tmYLLabelFontHeightF  = 0.012

  res@mpFillOn = True
  res@mpLandFillColor = "white"
  res@mpCenterLonF      = -160
  res@mpMinLatF         = lats
  res@mpMaxLatF         = latn

  cmap1  = read_colormap_file("MPL_jet")
  cmap2 = cmap1(49::,:)
  res@cnFillOn = True
  res@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
  res@cnFillPalette        = "WhiteBlueGreenYellowRed";cmap2
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "Manuallevels"
  ;res1@cnMinLevelValF       =  0.
  ;res1@cnMaxLevelValF       =  16.
  ;res1@cnLevelSpacingF      =  2.0

  res@lbBoxLinesOn          = True
  res@lbBoxSeparatorLinesOn = False
  res@lbTitleString         = "~F21~[%]"
  res@lbTitlePosition       = "Right"
  res@lbTitleDirection      = "Across"
  res@lbTitleFontHeightF    = 0.015
  res@lbOrientation         = "Horizontal"
  res@pmLabelBarHeightF     = 0.05
  res@pmLabelBarWidthF      = 0.6
  res@pmLabelBarOrthogonalPosF = 0.17

  res@tmXBLabelFont = 21
  res@tmYRLabelFont = 21

  res@gsnLeftString    = "~F21~(a) ~F21~Ratio"
  ;res@gsnRightString   = "~F21~PCC:"+sprintf("%5.2f", pcc)

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

  res@vpYF = ay
  res@vpXF = ax
  res@vpWidthF = 0.7
  res@vpHeightF= 0.4
  plot = gsn_csm_contour_map(wks, ust, res)

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
;-----------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/global_jets/data"
  ucdir = ddir+"/ERA5_daily_u250_1979-2022_r2.5.nc"
  sdir  = ddir+"/monte_carlo"+"/MLR_sig.nc"
;;----------------------------------------------------------
function siglevel(ovalue,mvalue)
begin
  rv = ovalue
  lati = rv&lat
  loni = rv&lon
  do i = 0, dimsizes(lati)-1
    do j = 0, dimsizes(loni)-1

      if(ovalue(i,j).gt.mvalue(i,j)) then
        rv(i,j) =  2.
      else
        rv(i,j) = -2.
      end if

    end do
  end do
  return rv
end
;-----------------------------------------------------
; Main code.
;-----------------------------------------------------
;; sig value.
  fs = addfile(sdir,"r")
  R = fs->R2
  R95 = R(474,:,:)
;contour plot
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
  resc@cnLineColor        = "magenta";magenta
  resc@cnLineThicknessF   = 2.5
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels = (/30.,45.,60.,75./)

  f1 = addfile(ddir+"/MLRegression.nc","r")
  coAmp1 = f1->coAmp({0:60},:)
  R2 = f1->R2({0:60},:)
;; draw plot  
;;res.
  undef("res")
  res = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False
  res@gsnAddCyclic = True
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30

  res@gsnLeftStringOrthogonalPosF  = 0.012
  res@gsnRightStringOrthogonalPosF = 0.012
  res@gsnLeftStringFontHeightF  = 0.016
  res@gsnRightStringFontHeightF = 0.016
  res@tmXBLabelFontHeightF  = 0.012
  res@tmYLLabelFontHeightF  = 0.012

  res@mpFillOn = True
  res@mpLandFillColor = "white"
  res@mpCenterLonF    = -160
  res@mpMinLatF       = 0
  res@mpMaxLatF       = 60

  res@cnFillOn = True
  res@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
  res@cnFillPalette        = "MPL_Blues";
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "Manuallevels"
  res@cnMinLevelValF       =  0;0.
  res@cnMaxLevelValF       =  3;3.5
  res@cnLevelSpacingF      =  0.25;0.25

  res@lbLabelBarOn = True
  res@lbBoxLinesOn          = True
  res@lbBoxSeparatorLinesOn = False
  res@lbTitleString         = "~F21~[m s~S~-1~N~]"
  res@lbTitlePosition       = "Right"
  res@lbTitleDirection      = "Across"
  res@lbTitleFontHeightF    = 0.015
  res@lbOrientation         = "Horizontal"
  res@pmLabelBarHeightF     = 0.05
  res@pmLabelBarWidthF      = 0.6
  res@pmLabelBarOrthogonalPosF = 0.22

  res@tmXBLabelFont = 21
  res@tmYRLabelFont = 21
  res@gsnRightString = ""

  sres   = True 
  sres@cnLineDrawOrder    = "PostDraw"                    
  sres@gsnDraw            = False                          
  sres@gsnFrame           = False                        
  sres@cnLineLabelsOn     = False                              
  sres@cnLinesOn          = False                           
  sres@cnInfoLabelOn      = False
  sres@cnFillOn           = False                       
  sres@cnInfoLabelOn      = False
  sres@cnFillScaleF       = 0.5
  sres@cnFillDotSizeF     = 0.003
  sres@gsnLeftString      =" "
  sres@gsnRightString     =" "
  sres@cnMonoFillColor    = True
  sres@cnFillColor        = "gray"

  opt = True
  opt@gsnShadeFillType   = "pattern"
  opt@gsnShadeHigh       = 17
  opt@gsnShadeFillScaleF = 1.5

  undef("u")
  u = dim_avg_n_Wrap(uc,0)

  res@vpYF = by
  res@vpXF = bx
  res@vpWidthF = 0.7
  res@vpHeightF= 0.3
  res@gsnLeftString  = "~F21~(b) coAmp"
  plot1  = gsn_csm_contour_map(wks, coAmp1, res) 
  plotc1 = gsn_csm_contour(wks, u, resc)

  xt = siglevel(R2,R95)
  dum1 = gsn_csm_contour(wks,xt,sres)
  dum1 = gsn_contour_shade(dum1,1.9,-999.,opt)

  overlay(plot1,dum1)
  overlay(plot1,plotc1)

  draw(plot1)  

;-----------------------------------------------------
; draw coPhase for All, EI, LA and NE.
; MJO phases color table
; colors = (/"indianred3","darkorange2","darkgoldenrod2","lightgoldenrod2",\
; "lightgoldenrodyellow","darkseagreen3","darkturquoise","deepskyblue3"/)
;
; time: 20240528
;-----------------------------------------------------
  color1 = (/"indianred3","darkorange2","darkgoldenrod2",\
            "lightgoldenrod2","lightgoldenrodyellow",\
            "darkseagreen3","darkturquoise","deepskyblue3"/)
  color2 = (/"indianred3","darkorange2","darkorange2",\
          "lightgoldenrod2","lightgoldenrod2",\
          "darkseagreen3","darkseagreen3","indianred3"/)
  color = color1
  ddir = "/home/sunming/data5/cuixy/global_jets/data"
;-----------------------------------------------------
; Main code.
;-----------------------------------------------------
  f1 = addfile(ddir+"/MLRegression.nc","r")
  coPha1 = f1->coPha
;; draw plot  
;;res.
  undef("res")
  res = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False
  res@gsnAddCyclic = True
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30

  res@gsnLeftStringOrthogonalPosF  = 0.012
  res@gsnRightStringOrthogonalPosF = 0.012
  res@gsnLeftStringFontHeightF  = 0.016
  res@gsnRightStringFontHeightF = 0.016
  res@tmXBLabelFontHeightF  = 0.012
  res@tmYLLabelFontHeightF  = 0.012

  res@mpFillOn = True
  res@mpLandFillColor = "white"
  res@mpCenterLonF    = -160
  res@mpMinLatF       = 0
  res@mpMaxLatF       = 60

  res@cnFillOn = True
  res@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
  res@cnFillMode = "CellFill"
  res@cnFillPalette        = color
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = ispan(1,8,1)+.5
  res@lbLabelBarOn = False

  res@tmXBLabelFont = 21
  res@tmYRLabelFont = 21
  res@gsnRightString = ""

;; NE
  u = dim_avg_n_Wrap(uc,0)

  res@gsnLeftString  = "~F21~(c) coPhase"
  res@vpYF = cy
  res@vpXF = cx
  res@vpWidthF = 0.7
  res@vpHeightF= 0.3
  plot2  = gsn_csm_contour_map(wks, coPha1, res) 
  plotc2 = gsn_csm_contour(wks, u, resc)
  overlay(plot2,plotc2)

;; add label bar.
  lbres = True
  lbres@lbAutoManage = False
  lbres@lbFillColors = color
  lbres@lbMonoFillPattern  = True
  lbres@lbLabelFontHeightF = 0.015
  lbres@lbOrientation      = "Horizontal"
  lbres@lbPerimOn          = False
  lbres@lbTitleOn = True
  lbres@lbTitleFontHeightF = 0.015
  lbres@lbTitleString = "MJO phase"
  lbres@lbTitleDirection  = "Across"
  lbres@lbBoxMinorExtentF = 0.2
  lbres@lbTitlePosition = "Right"
  lbres@lbLabelOffsetF  = -0.14
  lbres@vpWidthF  = 0.5
  lbres@vpHeightF = 0.2
  lbres@lbBoxLinesOn = False
  lbres@lbFillOpacityF = 0.8

  lblb = tostring(ispan(1,8,1))
  lbid = gsn_create_labelbar(wks,8,lblb,lbres)
  amres = True
  amres@amJust = "TopCenter"
  amres@amParallelPosF   =  0.0
  amres@amOrthogonalPosF =  0.12
  annoid = gsn_add_annotation(plot2,lbid,amres)

  draw(plot2)  

  ;i = 0
  ;do while (i.eq.0)
  ;  i = 0
  ;end do