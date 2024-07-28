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

  color = color2

  ddir = "/home/sunming/data5/cuixy/ENSO-MJO-Jet"
  ;ucdir = ddir+"/data/u_DJF_mean_1979-2022_r2.5.nc"
  ucdir = ddir+"/data/u250_daily_mean_1979-2022.nc"
;-----------------------------------------------------
; Main code.
;-----------------------------------------------------
  fuc = addfile(ucdir,"r")
  uc  = fuc->u(:,0,:,:)

  time = uc&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)
  year := year(n1:n2)
  mon  := mon(n1:n2)
  uc   := uc(n1:n2,:,:)

  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  year := year(nt)
  mon  := mon(nt)
  uc   := uc(nt,:,:)

  ni = ind(mon.eq.12)
  year(ni) = year(ni)+1
;; nino years.
  ninodir = "/home/sunming/data5/cuixy/data/clim_indices/Nino_3.4_index.csv"

  delim  = ","
  ninof  = asciiread(ninodir,-1,"string")
  iyear  = str_get_field(ninof,1,delim)
  ievent = str_get_field(ninof,6,delim)

  iyear  := toint(iyear(1::))
  ievent := tofloat(ievent(1::))

;;EI
  nei = ind(ievent.eq.1)
  yei = iyear(nei)
;;LA
  nla = ind(ievent.eq.-1)
  yla = iyear(nla)
;;NE
  nne = ind(ievent.eq.0)
  yne = iyear(nne)

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

  f1 = addfile(ddir+"/data/MLRegression_NE.nc","r")
  coPha1 = f1->coPha

  f2 = addfile(ddir+"/data/MLRegression_EI.nc","r")
  coPha2 = f2->coPha

  f3 = addfile(ddir+"/data/MLRegression_LA.nc","r")
  coPha3 = f3->coPha

  f4 = addfile(ddir+"/data/MLRegression_all.nc","r")
  coPha4 = f4->coPha

;; draw plot  
;;res.
  res = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False
  res@gsnAddCyclic = True
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30

  res@gsnLeftStringOrthogonalPosF  = 0.012
  res@gsnRightStringOrthogonalPosF = 0.012
  res@gsnLeftStringFontHeightF  = 0.02
  res@gsnRightStringFontHeightF = 0.02
  res@tmXBLabelFontHeightF  = 0.012
  res@tmYLLabelFontHeightF  = 0.012

  res@mpFillOn = True
  res@mpLandFillColor = "white"
  res@mpCenterLonF    = 180.
  res@mpMinLatF       = 5
  res@mpMaxLatF       = 55
  res@mpMinLonF       = 0
  res@mpMaxLonF       = 360

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

  wks  = gsn_open_wks("x11","maponly")

;; NE
  ny = yne
  nc = ind(year.eq.ny(0))
  do i = 1, dimsizes(ny)-1
    ntemp = ind(year.eq.ny(i))
    nc := array_append_record(nc, ntemp, 0)
    delete(ntemp)
  end do
  u = dim_avg_n_Wrap(uc(nc,:,:),0)

  res@gsnLeftString  = "~F21~(a) NE"
  res@vpYF = 0.8
  plot1  = gsn_csm_contour_map(wks, coPha1, res) 
  plotc1 = gsn_csm_contour(wks, u, resc)
  overlay(plot1,plotc1)
  delete([/nc,u,ny/])

;; EI
  ny = yei
  nc = ind(year.eq.ny(0))
  do i = 1, dimsizes(ny)-1
    ntemp = ind(year.eq.ny(i))
    nc := array_append_record(nc, ntemp, 0)
    delete(ntemp)
  end do
  u = dim_avg_n_Wrap(uc(nc,:,:),0)

  res@gsnLeftString  = "~F21~(b) EI"
  res@vpYF = 0.6
  plot2  = gsn_csm_contour_map(wks, coPha2, res)
  plotc2 = gsn_csm_contour(wks, u, resc)
  overlay(plot2,plotc2)
  delete([/nc,u,ny/])

;; EI
  ny = yla
  nc = ind(year.eq.ny(0))
  do i = 1, dimsizes(ny)-1
    ntemp = ind(year.eq.ny(i))
    nc := array_append_record(nc, ntemp, 0)
    delete(ntemp)
  end do
  u = dim_avg_n_Wrap(uc(nc,:,:),0)

  res@gsnLeftString  = "~F21~(c) LA"
  res@vpYF = 0.4
  plot3  = gsn_csm_contour_map(wks, coPha3, res)
  plotc3 = gsn_csm_contour(wks, u, resc)
  overlay(plot3,plotc3)
  delete([/nc,u,ny/])
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
  lbres@lbTitlePosition = "Left"
  lbres@lbLabelOffsetF  = -0.14
  lbres@vpWidthF  = 0.5
  lbres@vpHeightF = 0.25
  lbres@lbBoxLinesOn = False
  lbres@lbFillOpacityF = 0.8

  lblb = tostring(ispan(1,8,1))
  lbid = gsn_create_labelbar(wks,8,lblb,lbres)
  amres = True
  amres@amJust = "TopCenter"
  amres@amParallelPosF   =  0.0
  amres@amOrthogonalPosF =  0.0
  annoid = gsn_add_annotation(plot3,lbid,amres)

  draw(plot1)  
  draw(plot2) 
  draw(plot3) 

  i = 0
  do while (i.eq.0)
    i = 0
  end do