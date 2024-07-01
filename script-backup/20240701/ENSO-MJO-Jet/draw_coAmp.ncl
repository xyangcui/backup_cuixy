;-----------------------------------------------------
; draw coPhase for All, EI, LA and NE.
; MJO phases color table
; colors = (/"indianred3","darkorange2","darkgoldenrod2","lightgoldenrod2",\
; "lightgoldenrodyellow","darkseagreen3","darkturquoise","deepskyblue3"/)
;
; time: 20240528
;-----------------------------------------------------
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
  coAmp1 = f1->coAmp
  coAmp1 = coAmp1;*1.33
  pc1 = f1->pc

  f2 = addfile(ddir+"/data/MLRegression_EI.nc","r")
  coAmp2 = f2->coAmp
  coAmp2 = coAmp2;*1.37
  pc2 = f2->pc

  f3 = addfile(ddir+"/data/MLRegression_LA.nc","r")
  coAmp3 = f3->coAmp
  coAmp3 = coAmp3;*1.4
  pc3 = f3->pc

  f4 = addfile(ddir+"/data/MLRegression_all.nc","r")
  coAmpa = f4->coAmp
  coAmpa = coAmpa;*1.4
  pc4 = f4->pc


  ;coAmp1 = pc1(1,:,:)
  ;coAmp2 = pc2(1,:,:)
  ;coAmp3 = pc3(1,:,:)
;; 160-200E 25-35N
  ;dife = pc2 - pc1
  ;copy_VarCoords(coAmp1,dife(0,:,:))
  ;dif1 = dim_avg_n(dife(:,{15:25},{270:300}),(/1,2/))
  ;print(dif1)

  coAmp4 = coAmp2 - coAmp1
  coAmp5 = coAmp3 - coAmp1
  coAmp4 = dife(0,:,:)
  coAmp5 = dife(1,:,:)
  copy_VarCoords(coAmp2, coAmp4)
  copy_VarCoords(coAmp2, coAmp5)
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
  res@gsnLeftStringFontHeightF  = 0.016
  res@gsnRightStringFontHeightF = 0.016
  res@tmXBLabelFontHeightF  = 0.012
  res@tmYLLabelFontHeightF  = 0.012

  res@mpFillOn = True
  res@mpLandFillColor = "white"
  res@mpCenterLonF    = 180.
  res@mpMinLatF       = 5
  res@mpMaxLatF       = 55
  res@mpMinLonF       = 0
  res@mpMaxLonF       = 325

  res@cnFillOn = True
  res@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
  res@cnFillPalette        = ;"MPL_coolwarm";MPL_Blues";
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "Manuallevels"
  res@cnMinLevelValF       =  -3.5;0.
  res@cnMaxLevelValF       =   3.5;3.5
  res@cnLevelSpacingF      =   0.5;0.25
  res@lbLabelBarOn = False
  res@pmLabelBarOrthogonalPosF = 0.25

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
  res@vpYF = 0.95
  plot1  = gsn_csm_contour_map(wks, coAmp1, res) 
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
  res@vpYF = 0.78
  plot2  = gsn_csm_contour_map(wks, coAmp2, res)
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
  res@lbLabelBarOn = True
  res@gsnLeftString  = "~F21~(c) LA"
  res@vpYF = 0.61
  plot3  = gsn_csm_contour_map(wks, coAmp3, res)
  plotc3 = gsn_csm_contour(wks, u, resc)
  overlay(plot3,plotc3)

;;draw difference.
  res@cnFillPalette        = "MPL_RdBu";"WhiteBlueGreenYellowRed"
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "Manuallevels"
  res@cnMinLevelValF       =  -1.4
  res@cnMaxLevelValF       =   1.4
  res@cnLevelSpacingF      =   0.2
  res@lbLabelBarOn = False
  res@gsnLeftString  = "~F21~(d) EI-NE"
  res@vpYF = 0.4
  plot4  = gsn_csm_contour_map(wks, coAmp4, res)
  plotc4 = gsn_csm_contour(wks, u, resc)
  overlay(plot4,plotc4)

  res@lbLabelBarOn = True
  res@gsnLeftString  = "~F21~(e) LA-NE"
  res@vpYF = 0.2
  plot5  = gsn_csm_contour_map(wks, coAmp5, res)
  plotc5 = gsn_csm_contour(wks, u, resc)
  overlay(plot5,plotc5)

  draw(plot1)  
  draw(plot2) 
  draw(plot3) 
  draw(plot4) 
  draw(plot5) 

  i = 0
  do while (i.eq.0)
    i = 0
  end do