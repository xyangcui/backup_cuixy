;-----------------------------------------------------------------
;; MEJS
;; draw T tendency and T anomaly in Line chart.
;; draw OLR, udtdx and omega.
;; time: 2024-05-02
;; two part: (a) Line chart.
;-----------------------------------------------------------------
  rad   = 4.*atan(1.)/180. 
  tday  = 24.*60.*60.
  nlatn =  50
  nlats =  35
  nlonw = -105
  nlone = -75

  lats =  20
  latn =  60
  lonw = -120
  lone = -20

  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJS_temp_budget"
  tdir = ddir+"/"+"anom_dtdt_r2.5.nc"
  vdir = ddir+"/"+"anom_vdtdy_r2.5.nc"
  ttdir = "/home/sunming/data5/cuixy/global_jets/data/ERA5_daily_anom_tvert_1979-2022_2.5x2.5.nc"

  t1dir = ddir+"/"+"ERA5_daily_low_tvert_1979-2022_r2.5.nc"
  ucdir = ddir+"/"+"ERA5_daily_iso_uvert_1979-2022_r2.5.nc"
  vcdir = ddir+"/"+"ERA5_daily_iso_vvert_1979-2022_r2.5.nc"

  pdir   = "/home/sunming/data5/cuixy/global_jets/plot/Fig11_Linechart"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
  label = (/"day -5","day -2","day 1","day 4"/)
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
;; read in T ten and anom.
  ft  = addfile(tdir,"r")
  ftt = addfile(ttdir,"r")
  t  = ft->dtdt
  ta = ftt->t
  ta = lonFlip(ta)

  ft1 = addfile(t1dir,"r")
  t1  = ft1->t(:,{lats:latn},{lonw:lone})

  fv = addfile(vdir,"r")
  v = fv->vdtdy(:,{lats:latn},{lonw:lone})*tday*(-1)

  fu1 = addfile(ucdir,"r")
  fv1 = addfile(vcdir,"r")
  ua = fu1->u(:,{lats:latn},{lonw:lone})
  va = fv1->v(:,{lats:latn},{lonw:lone})

;--------------------------------------------------------------=
;; MJO phases.
  time = ft->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
  delim  = ","
  rmmf   = asciiread(rmmdir,-1,"string")
  iyear  = str_get_field(rmmf,1,delim)
  imon   = str_get_field(rmmf,2,delim)
  iday   = str_get_field(rmmf,3,delim)
  phaset = str_get_field(rmmf,6,delim)
  ampt   = str_get_field(rmmf,7,delim)

  iyear  := toint(iyear(1::))
  imon   := toint(imon(1::))
  iday   := toint(iday(1::))
  phaset := toint(phaset(1::))
  ampt   := tofloat(ampt(1::))

  amp   = new((/dimsizes(time)/),float)
  phase = new((/dimsizes(time)/),integer)
  amp   = 0.
  phase = 0

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
    amp(i)   = ampt(nt)
    phase(i) = phaset(nt)
    delete(nt)
  end do

  nn = ispan(0, dimsizes(time)-1, 1)
  amp!0   = "timeind"
  phase!0 = "timeind"
  amp&timeind   = nn
  phase&timeind = nn
  delete([/nn,ampt,phaset,iyear,imon,iday/])

  print("complete RMM selling "+ systemfunc("date"))

  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.2.and.day.eq.28)

  amp   := amp(n1:n2)
  phase := phase(n1:n2)
  mon   := mon(n1:n2)
  time  := time(n1:n2)

  nt     = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  amp   := amp(nt)
  phase := phase(nt)

  delete(nt)
  ind_DJF = amp&timeind   ; indice of DJF between 1979-2022 DJF.
;---------------------------------------------------------------
;; store lead-lag values for plotting. 
  tm = t(:,{nlats:nlatn},{nlonw:nlone})
  clatjs = cos(tm&lat*rad)
  tajs = wgt_areaave_Wrap(ta(:,{nlats:nlatn},{nlonw:nlone})\
                         , clatjs, 1.0, 1)
  tjs  = wgt_areaave_Wrap(t(:,{nlats:nlatn},{nlonw:nlone})\
                         , clatjs, 1.0, 1)
;;Phase 4.
  nt = ind(phase.eq.4.and.amp.gt.1.)
  nind = ind_DJF(nt)

  llp4 = new((/41,2/),float)
  llp4 = 0.

  do j = 0, 40
    ntt   = -10 + j
    nindd = nind + ntt
    llp4(j,0) = dim_avg_n(tjs(nindd), 0)
    llp4(j,1) = dim_avg_n(tajs(nindd), 0)
    delete([/ntt,nindd/])
  end do
  delete([/nt,nind/])
;;phase 8
  nt = ind(phase.eq.8.and.amp.gt.1.)
  nind = ind_DJF(nt)

  llp8 = new((/41,2/),float)
  llp8 = 0.

  do j = 0, 40
    ntt   = -10 + j
    nindd = nind + ntt
    llp8(j,0) = dim_avg_n(tjs(nindd), 0)
    llp8(j,1) = dim_avg_n(tajs(nindd), 0)
    delete([/ntt,nindd/])
  end do
  delete([/nt,nind/])
  print("calc data done. "+ systemfunc("date"))
;------------------------------------------------------
;; plot format: Line chart. about 5 lines. RGBA color.
  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@gsnRightStringFontHeightF  = 0.015
  res@gsnLeftStringFontHeightF   = 0.015
  res@gsnLeftStringOrthogonalPosF = -0.02
  res@gsnRightStringOrthogonalPosF = -0.02

  res@tmXTOn = False
  res@tmYROn = False

  res@tiXAxisString  = "~F21~Lag days"
  res@tiYAxisString  = "~F21~Temp tendency (K/day)"
  res@tmXBMode   = "Explicit"
  res@tmXBValues = (/ 0. , 5., 10., 15., 20., 25., 30., 35., 40. /)
  res@tmXBLabels = (/"-10","-5","0","5","10","15","20","25","30" /)
  res@tmXBLabelFontHeightF = 0.012
  res@tiXAxisFontHeightF = 0.012

  res@tiYAxisFontHeightF = 0.012

  res@tmYLMode        = "Manual"
  res@tmYLMinorPerMajor = 0
  res@tmYLTickSpacingF  = 0.3
  res@tmYLLabelFontHeightF = 0.012

  res@tmYLMajorLengthF   = -0.01
  res@gsnYRefLine = 0.
  res@gsnYRefLineDashPattern = 0

  res@xyLineThicknessF = 4.0

  color = (/"gold","red3"/)  ;"royalblue4"
  xlabel = ispan(0,40,1)
;----------------------------------------------------------------
;; Line plot
  wks = gsn_open_wks("x11", pdir)
  ;drawNDCGrid(wks)
  res@vpWidthF  = 0.38
  res@vpHeightF = 0.18
  res@vpXF = 0.1
  res@vpYF = 0.97
  res@trYMaxF  =  0.08
  res@trYMinF  = -0.08
  res@trXMaxF  =  20
  res@tmYLTickSpacingF  = 0.02
  res@tmYLPrecision = 1
  res@gsnLeftString  = "~F21~(a) (~F34~6~F21~T/~F34~6~F21~t~F21~)' & T' "
  res@gsnRightString = "~F21~ave: "+nlats+"-"+nlatn+"N, "+nlonw+"-"+nlone+"E"
  res@xyLineColor = color(0)
  res@xyDashPattern = 1

  resR = True
  resR@trYMaxF  = 0.8
  resR@trYMinF  = -0.8
  resR@tmYLTickSpacingF  = 0.2
  resR@xyLineThicknessF = 4.0
  resR@xyLineColor = color(0)
  resR@xyDashPattern = 0
  resR@tmYROn = False
  dum = gsn_csm_xy2(wks, xlabel, llp4(:,0),llp4(:,1), res,resR)

  restx = True
  restx@txFontHeightF = 0.02

  gsn_text_ndc(wks, "~F21~P4", 0.4, 0.95, restx)

  draw(dum)

  res@vpXF = 0.52
  res@vpYF = 0.97
  res@tmYLOn = False
  res@tiYAxisString  = ""
  res@gsnLeftString  = "~F21~(b) (~F34~6~F21~T/~F34~6~F21~t~F21~)' & T' "
  resR@tmYROn = True
  resR@tiYAxisString  = "~F21~Temp anomaly (K)"
  dum2 = gsn_csm_xy2(wks, xlabel, llp8(:,0),llp8(:,1), res,resR)

  gsn_text_ndc(wks, "~F21~P8", 0.8, 0.95, restx)

  draw(dum2)
;; for 2-D plot.
  nt = ind(phase.eq.4.and.amp.gt.1.)
  nind = ind_DJF(nt)
  llt  = new((/4,dimsizes(t1&lat),dimsizes(t1&lon)/),float)
  llt  = 0.
  llv  = llt
  llua = llt
  llva = llt
  printVarSummary(llt)
  printVarSummary(v)
  do j = 0, 3
    ntt   = -5 + j*3
    nindd = nind + ntt
    llt(j,:,:)  = dim_avg_n(t1(nindd,:,:), 0)
    llv(j,:,:)  = dim_avg_n(v(nindd,:,:), 0)
    llua(j,:,:) = dim_avg_n(ua(nindd,:,:), 0)
    llva(j,:,:) = dim_avg_n(va(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do

  copy_VarCoords(t1(0,:,:), llt(0,:,:))
  copy_VarCoords(t1(0,:,:), llua(0,:,:))
  copy_VarCoords(t1(0,:,:), llva(0,:,:))
  copy_VarCoords(t1(0,:,:), llv(0,:,:))
  delete([/nt,nind/])

;;Phase 8
  nt = ind(phase.eq.8.and.amp.gt.1.)
  nind = ind_DJF(nt)
  lltt  = new((/4,dimsizes(t1&lat),dimsizes(t1&lon)/),float)
  lltt  = 0.
  llvt  = lltt
  lluat = lltt
  llvat = lltt
  do j = 0, 3
    ntt   = -5 + j*3
    nindd = nind + ntt
    lltt(j,:,:)  = dim_avg_n(t1(nindd,:,:), 0)
    llvt(j,:,:)  = dim_avg_n(v(nindd,:,:), 0)
    lluat(j,:,:) = dim_avg_n(ua(nindd,:,:), 0)
    llvat(j,:,:) = dim_avg_n(va(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do

  copy_VarCoords(t1(0,:,:), lltt(0,:,:))
  copy_VarCoords(t1(0,:,:), lluat(0,:,:))
  copy_VarCoords(t1(0,:,:), llvat(0,:,:))
  copy_VarCoords(t1(0,:,:), llvt(0,:,:))
  delete([/nt,nind/])
;shading for T.
  cmap = read_colormap_file("MPL_RdBu")

  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = False
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 20
  rest@gsnLeftStringFontHeightF  = 0.035
  rest@gsnRightStringFontHeightF = 0.035

  rest@tmXBLabelFont = 21
  rest@tmYLLabelFont = 21

  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF   = lats
  rest@mpMaxLatF   = latn
  rest@mpMinLonF   = lonw-15
  rest@mpMaxLonF   = lone+15
  rest@mpLandFillColor = "white"

  rest@tmXBLabelFontHeightF = 0.025
  rest@tmYLLabelFontHeightF = 0.025

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -1.
  rest@cnMaxLevelValF         =  1.
  rest@cnLevelSpacingF        =  0.2
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = False

; contour phi plot for dtdt.
  resc = True
  resc@gsnDraw      = False
  resc@gsnMaximize  = False
  resc@gsnFrame     = False
  resc@gsnAddCyclic = False
  resc@gsnContourNegLineDashPattern = 1
  resc@cnMonoLineColor    = True
  resc@cnLinesOn          = True
  resc@cnLineLabelsOn     = False
  resc@cnLineColor        = "mediumorchid4"
  resc@cnLineThicknessF   = 3.0
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  ;resc@cnLevelSelectionMode = "Explicitlevels"
  ;resc@cnLevels  = (/-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0,\
  ;                        0.3,0.6,0.9,1.2,1.5,1.8,2.1/)

;; vector res.
  resw = True
  resw@gsnFrame     = False
  resw@gsnDraw      = False
  resw@gsnMaximize  = False
  resw@gsnAddCyclic = False
  resw@vcGlyphStyle            = "LineArrow"
  resw@vcLineArrowThicknessF   = 2
  resw@vcRefMagnitudeF         = 4
  resw@vcRefLengthF            = 0.1
  resw@vcLineArrowColor        = "green4"
  resw@vcRefAnnoOn             = False
  resw@vcRefAnnoString2On      = False
  resw@vcRefAnnoString1On      = False
  resw@vcRefAnnoPerimSpaceF    = 0.25
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoJust           = "BottomRight"
  resw@vcRefAnnoBackgroundColor = -1
  resw@vcRefAnnoString1        = "$VMG$" +" ~F21~m/s"

  respoly = True
  respoly@gsLineThicknessF = 4.0
  respoly@gsLineColor  = color(0)

  respoly1 = respoly
  respoly1@gsLineDashPattern = 0
  respoly1@gsLineColor  = color(1)

  resP = True
  resP@gsnMaximize      = True
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelRowSpec  = True
  resP@gsnPanelTop      = 0.73
  resP@gsnPanelBottom   = 0.02
  resP@gsnPanelYWhiteSpacePercent = 3
  resP@gsnPanelXWhiteSpacePercent = 5

  resP@lbBoxLinesOn       = True
  resP@lbLabelFontHeightF = 0.012

  resP@pmLabelBarOrthogonalPosF = 0.
  resP@pmLabelBarWidthF         = 0.5
  resP@pmLabelBarHeightF        = 0.05

  resP@lbTitleOn = True
  resP@lbTitleString = "[K day~S~-1~N~]"
  resP@lbTitlePosition  = "Right"
  resP@lbTitleDirection = "Across"
  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.013
  resP@lbTitleStringParallelPosF = -0.1

  print("complete res. "+ systemfunc("date"))
;;----------------------------------------------------------
;; plot.
  resw@vcRefAnnoString1        = "$VMG$" +" ~F21~m s~S~-1~N~"
  resw@vcRefAnnoFontHeightF    = 0.02
  resw@vcRefAnnoOrthogonalPosF = 0.23;-1.5
  resw@vcRefAnnoParallelPosF   = 0.4

  plotv = new(4,graphic)
  plott = new(4,graphic)
  plotV = new(4,graphic)

  plotvt = new(4,graphic)
  plottt = new(4,graphic)
  plotVt = new(4,graphic)

  text2  = "~F21~-(v~F34~6~F21~T~F21~/~F34~6~F21~y~F21~)'"
  text4  = "~F21~(~F33~w~F33~s~F21~)'"

  l1 = (/"~F21~(c) ","~F21~(d) ","~F21~(e) ","~F21~(f) "/)
  l2 = (/"~F21~(g) ","~F21~(h) ","~F21~(i) ","~F21~(j) "/)
;u
  do n = 0,3
    resw@vcRefAnnoOn        = False
    resw@vcRefAnnoString1On = False

    rest@gsnLeftString  = l1(n)+text2+"@"+label(n)
    rest@gsnRightString = "P4"
    plotv(n) = gsn_csm_contour_map(wks, llv(n,:,:), rest)
    plott(n) = gsn_csm_contour(wks, llt(n,:,:), resc)
    plotV(n) = gsn_csm_vector(wks, llua(n,::2,::2),llva(n,::2,::2), resw)
    dum1  = gsn_add_polyline(wks, plotv(n), (/nlonw,nlone,nlone,nlonw,nlonw/), \
                                (/nlats,nlats,nlatn,nlatn,nlats/), respoly)
    overlay(plotv(n), dum1)
    overlay(plotv(n),plott(n))
    overlay(plotv(n),plotV(n))
;;Phase 8
    if (n.eq.3) then
      resw@vcRefAnnoOn         = True
      resw@vcRefAnnoString1On  = True
      resw@vcRefAnnoFontHeightF = 0.03
      resw@vcRefAnnoOrthogonalPosF = 0.16
      resw@vcRefAnnoParallelPosF   = 1.0
    end if
    rest@gsnLeftString  = l2(n)+text2+"@"+label(n)
    rest@gsnRightString = "P8"
    plotvt(n) = gsn_csm_contour_map(wks, llvt(n,:,:), rest)
    plottt(n) = gsn_csm_contour(wks, lltt(n,:,:), resc)
    plotVt(n) = gsn_csm_vector(wks, lluat(n,::2,::2),llvat(n,::2,::2), resw)
    dum1t  = gsn_add_polyline(wks, plotvt(n), (/nlonw,nlone,nlone,nlonw,nlonw/), \
                                (/nlats,nlats,nlatn,nlatn,nlats/), respoly)
    overlay(plotvt(n), dum1t)
    overlay(plotvt(n),plottt(n))
    overlay(plotvt(n),plotVt(n))
  end do

  gsn_panel(wks, (/plotv(0),plotvt(0),\
                   plotv(1),plotvt(1),\
                   plotv(2),plotvt(2),\
                   plotv(3),plotvt(3)/),(/2,2,2,2/), resP)

  i = 1
  do while (i.eq.1)
    i = 1
  end do