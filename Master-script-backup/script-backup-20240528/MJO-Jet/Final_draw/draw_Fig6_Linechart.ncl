;-----------------------------------------------------------------
;; MEJS
;; draw T tendency and T anomaly in Line chart.
;; draw OLR, udtdx and omega.
;; time: 2024-05-02
;; two part: (a) Line chart.
;-----------------------------------------------------------------
  rad   = 4.*atan(1.)/180. 
  tday  = 24.*60.*60.
  mlatn = 30
  mlats = 20
  mlonw = 50
  mlone = 90

  elatn = 35
  elats = 25
  elonw = 160
  elone = 200

  lats = 10
  latn = 50
  lonw = 40
  lone = 210

  ddir = "/home/sunming/data5/cuixy/global_jets/data"
  tdir = ddir+"/temp_budget/anom_"+"dtdt"+"_r2.5.nc"
  udir = ddir+"/temp_budget/"+"anom_udtdx_r2.5.nc"
  wdir = ddir+"/temp_budget/"+"anom_wdtdp_r2.5.nc"

  u1dir  = ddir+"/"+"ERA5_daily_low_uvert_1979-2022_r2.5.nc"
  t1dir  = ddir+"/"+"ERA5_daily_iso_tvert_1979-2022_r2.5.nc"

  ttdir = ddir+"/ERA5_daily_anom_tvert_1979-2022_2.5x2.5.nc"
  pdir = "/home/sunming/data5/cuixy/global_jets/plot/Fig4_bar_MEJS"

  label = (/"day -15","day -10","day -5","day 0"/)
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
;; read in T ten and anom.
  ft  = addfile(tdir,"r")
  ftt = addfile(ttdir,"r")

  t = ft->dtdt
  ta  = ftt->t
;; read in udtdx and omega.
  fu  = addfile(udir,"r")
  fw  = addfile(wdir,"r")

  u = fu->udtdx(:,{lats:latn},{lonw:lone})*tday*(-1)
  w = fw->wdtdp(:,{lats:latn},{lonw:lone})*tday*(-1)

  ft = addfile(t1dir,"r")
  t1 = ft->t(:,{lats:latn},{lonw:lone})

  fu = addfile(u1dir,"r")
  u1 = fu->u(:,{lats:latn},{lonw:lone})

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
;;MEJS.
  nt = ind(phase.eq.2.and.amp.gt.1.)
  nind = ind_DJF(nt)

  llmejs = new((/41,2/),float)
  llmejs = 0.
  tm = t(:,{mlats:mlatn},{mlonw:mlone})
  clatjs = cos(tm&lat*rad)
  tajs = wgt_areaave_Wrap(ta(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  tjs  = wgt_areaave_Wrap(t(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  do j = 0, 40
    ntt   = -15 + j
    nindd = nind + ntt
    llmejs(j,0) = dim_avg_n(tjs(nindd), 0)
    llmejs(j,1) = dim_avg_n(tajs(nindd), 0)
    delete([/ntt,nindd/])
  end do
;; line for dtdt
  lme = new((/2,41/),float)
  lme(0,:) = llmejs(:,0)
  lme(1,:) = llmejs(:,1)
  delete([/tajs,tjs,clatjs,llmejs/])
;;EAJS.
  lleajs = new((/41,2/),float)
  lleajs = 0.
  tm = t(:,{elats:elatn},{elonw:elone})
  clatjs = cos(tm&lat*rad)
  tajs = wgt_areaave_Wrap(ta(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  tjs  = wgt_areaave_Wrap(t(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  do j = 0, 40
    ntt   = -15 + j
    nindd = nind + ntt
    lleajs(j,0) = dim_avg_n(tjs(nindd), 0)
    lleajs(j,1) = dim_avg_n(tajs(nindd), 0)
    delete([/ntt,nindd/])
  end do
;; line for dtdt
  lea = new((/2,41/),float)
  lea(0,:) = lleajs(:,0)
  lea(1,:) = lleajs(:,1)
  delete([/tajs,tjs,clatjs,lleajs/])

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
  res@tmXBLabels = (/"-15","-10","-5","0","5","10","15","20","25" /)
  res@tmXBLabelFontHeightF = 0.012
  res@tiXAxisFontHeightF = 0.012

  res@tiYAxisFontHeightF = 0.012

  res@tmYLMode        = "Manual"
  res@tmYLMinorPerMajor = 0
  res@tmYLTickSpacingF  = 0.3
  res@tmYLLabelFontHeightF = 0.012

  res@tmYLMajorLengthF   = -0.01
  res@trYMaxF  = 1.5
  res@trYMinF  = -1.5
  res@gsnYRefLine = 0.
  res@gsnYRefLineDashPattern = 0

  res@xyLineThicknessF = 4.0

  color = (/"blue3","red3"/)  ;"royalblue4"

  xlabel = ispan(0,40,1)
;----------------------------------------------------------------
;; Line plot
  wks = gsn_open_wks("x11", pdir)
  ;drawNDCGrid(wks)
  res@vpWidthF  = 0.5
  res@vpHeightF = 0.2
  res@vpXF = 0.25
  res@vpYF = 0.97
  res@trYMaxF  =  0.5
  res@trYMinF  = -0.5
  res@trXMaxF  =  25
  res@tmYLTickSpacingF  = 0.1
  res@gsnLeftString  = "~F21~(a) (~F34~6~F21~T/~F34~6~F21~t~F21~)' & T' "
  ;res@gsnRightString = "--- (~F34~6~F21~T/~F34~6~F21~t~F21~)' - T'"
  res@xyLineColor = color(0)
  res@xyDashPattern = 1

  resR = True
  resR@trYMaxF  =  4.0
  resR@trYMinF  = -4.0
  resR@xyLineThicknessF = 4.0
  resR@xyLineColor = color(0)
  resR@xyDashPattern = 0
  resR@tiYAxisString  = "~F21~Temp anomaly (K)"
  dum = gsn_csm_xy2(wks, xlabel, lme(0,:),lme(1,:), res,resR)

  res@xyLineColor = color(1)
  resR@xyLineColor = color(1) 
  dumt = gsn_csm_xy2(wks, xlabel, lea(0,:),lea(1,:), res,resR) 

  overlay(dum,dumt)

  resl = True
  resl@gsLineColor = "black"
  resl@gsLineThicknessF  = 4
  resl@gsLineDashPattern = 1

  restx = True
  restx@txFontHeightF = 0.015

  gsn_polyline_ndc(wks, (/0.52,0.58/), (/0.98,0.98/), resl)
  gsn_text_ndc(wks,"~F21~(~F34~6~F21~T/~F34~6~F21~t~F21~)'", 0.62,0.98,restx)

  resl@gsLineDashPattern = 0
  gsn_polyline_ndc(wks, (/0.68,0.73/), (/0.98,0.98/), resl)
  gsn_text_ndc(wks,"~F21~T'", 0.75,0.98,restx)

  draw(dum)
;; for 2-D plot.
  llt  = new((/4,dimsizes(t1&lat),dimsizes(t1&lon)/),float)
  llt  = 0.
  llu  = llt
  llw  = llt
  llua = llt
  llva = llt(0,:,:)

  do j = 0, 3
    ntt   = -10 + j*3
    nindd = nind + ntt
    llt(j,:,:) = dim_avg_n(t1(nindd,:,:), 0)
    llu(j,:,:) = dim_avg_n(u(nindd,:,:), 0)
    llw(j,:,:) = dim_avg_n(w(nindd,:,:), 0)
    llua(j,:,:)= dim_avg_n(u1(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do
  llw = smth9(llw, 0.5, 0.25, False)

  copy_VarCoords(t1(0,:,:), llt(0,:,:))
  copy_VarCoords(t1(0,:,:), llua(0,:,:))
  copy_VarCoords(t1(0,:,:), llva)

  copy_VarCoords(t1(0,:,:), llw(0,:,:))
  copy_VarCoords(t1(0,:,:), llu(0,:,:))
;shading for T.
  cmap = read_colormap_file("MPL_RdBu")

  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = False
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 20
  rest@gsnLeftStringFontHeightF  = 0.025
  rest@gsnRightStringFontHeightF = 0.025

  rest@tmXBLabelFont = 21
  rest@tmYLLabelFont = 21

  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF   = lats
  rest@mpMaxLatF   = latn
  rest@mpMinLonF   = lonw
  rest@mpMaxLonF   = lone
  rest@mpLandFillColor = "white"

  rest@tmXBLabelFontHeightF = 0.015
  rest@tmYLLabelFontHeightF = 0.015

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -1.5
  rest@cnMaxLevelValF         =  1.5
  rest@cnLevelSpacingF        =  0.3
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
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels  = (/-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0,\
                             0.3,0.6,0.9,1.2,1.5,1.8,2.1/)

;; vector res.
  resw = True
  resw@gsnFrame     = False
  resw@gsnDraw      = False
  resw@gsnMaximize  = False
  resw@gsnAddCyclic = False
  resw@vcGlyphStyle            = "LineArrow"
  resw@vcLineArrowThicknessF   = 1.5
  resw@vcRefMagnitudeF         = 40
  resw@vcRefLengthF            = 0.05
  resw@vcLineArrowColor        = "black"
  resw@vcRefAnnoOn             = False
  resw@vcRefAnnoString2On      = False
  resw@vcRefAnnoString1On      = False
  resw@vcRefAnnoPerimSpaceF    = 0.4
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoFontHeightF    = 0.01
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
  resP@gsnMaximize      = False
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelRowSpec  = True
  resP@gsnPanelTop      = 0.73
  resP@gsnPanelBottom   = 0.02
  resP@gsnPanelYWhiteSpacePercent = 3

  resP@lbBoxLinesOn       = True
  resP@lbLabelFontHeightF = 0.012

  resP@pmLabelBarOrthogonalPosF = 0.
  resP@pmLabelBarWidthF         = 0.5
  resP@pmLabelBarHeightF        = 0.05

  resP@lbTitleOn = True
  resP@lbTitleString = "[K/day]"
  resP@lbTitlePosition  = "Right"
  resP@lbTitleDirection = "Extend"
  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.013
  resP@lbTitleStringParallelPosF = -0.1

  print("complete res. "+ systemfunc("date"))
;;----------------------------------------------------------
;; plot.
  resw@vcRefAnnoString1        = "$VMG$" +" ~F21~m/s"
  resw@vcRefAnnoFontHeightF    = 0.02
  resw@vcRefAnnoOrthogonalPosF = 0.23;-1.5
  resw@vcRefAnnoParallelPosF   = 0.4

  plotu   = new(4,graphic)
  plott_u = new(4,graphic)
  plott_t = new(4,graphic)
  plotw   = new(4,graphic)
  plott_w = new(4,graphic)

  text2  = "~F21~-(u~F34~6~F21~T~F21~/~F34~6~F21~x~F21~)'"
  text4  = "~F21~(~F33~w~F33~s~F21~)'"

  label = (/"day-10","day-7","day-4","day-1"/)
  l1 = (/"~F21~(b) ","~F21~(c) ","~F21~(d) ","~F21~(e) "/)
  l2 = (/"~F21~(f) ","~F21~(g) ","~F21~(h) ","~F21~(i) "/)
;u
  do n = 0,3
    resw@vcRefAnnoOn             = False
    resw@vcRefAnnoString1On      = False
    if (n.eq.3) then
      resw@vcRefAnnoOn         = True
      resw@vcRefAnnoString1On  = True
    end if
    rest@gsnLeftString  = l1(n)+text2+"@"+label(n)
    ;rest@gsnRightString = text2
    plotu(n)   = gsn_csm_contour_map(wks, llu(n,:,:), rest)
    plott_t(n) = gsn_csm_contour(wks, llt(n,:,:), resc)
    plott_u(n) = gsn_csm_vector(wks, llua(n,::3,::3),llva(::3,::3), resw)
    dum1  = gsn_add_polyline(wks, plotu(n), (/mlonw,mlone,mlone,mlonw,mlonw/), \
                                (/20,20,30,30,20/), respoly)
    dum2  = gsn_add_polyline(wks, plotu(n), (/160,200,200,160,160/), \
                          (/25,25,35,35,25/), respoly1)
    overlay(plotu(n), dum1)
    overlay(plotu(n), dum2)
    overlay(plotu(n),plott_t(n))
    overlay(plotu(n),plott_u(n))
  end do

;; w
  do n = 0,3
    rest@gsnLeftString  = l2(n)+text4+"@"+label(n)
    ;rest@gsnRightString = text4
    plotw(n)   = gsn_csm_contour_map(wks, llw(n,:,:), rest)
    plott_w(n) = gsn_csm_contour(wks, llt(n,:,:), resc)
    dum1  = gsn_add_polyline(wks, plotw(n), (/mlonw,mlone,mlone,mlonw,mlonw/), \
                                (/20,20,30,30,20/), respoly)
    dum2  = gsn_add_polyline(wks, plotw(n), (/160,200,200,160,160/), \
                          (/25,25,35,35,25/), respoly1)
    overlay(plotw(n), dum1)
    overlay(plotw(n), dum2)
    overlay(plotw(n),plott_w(n))
  end do

  gsn_panel(wks, (/plotu(0),plotw(0),\
                   plotu(1),plotw(1),\
                   plotu(2),plotw(2),\
                   plotu(3),plotw(3)/),(/2,2,2,2/), resP)
  i = 1
  do while (i.eq.1)
    i = 1
  end do


