;--------------------------------------------------------
; draw temperature anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour.
; time: 2023-11-22
;; for MEJS: 50-110
;; for EAJS: 110-200
;--------------------------------------------------------

;-----------------------------------------------------------

  ddir = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  tdir   = ddir+"/NAJS_temp_budget/anom_"+"dtdt"+"_r2.5.nc"
  u1dir   = ddir+"/NAJS_temp_budget/"+"ERA5_daily_low_uvert_1979-2022_r2.5.nc"
  t1dir   = ddir+"/NAJS_temp_budget/"+"ERA5_daily_iso_tvert_1979-2022_r2.5.nc"

  utdir   = ddir+"/NAJS_temp_budget/anom_"+"udtdx"+"_r2.5.nc"
  vtdir   = ddir+"/NAJS_temp_budget/anom_"+"vdtdy"+"_r2.5.nc"
  wtdir   = ddir+"/NAJS_temp_budget/anom_"+"wdtdp"+"_r2.5.nc"

  udir   = ddir+"/"+"ERA5_daily_anom_u250_1979-2022_2.5x2.5.nc"
  vdir   = ddir+"/"+"ERA5_daily_anom_v250_1979-2022_2.5x2.5.nc"

  label = (/"(a) phase 8 lag -12","(b) phase 8 lag -10",\
           "(c) phase 8 lag -8","(d) phase 8 lag -6","(e) phase 8 lag -4"/)

  lats = 10
  latn = 50
  lonw = -110
  lone = -20

  tday = 24.*60.*60.
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  ft = addfile(tdir,"r")
  t  = ft->dtdt(:,{lats:latn},{lonw:lone})

  fut = addfile(utdir,"r")
  fvt = addfile(vtdir,"r")
  fwt = addfile(wtdir,"r")

  u  = fut->udtdx(:,{lats:latn},{lonw:lone})*tday*(-1)
  v  = fvt->vdtdy(:,{lats:latn},{lonw:lone})*tday*(-1)
  w  = fwt->wdtdp(:,{lats:latn},{lonw:lone})*tday*(-1)
  q  = t;-(u+w+v)
  copy_VarCoords(t,q)
  copy_VarCoords(t,u)
  copy_VarCoords(t,v)
  copy_VarCoords(t,w)

  time = t&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  delete([/t,ft/])

  ft = addfile(t1dir,"r")
  t1  = ft->t(:,{lats:latn},{lonw:lone})

  fu = addfile(udir,"r")
  u1 = fu->u;(:,{lats:latn},{lonw:lone})

  fv = addfile(vdir,"r")
  v1 = fv->v;(:,{lats:latn},{lonw:lone})

  print("complete data selling "+ systemfunc("date"))

;-------------------------------------------------------
;; selling MJO.

  delim = ","
  rmmf  = asciiread(rmmdir,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
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

;-------------------------------------------------------------

  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)

  amp   := amp(n1:n2)
  phase := phase(n1:n2)
  mon   := mon(n1:n2)
  time  := time(n1:n2)

  nt     = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  amp   := amp(nt)
  phase := phase(nt)

  delete(nt)
  ind_DJF = amp&timeind   ; indice of DJF between 1979-2022 DJF.

;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
  llt = new((/6,dimsizes(u&lat),dimsizes(u&lon)/),float)
  llt = 0.
  llv1 = new((/6,dimsizes(v1&lat),dimsizes(v1&lon)/),float)
  llv1 = 0.

  llu1 = new((/6,dimsizes(u1&lat),dimsizes(u1&lon)/),float)
  llu1 = 0.
  ;llu1 = llt
  ;llv1 = llt
  llu  = llt
  llv  = llt
  llw  = llt
  llq  = llt

  ip = 8
  nt = ind(phase.eq.ip.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 5
    ntt = 5 + j*5
    nindd = nind + ntt
    llt(j,:,:) = dim_avg_n(t1(nindd,:,:), 0)
    llu1(j,:,:) = dim_avg_n(u1(nindd,:,:), 0)
    llv1(j,:,:) = dim_avg_n(v1(nindd,:,:), 0)
    llu(j,:,:) = dim_avg_n(u(nindd,:,:), 0)
    llv(j,:,:) = dim_avg_n(v(nindd,:,:), 0)
    llw(j,:,:) = dim_avg_n(w(nindd,:,:), 0)
    llq(j,:,:) = dim_avg_n(q(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do

  llw = smth9(llw, 0.5, 0.25, False)
  llq = smth9(llq, 0.5, 0.25, False)

  copy_VarCoords(u(0,:,:), llt(0,:,:))
  copy_VarCoords(u(0,:,:), llu(0,:,:))
  copy_VarCoords(u(0,:,:), llv(0,:,:))
  copy_VarCoords(u(0,:,:), llw(0,:,:))
  copy_VarCoords(u(0,:,:), llq(0,:,:))
  copy_VarCoords(u1(0,:,:), llu1(0,:,:))
  copy_VarCoords(v1(0,:,:), llv1(0,:,:))

  ;llt  := dim_avg_n_Wrap(llt, 0)
  ;llu1 := dim_avg_n_Wrap(llu1, 0)
  ;llv1 := dim_avg_n_Wrap(llv1, 0)
  ;llu  := dim_avg_n_Wrap(llu, 0)
  ;llv  := dim_avg_n_Wrap(llv, 0)
  ;llw  := dim_avg_n_Wrap(llw, 0)
  ;llq  := dim_avg_n_Wrap(llq, 0)

  printVarSummary(llu1)

  print("complete calc leadlag data. "+ systemfunc("date"))

;; resources

;shading for T.
  cmap = read_colormap_file("MPL_RdBu")

  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = False
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 20
  rest@gsnLeftStringFontHeightF = 0.01
  rest@gsnRightStringFontHeightF = 0.01

  rest@tmXBLabelFont = 21
  rest@tmYLLabelFont = 21

  rest@vpWidthF  = 0.23
  rest@vpHeightF = 0.15

  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF   = lats
  rest@mpMaxLatF   = latn
  rest@mpMinLonF   = lonw
  rest@mpMaxLonF   = lone
  rest@mpLandFillColor = "white"

  rest@tmXBLabelFontHeightF = 0.008
  rest@tmYLLabelFontHeightF = 0.008

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -0.08
  rest@cnMaxLevelValF         =  0.08
  rest@cnLevelSpacingF        =  0.01
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = False

; contour phi plot for dtdt.
  resphi = True
  resphi@gsnDraw      = False
  resphi@gsnMaximize  = False
  resphi@gsnFrame     = False
  resphi@gsnAddCyclic = False
  resphi@gsnContourNegLineDashPattern = 1
  resphi@cnMonoLineColor    = True
  resphi@cnLinesOn          = True
  resphi@cnLineLabelsOn     = False
  resphi@cnLineColor        = "purple"
  resphi@cnLineThicknessF   = 2.0
  resphi@cnInfoLabelPerimOn = False
  resphi@cnInfoLabelOn      = False
  resphi@cnLevelSelectionMode = "Explicitlevels"
  ;resphi@cnLevels      = (/-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0.3,0.6,0.9,1.2,1.5,1.8/)

;; vector res.
  resw = True
  resw@gsnFrame     = False
  resw@gsnDraw      = False
  resw@gsnMaximize  = False
  resw@gsnAddCyclic = False
  resw@vcGlyphStyle            = "LineArrow"
  resw@vcLineArrowThicknessF   = 1.5
  resw@vcRefMagnitudeF         = 2
  resw@vcRefLengthF            = 0.05
  resw@vcLineArrowColor        = "black"
  resw@vcRefAnnoOn             = False
  resw@vcRefAnnoString2On      = False
  resw@vcRefAnnoString1On      = False
  resw@vcRefAnnoPerimSpaceF    = 0.4
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoFontHeightF    = 0.015
  resw@vcRefAnnoJust           = "BottomRight"
  resw@vcRefAnnoBackgroundColor = -1
  resw@vcRefAnnoString1        = "$VMG$" +" ~F21~m/s"

  resP = True
  resP@gsnMaximize      = False
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelRowSpec  = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02
  resP@gsnPanelYWhiteSpacePercent = 3

  resP@pmLabelBarOrthogonalPosF = -0.02
  resP@pmLabelBarWidthF         = 0.7
  resP@pmLabelBarHeightF        = 0.05

  resP@lbTitleFontHeightF= 0.012
  resP@lbBoxLinesOn      = True
  resP@lbLabelFontHeightF = 0.012

  print("complete res. "+ systemfunc("date"))

  respoly = True
  respoly@gsLineThicknessF = 3.0
  respoly@gsLineColor      = "green4"
;; plot.

  wks = gsn_open_wks("x11", pdir+"/"+"Fig8.budget_shading_MEJS")

  text2  = "~F21~-(u~F34~6~F21~T~F21~/~F34~6~F21~x~F21~)'"
  text3  = "~F21~-(v~F34~6~F21~T~F21~/~F34~6~F21~y~F21~)'"
  text4  = "~F21~(~F33~w~F33~s~F21~)'"
  text5  = "~F21~Q'/c~B~p~N~"

;u
  resw@vcRefAnnoOn             = False
  resw@vcRefAnnoString1On      = True
  resw@vcRefAnnoString1        = "$VMG$" +" ~F21~m/s"
  resw@vcRefAnnoFontHeightF    = 0.03
  resw@vcRefAnnoOrthogonalPosF = 0.44
  resw@vcRefAnnoParallelPosF   = 0.44 

  rest@gsnLeftString  = ""
  rest@gsnRightString = text2
  plotu   = gsn_csm_contour_map(wks, llu(0,:,:), rest)
  plott_t = gsn_csm_contour(wks, llt(0,:,:), resphi)
  plott_u = gsn_csm_vector(wks, llu1(0,::3,::3),llv1(0,::3,::3), resw)
  dum1  = gsn_add_polyline(wks, plotu, (/-70,-40,-40,-70,-70/), \
                          (/20,20,28,28,20/), respoly)
  overlay(plotu, dum1)
  overlay(plotu,plott_t)
  overlay(plotu,plott_u)
  respoly@gsLineColor  = "black"
; v
  rest@gsnLeftString  = ""
  rest@gsnRightString = text3
  plotv   = gsn_csm_contour_map(wks, llv(0,:,:), rest)
  plott_t = gsn_csm_contour(wks, llt(0,:,:), resphi)
  plott_v = gsn_csm_vector(wks, llu1(0,::2,::2),llv1(0,::2,::2), resw)
  dum1  = gsn_add_polyline(wks, plotv, (/-70,-40,-40,-70,-70/), \
                        (/20,20,28,28,20/), respoly)

  overlay(plotv,plott_t)
  overlay(plotv, dum1)
  overlay(plotv,plott_v)
;; w
  rest@gsnLeftString  = ""
  rest@gsnRightString = text4
  plotw   = gsn_csm_contour_map(wks, llq(0,:,:), rest)
  plott_t = gsn_csm_contour(wks, llt(0,:,:), resphi)
  plott_w = gsn_csm_vector(wks, llu1(0,::2,::2),llv1(0,::2,::2), resw)
  dum1  = gsn_add_polyline(wks, plotw, (/-70,-40,-40,-70,-70/), \
                        (/20,20,28,28,20/), respoly)

  overlay(plotw,plott_t)
  overlay(plotw, dum1)
  overlay(plotw,plott_w)
;; q
  rest@gsnLeftString  = ""
  rest@gsnRightString = text5
  plotq   = gsn_csm_contour_map(wks, llq(0,:,:), rest)
  plott_t = gsn_csm_contour(wks, llt(0,:,:), resphi)
  plott_q = gsn_csm_contour(wks, llt(0,:,:), resphi)
  dum1 = gsn_add_polyline(wks, plotq, (/-70,-40,-40,-70,-70/), \
                        (/20,20,28,28,20/), respoly)
  overlay(plotq,plott_t)
  overlay(plotq, dum1)

  gsn_panel(wks, (/plotv,plotw/), (/1,1/), resP)

  i = 1
  do while (i.eq.1)
    i = 1
  end do
