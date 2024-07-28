;--------------------------------------------------------
; draw temperature anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour.
; time: 2023-11-22
;; for MEJS: 50-110
;; for EAJS: 110-200
;--------------------------------------------------------

;-----------------------------------------------------------

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  tdir   = ddir+"/temp_budget/anom_"+"dtdt"+"_r2.5.nc"

  u1dir   = ddir+"/"+"ERA5_daily_low_uvert_1979-2022_r2.5.nc"
  t1dir   = ddir+"/"+"ERA5_daily_iso_tvert_1979-2022_r2.5.nc"

  utdir   = ddir+"/temp_budget/anom_"+"udtdx"+"_r2.5.nc"
  vtdir   = ddir+"/temp_budget/anom_"+"vdtdy"+"_r2.5.nc"
  wtdir   = ddir+"/temp_budget/anom_"+"wdtdp"+"_r2.5.nc"

  lonwm  = 60
  lonem  = 150

  stdate = 2

  label2 = (/"phase 1 day 2","phase 1 day 4",\
            "phase 1 day 6","phase 1 day 8",\
            "phase 1 day 10"/)

  label = label2

  lats = 10
  latn = 50
  lonw = lonwm
  lone = lonem

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
  q  = t-(u+w+v)
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

  fu = addfile(u1dir,"r")
  u1  = fu->u(:,{lats:latn},{lonw:lone})

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
  llt1  = new((/6,dimsizes(t1&lat),dimsizes(t1&lon)/),float)
  llt1  = 0.
  llu11 = llt1
  llua  = llt1
  llva  = llt1
  llw1  = llt1
  llq1  = llt1
  llv1  = llt1(0,:,:)

  ip = 2
  nt = ind(phase.eq.ip.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 6
    ntt   = -6 + j
    nindd = nind + ntt
    llt1(j,:,:)  = dim_avg_n(t1(nindd,:,:), 0)
    llu11(j,:,:) = dim_avg_n(u1(nindd,:,:), 0)
    llua(j,:,:)  = dim_avg_n(u(nindd,:,:), 0)
    llva(j,:,:)  = dim_avg_n(v(nindd,:,:), 0)
    llw1(j,:,:)  = dim_avg_n(w(nindd,:,:), 0)
    llq1(j,:,:)  = dim_avg_n(q(nindd,:,:), 0)
    delete([/ntt,nindd/])
    end do

  llw1 = smth9(llw1, 0.5, 0.25, False)
  llq1 = smth9(llq1, 0.5, 0.25, False)

  llt  = dim_avg_n(llt1, 0)
  llu1 = dim_avg_n(llu11,0)
  llu  = dim_avg_n(llua, 0)
  llv  = dim_avg_n(llva, 0)
  llw  = dim_avg_n(llw1, 0)
  llq  = dim_avg_n(llq1, 0)

  copy_VarCoords(t1(0,:,:), llt)
  copy_VarCoords(t1(0,:,:), llu)
  copy_VarCoords(t1(0,:,:), llv)
  copy_VarCoords(t1(0,:,:), llw)
  copy_VarCoords(t1(0,:,:), llq)
  copy_VarCoords(t1(0,:,:), llu1)
  ;copy_VarCoords(t1(0,:,:), llv1)

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
  rest@gsnLeftStringFontHeightF = 0.025
  rest@gsnRightStringFontHeightF = 0.025

  rest@tmXBLabelFont = 21
  rest@tmYLLabelFont = 21

  ;rest@vpWidthF  = 0.23
  ;rest@vpHeightF = 0.15

  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF   = lats
  rest@mpMaxLatF   = latn
  rest@mpMinLonF   = lonw
  rest@mpMaxLonF   = lone
  rest@mpLandFillColor = "white"

  rest@tmXBLabelFontHeightF = 0.02
  rest@tmYLLabelFontHeightF = 0.02

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -2.1
  rest@cnMaxLevelValF         =  2.1
  rest@cnLevelSpacingF        =  0.3
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
  resphi@cnLineColor        = "mediumorchid4"
  resphi@cnLineThicknessF   = 3.0
  resphi@cnInfoLabelPerimOn = False
  resphi@cnInfoLabelOn      = False
  resphi@cnLevelSelectionMode = "Explicitlevels"
  resphi@cnLevels      = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0,1.2/)

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
  respoly@gsLineThicknessF = 4.0
  respoly@gsLineColor      = "green3"

  respoly1 = respoly
  respoly1@gsLineDashPattern = 1

;; plot.

  wks = gsn_open_wks("x11", pdir+"/"+"Fig5_budget_shading_MEJS")

  plotu   = new(1,graphic)
  plott_u = new(1,graphic)
  plott_t = new(1,graphic)
  plotv   = new(1,graphic)
  plott_v = new(1,graphic)
  plotw   = new(1,graphic)
  plott_w = new(1,graphic)
  plotq   = new(1,graphic)
  plott_q = new(1,graphic)

  text2  = "~F21~-(u~F34~6~F21~T~F21~/~F34~6~F21~x~F21~)'"
  text3  = "~F21~-(v~F34~6~F21~T~F21~/~F34~6~F21~y~F21~)'"
  text4  = "~F21~(~F33~w~F33~s~F21~)'"
  text5  = "~F21~Q'/c~B~p~N~"

  l1 = (/"~F21~(a) ","~F21~(b) ","~F21~(c) ","~F21~(d) ","~F21~(e) "/)
  l2 = (/"~F21~(f) ","~F21~(g) ","~F21~(h) ","~F21~(i) ","~F21~(j) "/)
  l3 = (/"~F21~(k) ","~F21~(l) ","~F21~(m) ","~F21~(n) ","~F21~(o) "/)
  l4 = (/"~F21~(p) ","~F21~(q) ","~F21~(r) ","~F21~(s) ","~F21~(t) "/)
;u
  resw@vcRefAnnoOn             = True
  resw@vcRefAnnoString1On      = True
  resw@vcRefAnnoString1        = "$VMG$" +" ~F21~m/s"
  resw@vcRefAnnoFontHeightF    = 0.02
  resw@vcRefAnnoOrthogonalPosF = -1.27
  resw@vcRefAnnoParallelPosF   = 1.0

  rest@gsnLeftString  = l1(0)+text2
  plotu   = gsn_csm_contour_map(wks, llu, rest)
  plott_t = gsn_csm_contour(wks, llt, resphi)
  plott_u = gsn_csm_vector(wks, llu1(::3,::3),llv1(::3,::3), resw)
  dum1  = gsn_add_polyline(wks, plotu, (/70,100,100,70,70/), \
                            (/25,25,30,30,25/), respoly)
  overlay(plotu,plott_t)

; v
  rest@gsnLeftString  = l1(1)+text3
  plotv   = gsn_csm_contour_map(wks, llv, rest)
  plott_v = gsn_csm_contour(wks, llt, resphi)
  dum2  = gsn_add_polyline(wks, plotv, (/70,100,100,70,70/), \
                            (/25,25,30,30,25/), respoly)
  overlay(plotv,plott_v)

;; w
  rest@gsnLeftString  = l1(2)+text4
  plotw   = gsn_csm_contour_map(wks, llw, rest)
  plott_w = gsn_csm_contour(wks, llt, resphi)
  dum3  = gsn_add_polyline(wks, plotw, (/70,100,100,70,70/), \
                            (/25,25,30,30,25/), respoly)
  overlay(plotw,plott_w)

;; q
  rest@gsnLeftString  = l1(3)+text5
  plotq   = gsn_csm_contour_map(wks, llq, rest)
  plott_q = gsn_csm_contour(wks, llt, resphi)
  dum4  = gsn_add_polyline(wks, plotq, (/70,100,100,70,70/), \
                            (/25,25,30,30,25/), respoly)
  overlay(plotq,plott_q)

  gsn_panel(wks, (/plotu,plotv,plotw,plotq/), (/2,2/), resP)

  i = 1
  do while (i.eq.1)
    i = 1
  end do