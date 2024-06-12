;--------------------------------------------------------
; draw temperature anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour.
; time: 2023-11-22
;; for MEJS: 50-110
;; for EAJS: 110-200
;--------------------------------------------------------
  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  tdir   = ddir+"/temp_budget/anom_"+"dtdt"+"_r2.5.nc"

  u1dir  = ddir+"/"+"ERA5_daily_low_uvert_1979-2022_r2.5.nc"
  t1dir  = ddir+"/"+"ERA5_daily_iso_tvert_1979-2022_r2.5.nc"

  utdir  = ddir+"/temp_budget/anom_"+"udtdx"+"_r2.5.nc"
  vtdir  = ddir+"/temp_budget/anom_"+"vdtdy"+"_r2.5.nc"
  wtdir  = ddir+"/temp_budget/anom_"+"wdtdp"+"_r2.5.nc"

  label = (/"day -15","day -10","day -5","day 0"/)

  lats = 10
  latn = 50
  lonw = 50
  lone = 220

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
  q  = (w+v);t-(u+w+v)
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
  t1 = ft->t(:,{lats:latn},{lonw:lone})

  fu = addfile(u1dir,"r")
  u1 = fu->u(:,{lats:latn},{lonw:lone})

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
  llt  = new((/4,dimsizes(t1&lat),dimsizes(t1&lon)/),float)
  llt  = 0.
  llu  = llt
  llw  = llt
  llq  = llt
  llv  = llt

  llua = llt
  llva = llt(0,:,:)

  nt = ind(phase.eq.2.or.phase.eq.3.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 3
    ntt   = -15 + j*5
    nindd = nind + ntt
    llt(j,:,:) = dim_avg_n(t1(nindd,:,:), 0)
    llua(j,:,:)= dim_avg_n(u1(nindd,:,:), 0)
    llu(j,:,:) = dim_avg_n(u(nindd,:,:), 0)
    llv(j,:,:) = dim_avg_n(v(nindd,:,:), 0)
    llw(j,:,:) = dim_avg_n(w(nindd,:,:), 0)
    llq(j,:,:) = dim_avg_n(q(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do

  llw = smth9(llw, 0.5, 0.25, False)
  llq = smth9(llq, 0.5, 0.25, False)

; (-7 - -3 8:12  )

  ;llt  = dim_avg_n(llt1(8:12,:,:), 0)
  ;llu1 = dim_avg_n(llu11(8:12,:,:), 0)
  ;llu  = dim_avg_n(llua(8:12,:,:), 0)
  ;llv  = dim_avg_n(llva(8:12,:,:), 0)
  ;llw  = dim_avg_n(llw1(8:12,:,:), 0)
  ;llq  = dim_avg_n(llq1(8:12,:,:), 0)
;; day -10 5 -5 10 0 15 5 20
  ;sdate = 20
  ;llt  = llt1(sdate,:,:)
  ;llu1 = llu11(sdate,:,:)
  ;llu  = llua(sdate,:,:)
  ;llv  = llva(sdate,:,:)
  ;llw  = llw1(sdate,:,:)
  ;llq  = llq1(sdate,:,:)

  copy_VarCoords(t1(0,:,:), llt(0,:,:))
  copy_VarCoords(t1(0,:,:), llua(0,:,:))
  copy_VarCoords(t1(0,:,:), llva)

  copy_VarCoords(t1(0,:,:), llw(0,:,:))
  copy_VarCoords(t1(0,:,:), llq(0,:,:))
  copy_VarCoords(t1(0,:,:), llu(0,:,:))
  copy_VarCoords(t1(0,:,:), llv(0,:,:))

  print("complete calc leadlag data. "+ systemfunc("date"))

;shading for T.
  cmap = read_colormap_file("MPL_RdBu")

  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = False
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 20
  rest@gsnLeftStringFontHeightF  = 0.03
  rest@gsnRightStringFontHeightF = 0.03

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
  resw@vcRefMagnitudeF         = 80
  resw@vcRefLengthF            = 0.1
  resw@vcLineArrowColor        = "black"
  resw@vcRefAnnoOn             = False
  resw@vcRefAnnoString2On      = False
  resw@vcRefAnnoString1On      = False
  resw@vcRefAnnoPerimSpaceF    = 0.4
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoFontHeightF    = 0.025
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

  resP@lbBoxLinesOn       = True
  resP@lbLabelFontHeightF = 0.012

  resP@pmLabelBarOrthogonalPosF = -0.04
  resP@pmLabelBarWidthF         = 0.7
  resP@pmLabelBarHeightF        = 0.05

  respoly = True
  respoly@gsLineThicknessF = 4.0
  respoly@gsLineColor      = "green3"

  respoly1 = respoly
  respoly1@gsLineDashPattern = 1

  print("complete res. "+ systemfunc("date"))
;;----------------------------------------------------------
;; plot.
  resw@vcRefAnnoString1        = "$VMG$" +" ~F21~m/s"
  resw@vcRefAnnoFontHeightF    = 0.02
  resw@vcRefAnnoOrthogonalPosF = 0.23;-1.5
  resw@vcRefAnnoParallelPosF   = 1.0

  wks = gsn_open_wks("x11", pdir+"/"+"Fig5_budget_shading_MEJS")

  plotu   = new(4,graphic)
  plott_u = new(4,graphic)
  plott_t = new(4,graphic)
  plotv   = new(4,graphic)
  plott_v = new(4,graphic)
  plotw   = new(4,graphic)
  plott_w = new(4,graphic)
  plotq   = new(4,graphic)
  plott_q = new(4,graphic)

  text2  = "~F21~-(u~F34~6~F21~T~F21~/~F34~6~F21~x~F21~)'"
  text3  = "~F21~-(v~F34~6~F21~T~F21~/~F34~6~F21~y~F21~)'"
  text4  = "~F21~(~F33~w~F33~s~F21~)'"
  text5  = "~F21~Q'/c~B~p~N~"

  l1 = (/"~F21~(a) ","~F21~(b) ","~F21~(c) ","~F21~(d) ","~F21~(e) "/)
  l2 = (/"~F21~(e) ","~F21~(f) ","~F21~(g) ","~F21~(h) ","~F21~(j) "/)
  l3 = (/"~F21~(i) ","~F21~(j) ","~F21~(k) ","~F21~(l) ","~F21~(o) "/)
  l4 = (/"~F21~(m) ","~F21~(n) ","~F21~(o) ","~F21~(p) ","~F21~(t) "/)
;u
  do n = 0,3
    resw@vcRefAnnoOn             = False
    resw@vcRefAnnoString1On      = False
    if (n.eq.3) then
      resw@vcRefAnnoOn         = True
      resw@vcRefAnnoString1On  = True
    end if
    rest@gsnLeftString  = l1(n)+label(n)
    rest@gsnRightString = text2
    plotu(n)   = gsn_csm_contour_map(wks, llu(n,:,:), rest)
    plott_t(n) = gsn_csm_contour(wks, llt(n,:,:), resc)
    plott_u(n) = gsn_csm_vector(wks, llua(n,::3,::3),llva(::3,::3), resw)
    dum1  = gsn_add_polyline(wks, plotu(n), (/60,100,100,60,60/), \
                                (/23,23,31,31,23/), respoly)
    ;dum2  = gsn_add_polyline(wks, plotu(n), (/105,140,140,105,105/), \
    ;                      (/25,25,35,35,25/), respoly1)
    overlay(plotu(n), dum1)
    ;overlay(plotu(n), dum2)
    overlay(plotu(n),plott_t(n))
    overlay(plotu(n),plott_u(n))
  end do
; v
  do n = 0,3
    rest@gsnLeftString  = l2(n)+label(n)
    rest@gsnRightString = text3
    plotv(n)   = gsn_csm_contour_map(wks, llv(n,:,:), rest)
    plott_v(n) = gsn_csm_contour(wks, llt(n,:,:), resc)
    dum1  = gsn_add_polyline(wks, plotv(n), (/60,100,100,60,60/), \
                                (/23,23,31,31,23/), respoly)
    ;dum2  = gsn_add_polyline(wks, plotv(n), (/105,140,140,105,105/), \
    ;                      (/25,25,35,35,25/), respoly1)
    overlay(plotv(n), dum1)
    ;overlay(plotv(n), dum2)
    overlay(plotv(n),plott_v(n))
  end do
;; w
  do n = 0,3
    rest@gsnLeftString  = l3(n)+label(n)
    rest@gsnRightString = text4
    plotw(n)   = gsn_csm_contour_map(wks, llw(n,:,:), rest)
    plott_w(n) = gsn_csm_contour(wks, llt(n,:,:), resc)
    dum1  = gsn_add_polyline(wks, plotw(n), (/60,100,100,60,60/), \
                                (/23,23,31,31,23/), respoly)
    ;dum2  = gsn_add_polyline(wks, plotw(n), (/105,140,140,105,105/), \
    ;                      (/25,25,35,35,25/), respoly1)
    overlay(plotw(n), dum1)
    ;overlay(plotw(n), dum2)
    overlay(plotw(n),plott_w(n))
  end do
;; q
  do n = 0,3
    rest@gsnLeftString  = l4(n)+label(n)
    rest@gsnRightString = text5
    plotq(n)   = gsn_csm_contour_map(wks, llq(n,:,:), rest)
    plott_q(n) = gsn_csm_contour(wks, llt(n,:,:), resc)
    dum1  = gsn_add_polyline(wks, plotq(n), (/60,100,100,60,60/), \
                                (/23,23,31,31,23/), respoly)
    ;dum2  = gsn_add_polyline(wks, plotq(n), (/105,140,140,105,105/), \
    ;                      (/25,25,35,35,25/), respoly1)
    overlay(plotq(n), dum1)
    ;overlay(plotq(n), dum2)
    overlay(plotq(n),plott_q(n))
  end do

  gsn_panel(wks, (/plotu(0),plotw(0),\
            plotu(1),plotw(1),\
            plotu(2),plotw(2),\
            plotu(3),plotw(3)/),(/2,2,2,2/), resP)

  i = 1
  do while (i.eq.1)
    i = 1
  end do