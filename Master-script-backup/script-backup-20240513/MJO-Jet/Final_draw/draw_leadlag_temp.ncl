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
  udir   = ddir+"/"+"ERA5_daily_u250_1979-2022.nc"

  u1dir   = ddir+"/"+"ERA5_daily_low_uvert_1979-2022_r2.5.nc"
  t1dir   = ddir+"/"+"ERA5_daily_iso_tvert_1979-2022_r2.5.nc"

  utdir   = ddir+"/temp_budget/anom_"+"udtdx"+"_r2.5.nc"
  vtdir   = ddir+"/temp_budget/anom_"+"vdtdy"+"_r2.5.nc"
  wtdir   = ddir+"/temp_budget/anom_"+"wdtdp"+"_r2.5.nc"

  lonwm = 60
  lonem = 150

  lonwe = 110
  lonee = 200

  stdate = -7

  label1 = (/"phase 3 day -7","phase 3 day -5",\
            "phase 3 day -3","phase 3 day 20",\
            "phase 1 day 25"/)

  label2 = (/"phase 1 day -5","phase 1 day 0",\
            "phase 1 day 5","phase 1 day 10",\
            "phase 1 day 15"/)

  label = label1

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
  q  = v+w
  copy_VarCoords(t,q)
  copy_VarCoords(t,u)
  copy_VarCoords(t,v)
  copy_VarCoords(t,w)

  print("complete data selling "+ systemfunc("date"))

  time = t&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  delete([/ft/])

  ft = addfile(t1dir,"r")
  t1  = ft->t(:,{lats:latn},{lonw:lone})

  fu = addfile(u1dir,"r")
  u1  = fu->u(:,{lats:latn},{lonw:lone})

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
  llt = new((/3,dimsizes(t&lat),dimsizes(t&lon)/),float)
  llt = 0.

  llt1 = new((/3,dimsizes(t1&lat),dimsizes(t1&lon)/),float)
  llt1 = 0.

  llu1 = llt1
  llv1 = llt1(0,:,:)

  llu = llt
  llv = llt
  llw = llt
  llq = llt

  ip = 3
  nt = ind(phase.eq.ip.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 2
    ntt = stdate+j*2
    nindd = nind + ntt
    llt(j,:,:) = dim_avg_n(t(nindd,:,:), 0)
    llu(j,:,:) = dim_avg_n(u(nindd,:,:), 0)
    llv(j,:,:) = dim_avg_n(v(nindd,:,:), 0)
    llw(j,:,:) = dim_avg_n(w(nindd,:,:), 0)
    llq(j,:,:) = dim_avg_n(q(nindd,:,:), 0)
    llt1(j,:,:) = dim_avg_n(t1(nindd,:,:), 0)
    llu1(j,:,:) = dim_avg_n(u1(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do

  llw = smth9(llw, 0.5, 0.25, False)
  llq = smth9(llq, 0.5, 0.25, False)

  copy_VarCoords(t(0,:,:), llt1(0,:,:))
  copy_VarCoords(t(0,:,:), llu(0,:,:))
  copy_VarCoords(t(0,:,:), llv(0,:,:))
  copy_VarCoords(t(0,:,:), llw(0,:,:))
  copy_VarCoords(t(0,:,:), llq(0,:,:))
  copy_VarCoords(t1(0,:,:), llu1(0,:,:))
  copy_VarCoords(t1(0,:,:), llv1)
  print("complete calc leadlag data. "+ systemfunc("date"))

;; resources

;shading T plot
  cmap = read_colormap_file("MPL_RdBu")

  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = False
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 15
  rest@gsnLeftStringFontHeightF = 0.01
  rest@gsnRightStringFontHeightF = 0.01

  ;;font
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
  rest@cnMinLevelValF         = -1.8
  rest@cnMaxLevelValF         =  1.8
  rest@cnLevelSpacingF        =  0.3
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = False

; contour phi plot.
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
  resphi@cnLevels      = (/-1.0,-0.9,-0.8,-0.7,-0.6,-0.5\
                        ,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)

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

  resP@pmLabelBarOrthogonalPosF = 0.0
  resP@pmLabelBarWidthF         = 0.6

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

  wks = gsn_open_wks("x11", pdir+"/"+"dtdt_contour"+"_shading"+ip)

  plotu   = new(3,graphic)
  plott_u = new(3,graphic)
  plott_t = new(3,graphic)
  plotv   = new(3,graphic)
  plott_v = new(3,graphic)
  plotw   = new(3,graphic)
  plott_w = new(3,graphic)
  plotq   = new(3,graphic)
  plott_q = new(3,graphic)

  text2  = "~F21~-(u~F34~6~F21~T~F21~/~F34~6~F21~x~F21~)'"
  text3  = "~F21~-(v~F34~6~F21~T~F21~/~F34~6~F21~y~F21~)'"
  text4  = "~F21~(~F33~w~F33~s~F21~)'"
  text5  = "~F21~Q'/c~B~p~N~"

  l1 = (/"~F21~(a) ","~F21~(b) ","~F21~(c) ","~F21~(d) ","~F21~(e) "/)
  l2 = (/"~F21~(f) ","~F21~(g) ","~F21~(h) ","~F21~(i) ","~F21~(j) "/)
  l3 = (/"~F21~(k) ","~F21~(l) ","~F21~(m) ","~F21~(n) ","~F21~(o) "/)
  l4 = (/"~F21~(p) ","~F21~(q) ","~F21~(r) ","~F21~(s) ","~F21~(t) "/)
;u
  do n = 0,2
    rest@gsnLeftString  = l1(n)+label(n)
    rest@gsnRightString = text2
    plotu(n)   = gsn_csm_contour_map(wks, llu(n,:,:), rest)
    plott_t(n) = gsn_csm_contour(wks, llt1(n,:,:), resphi)
    plott_u(n) = gsn_csm_vector(wks, llu1(n,::3,::3),llv1(::3,::3), resw)
    dum1  = gsn_add_polyline(wks, plotu(n), (/60,90,90,60,60/), \
                            (/20,20,27,27,20/), respoly)
    dum2  = gsn_add_polyline(wks, plotu(n), (/160,190,190,160,160/), \
                            (/25,25,33,33,25/), respoly)
    overlay(plotu(n), dum1)
    overlay(plotu(n), dum2)
    overlay(plotu(n),plott_t(n))
    overlay(plotu(n),plott_u(n))
  end do
; v
  do n = 0,2
    rest@gsnLeftString  = l2(n)+label(n)
    rest@gsnRightString = text3
    plotv(n)   = gsn_csm_contour_map(wks, llv(n,:,:), rest)
    plott_v(n) = gsn_csm_contour(wks, llt1(n,:,:), resphi)
    dum1  = gsn_add_polyline(wks, plotv(n), (/60,90,90,60,60/), \
                            (/20,20,27,27,20/), respoly)
    dum2  = gsn_add_polyline(wks, plotv(n), (/160,190,190,160,160/), \
                            (/25,25,33,33,25/), respoly)
    overlay(plotv(n), dum1)
    overlay(plotv(n), dum2)
    overlay(plotv(n),plott_v(n))
  end do
;; w
  do n = 0,2
    rest@gsnLeftString  = l3(n)+label(n)
    rest@gsnRightString = text4
    plotw(n)   = gsn_csm_contour_map(wks, llw(n,:,:), rest)
    plott_w(n) = gsn_csm_contour(wks, llt(n,:,:), resphi)
    dum1  = gsn_add_polyline(wks, plotw(n), (/60,90,90,60,60/), \
                            (/20,20,27,27,20/), respoly)
    dum2  = gsn_add_polyline(wks, plotw(n), (/160,190,190,160,160/), \
                            (/25,25,33,33,25/), respoly)
    overlay(plotw(n), dum1)
    overlay(plotw(n), dum2)
    overlay(plotw(n),plott_w(n))
  end do
;; q
  do n = 0,2
    rest@gsnLeftString  = l4(n)+label(n)
    rest@gsnRightString = text5
    plotq(n)   = gsn_csm_contour_map(wks, llq(n,:,:), rest)
    plott_q(n) = gsn_csm_contour(wks, llt(n,:,:), resphi)
    dum1  = gsn_add_polyline(wks, plotq(n), (/60,90,90,60,60/), \
                            (/20,20,27,27,20/), respoly)
    dum2  = gsn_add_polyline(wks, plotq(n), (/160,190,190,160,160/), \
                            (/25,25,33,33,25/), respoly)
    overlay(plotq(n), dum1)
    overlay(plotq(n), dum2)
    overlay(plotq(n),plott_q(n))
  end do

  gsn_panel(wks, (/plotu(0),plotv(0),plotw(0),plotq(0),\
                plotu(1),plotv(1),plotw(1),plotq(1),\
                plotu(2),plotv(2),plotw(2),plotq(2)/),\
                (/4,4,4/), resP)

  i = 1
  do while (i.eq.1)
    i = 1
  end do
