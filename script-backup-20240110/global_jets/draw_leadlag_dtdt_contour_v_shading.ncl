;--------------------------------------------------------
; draw temperature anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour.
; time: 2023-11-22
;--------------------------------------------------------

;-----------------------------------------------------------

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  tdir   = ddir+"/temp_budget/anom_"+"dtdt"+"_r2.5.nc"
  vdir   = ddir+"/temp_budget/anom_"+var+"_r2.5.nc"
  udir   = ddir+"/"+"ERA5_daily_u250_1979-2022.nc"

  utdir   = ddir+"/temp_budget/anom_"+"udtdx"+"_r2.5.nc"
  vtdir   = ddir+"/temp_budget/anom_"+"vdtdy"+"_r2.5.nc"
  wtdir   = ddir+"/temp_budget/anom_"+"wdtdp"+"_r2.5.nc"
;--------------------------------------------------------
; u250
;--------------------------------------------------------
  f    = addfile(udir,"r")

  time = f->time
  TIME = cd_calendar(time, 0)
  mon  = toint(TIME(:,1))
  nt   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)

  u    = dim_avg_n_Wrap(f->u(nt,0,:,:),0)

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
  resc@cnLineColor        = "magenta"
  resc@cnLineThicknessF   = 2.0
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels      = (/30.,40.,50./)

  delete([/time,f,TIME,mon,nt/])

  print("draw clim.u complete. "+ systemfunc("date"))
;--------------------------------------------------------

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  ft = addfile(tdir,"r")
  fv = addfile(vdir,"r")
  t  = ft->dtdt
  v  = (fv->$var$)*24.*60.*60.*(-1.)

  fut = addfile(utdir,"r")
  fvt = addfile(vtdir,"r")
  fwt = addfile(wtdir,"r")

  ut  = fut->udtdx
  vt  = fvt->vdtdy
  wt  = fwt->wdtdp
  ;v   = t + (ut+vt+wt)*24.*60.*60.
  copy_VarCoords(t,v)

  print("complete data selling "+ systemfunc("date"))

  time = t&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
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
  llt = new((/8,dimsizes(t&lat),dimsizes(t&lon)/),float)
  llt = 0.

  llv = new((/8,dimsizes(v&lat),dimsizes(v&lon)/),float)
  llv = 0.

  ip = 1
  nt = ind(phase.eq.ip.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 7
    ntt   = -15 + j*5
    nindd = nind + ntt
    llt(j,:,:) = dim_avg_n(t(nindd,:,:), 0)
    llv(j,:,:)   = dim_avg_n(v(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do
  copy_VarCoords(t(0,:,:), llt(0,:,:))
  copy_VarCoords(v(0,:,:), llv(0,:,:))
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
  rest@gsnMajorLonSpacing = 30
  ;rest@gsnCenterStringFontHeightF = 0.025 

  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF   = -20
  rest@mpMaxLatF   = 60
  rest@mpMinLonF   = 0
  rest@mpCenterLonF= 120
  rest@mpMaxLonF   = 240
  rest@mpLandFillColor = "white"

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
  resphi@cnLineColor        = "black"
  resphi@cnLineThicknessF   = 2.5
  resphi@cnInfoLabelPerimOn = False
  resphi@cnInfoLabelOn      = False
  resphi@cnLevelSelectionMode = "Explicitlevels"
  resphi@cnLevels      = (/-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4/)

; contour olr plot.
  resolr = True
  resolr@gsnDraw      = False
  resolr@gsnMaximize  = False
  resolr@gsnFrame     = False
  resolr@gsnAddCyclic = True
  resolr@gsnContourNegLineDashPattern = 1
  resolr@cnMonoLineColor    = True
  resolr@cnLinesOn          = True
  resolr@cnLineLabelsOn     = False
  resolr@cnLineColor        = "green4"
  resolr@cnLineThicknessF   = 4.0
  resolr@cnInfoLabelPerimOn = False
  resolr@cnInfoLabelOn      = False
  resolr@cnLevelSelectionMode = "Explicitlevels"
  resolr@cnLevels      = (/-20,-10,10,20/)

  resP = True
  resP@gsnMaximize      = True
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02

  resP@pmLabelBarOrthogonalPosF = 0.0
  resP@pmLabelBarWidthF         = 0.85

  resP@lbTitleString     = "Temp anom (K)"
  resP@lbTitlePosition   = "Right"
  resP@lbTitleDirection  = "Across"
  resP@lbTitleFontHeightF= 0.015
  resP@lbBoxLinesOn      = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012
  print("complete res. "+ systemfunc("date"))

  respoly = True
  respoly@gsLineThicknessF = 2.5
  respoly@gsLineColor      = "black"
;; plot.
  wks = gsn_open_wks("x11", pdir+"/"+"dtdt_contour"+var+"_shading"+ip)

  plott   = new(8,graphic)
  plotphi = plott

  do n = 0,7
    plott(n)   = gsn_csm_contour_map(wks, llv(n,:,:), rest)
    plotphi(n) = gsn_csm_contour(wks, llt(n,:,:), resphi)
    plotu = gsn_csm_contour(wks, u, resc)
    dum1  = gsn_add_polyline(wks, plott(n), (/60,80,80,60,60/), \
                                (/20,20,30,30,20/), respoly)
    dum2  = gsn_add_polyline(wks, plott(n), (/170,190,190,170,170/), \
                                (/25,25,35,35,25/), respoly)
    overlay(plott(n), dum1)
    overlay(plott(n), dum2)
    overlay(plott(n),plotu)
    overlay(plott(n),plotphi(n))
  end do

  resP@gsnPanelFigureStringsPerimOn = True
  resP@gsnPanelFigureStringsJust    = "Bottom Right"
  resP@gsnPanelFigureStringsFontHeightF = 0.012
  resP@gsnPanelFigureStrings = (/"(a) phase 1 lag -15","(e) phase 1 lag 5",\
                              "(b) phase 1 lag -10","(f) phase 1 lag 10",\
                              "(c) phase 1 lag -5","(g) phase 1 lag 15",\
                              "(d) phase 1 lag 0","(h) phase 1 lag 20"/)
  gsn_panel(wks, (/plott(0),plott(4),plott(1),plott(5),\
                plott(2),plott(6),plott(3),plott(7)/), (/2,2,2,2/), resP)

  i = 1
  do while (i.eq.1)
    i = 1
  end do
