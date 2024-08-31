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
  utdir  = ddir+"/temp_budget/anom_"+"vdtdy"+"_r2.5.nc"
  adir = ddir+"/temp_budget/anom_vdtdy_seperate_r2.5.nc"

  lonwm = 10
  lonem = 110

  lonwe = 110
  lonee = 200

  stdate = -5

  label = (/"~F21~(a) phase 1 day -5","~F21~(b) phase 1 day 0",\
            "~F21~(c) phase 1 day 5","~F21~(d) phase 1 day 10",\
            "~F21~(e) phase 1 day 15","~F21~(e) phase 1 day 20",\
            "~F21~(e) phase 1 day 25"/)

  lats = 10
  latn = 50
  lonw = 10
  lone = 200

  tday = 24.*60.*60.

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  ft = addfile(tdir,"r")
  t  = ft->dtdt(:,{lats:latn},{lonw:lone})

  fut = addfile(utdir,"r")
  u  = fut->vdtdy(:,{lats:latn},{lonw:lone})*tday*(-1)

  f  = addfile(adir,"r")
  ultl = f->vltl(:,{lats:latn},{lonw:lone})*tday*(-1)
  ulti = f->vlti(:,{lats:latn},{lonw:lone})*tday*(-1)
  ulth = f->vlth(:,{lats:latn},{lonw:lone})*tday*(-1)
  uitl = f->vitl(:,{lats:latn},{lonw:lone})*tday*(-1)
  uiti = f->viti(:,{lats:latn},{lonw:lone})*tday*(-1)
  uith = f->vith(:,{lats:latn},{lonw:lone})*tday*(-1)
  uhtl = f->vhtl(:,{lats:latn},{lonw:lone})*tday*(-1)
  uhti = f->vhti(:,{lats:latn},{lonw:lone})*tday*(-1)
  uhth = f->vhth(:,{lats:latn},{lonw:lone})*tday*(-1)

  copy_VarCoords(t,u)
  copy_VarCoords(t,ultl)
  copy_VarCoords(t,ulti)
  copy_VarCoords(t,ulth)
  copy_VarCoords(t,uitl)
  copy_VarCoords(t,uiti)
  copy_VarCoords(t,uith)
  copy_VarCoords(t,uhtl)
  copy_VarCoords(t,uhti)
  copy_VarCoords(t,uhth)

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
  llt = new((/7,dimsizes(t&lat),dimsizes(t&lon)/),float)
  llt = 0.

  llv = new((/7,dimsizes(t&lat),dimsizes(t&lon)/),float)
  llv = 0.

  ip = 1
  nt = ind(phase.eq.ip.and.amp.gt.1.)
  nind = ind_DJF(nt)

  v = uitl

  do j = 0, 6
    ntt   = stdate + j*5
    nindd = nind + ntt
    llt(j,:,:) = dim_avg_n(u(nindd,:,:), 0)
    llv(j,:,:) = dim_avg_n(v(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do
  copy_VarCoords(t(0,:,:), llt(0,:,:))
  copy_VarCoords(t(0,:,:), llv(0,:,:))
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
  rest@gsnLeftStringFontHeightF = 0.006
  rest@gsnRightStringFontHeightF = 0.006

  rest@vpWidthF  = 0.23
  rest@vpHeightF = 0.15

  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF   = lats
  rest@mpMaxLatF   = latn
  rest@mpMinLonF   = lonw
  rest@mpMaxLonF   = lone
  rest@mpLandFillColor = "white"

  rest@tmXBLabelFontHeightF = 0.004
  rest@tmYLLabelFontHeightF = 0.004

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
  resphi@cnLineColor        = "black"
  resphi@cnLineThicknessF   = 2.5
  resphi@cnInfoLabelPerimOn = False
  resphi@cnInfoLabelOn      = False
  resphi@cnLevelSelectionMode = "Explicitlevels"
  resphi@cnLevels      = (/-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,\
                          0.3,0.6,0.9,1.2,1.5,1.8/)

  resP = True
  resP@gsnMaximize      = False
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02
  resP@gsnPanelYWhiteSpacePercent = 3

  resP@pmLabelBarOrthogonalPosF = 0.0
  resP@pmLabelBarWidthF         = 0.6

  resP@lbTitleString     = ""
  resP@lbTitlePosition   = "Right"
  resP@lbTitleDirection  = "Across"
  resP@lbTitleFontHeightF= 0.015
  resP@lbBoxLinesOn      = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012
  print("complete res. "+ systemfunc("date"))

  respoly = True
  respoly@gsLineThicknessF = 2.5
  respoly@gsLineColor      = "green4"


;; plot.
  wks = gsn_open_wks("x11", pdir+"/"+"dtdt_contour"+"_shading"+ip)

  plott   = new(7,graphic)
  plotphi = plott

  do n = 0,6
    rest@gsnLeftString = label(n)
    plott(n)   = gsn_csm_contour_map(wks, llv(n,:,:), rest)
    plotphi(n) = gsn_csm_contour(wks, llt(n,:,:), resphi)
    dum1  = gsn_add_polyline(wks, plott(n), (/60,90,90,60,60/), \
                            (/20,20,27,27,20/), respoly)
    dum2  = gsn_add_polyline(wks, plott(n), (/160,190,190,160,160/), \
                            (/25,25,33,33,25/), respoly)
    overlay(plott(n), dum1)
    overlay(plott(n), dum2)
    overlay(plott(n),plotphi(n))
  end do

  gsn_panel(wks, (/plott(0),plott(3),plott(1),plott(4),\
                plott(2),plott(5),plott(6)/), (/2,2,2,1/), resP)

  i = 1
  do while (i.eq.1)
    i = 1
  end do
