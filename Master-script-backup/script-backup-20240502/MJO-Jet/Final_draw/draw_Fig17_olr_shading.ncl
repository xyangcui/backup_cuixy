;--------------------------------------------------------
; draw uwind anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour. with WAF.
; time: 2023-12-18
;--------------------------------------------------------
  lats =  0
  latn =  80
  lonw =  -150
  lone =  0

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  olrdir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
  wdir = ddir+"/NAJS_temp_budget/ERA5_daily_anom_w500_1979-2022_r2.5.nc"

  udir = ddir+"/"+"NAJS_temp_budget/ERA5_daily_anom_ulevel_1979-2022_r2.5.nc"
  vdir = ddir+"/"+"NAJS_temp_budget/ERA5_daily_anom_vlevel_1979-2022_r2.5.nc"
  tdir = ddir+"/NAJS_temp_budget/"+"ERA5_daily_low_tvert_1979-2022_r2.5.nc"

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  folr = addfile(olrdir,"r")

  olr = lonFlip(folr->olr)

  olr := olr(:,{-30:90},:)

  lat = olr&lat
  lon = olr&lon

  time = olr&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  fw = addfile(wdir,"r")

  w  = fw->w(:,{0:50},:)

;; vdtdy
  ;ft = addfile(tdir,"r")
  ;t = ft->t(:,{lats:latn},{lonw:lone})

  fu = addfile(udir,"r")
  u  = fu->u(:,{250},{lats:latn},{lonw:lone})

  fv = addfile(vdir,"r")
  v  = fv->v(:,{250},{lats:latn},{lonw:lone})

  print("complete data selling "+ systemfunc("date"))

  printVarSummary(u)
  printVarSummary(olr)
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
  llolr = new((/6,dimsizes(olr&lat),dimsizes(olr&lon)/),float)
  llw = new((/6,dimsizes(w&lat),dimsizes(w&lon)/),float)

  llu = new((/6,dimsizes(u&lat),dimsizes(u&lon)/),float)
  llv = new((/6,dimsizes(v&lat),dimsizes(v&lon)/),float)

  nt = ind(phase.eq.8.and.amp.gt.1.)
  nind = ind_DJF(nt)

  sdate = 5
  do j = 0, 5
    ntt   = sdate + j*5
    nindd = nind + ntt
    llolr(j,:,:) = dim_avg_n(olr(nindd,:,:), 0)
    llw(j,:,:) = dim_avg_n(w(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do
  copy_VarCoords(olr(0,:,:), llolr(0,:,:))
  copy_VarCoords(w(0,:,:), llw(0,:,:))

  llolr = smth9_Wrap(llolr, 0.5, 0.25, False)
  llw = smth9_Wrap(llw, 0.5, 0.25, False)

  do j = 0, 5
    ntt   = sdate + j*5
    nindd = nind + ntt
    llu(j,:,:) = dim_avg_n(u(nindd,:,:), 0)
    llv(j,:,:) = dim_avg_n(v(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do

  copy_VarCoords(u(0,:,:), llu(0,:,:))
  copy_VarCoords(v(0,:,:), llv(0,:,:))

  print("complete calc leadlag data. "+ systemfunc("date"))

;shading omega plot.
  cmap = read_colormap_file("BlueWhiteOrangeRed")

  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = False
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 30
  rest@gsnLeftStringFontHeightF  = 0.015 
  rest@gsnRightStringFontHeightF = 0.015
  rest@gsnLeftStringOrthogonalPosF  = 0.01
  rest@gsnRightStringOrthogonalPosF = 0.01

  rest@tmXTOn = False
  rest@tmYROn = False
  rest@tmYLLabelFontHeightF = 0.012
  rest@tmXBLabelFontHeightF = 0.012

  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF   = -25;10;
  rest@mpMaxLatF   = 90;65; 
  rest@mpMinLonF   = 80;-110; 
  rest@mpMaxLonF   = 360;-30; 
  rest@mpCenterLonF = 180.
  rest@mpLandFillColor = "white"

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap;(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -5;-0.024;;
  rest@cnMaxLevelValF         =  5;0.024;;
  rest@cnLevelSpacingF        =  1;0.004;;
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = True
  rest@pmLabelBarWidthF       = 0.5
  rest@pmLabelBarOrthogonalPosF = 0.3

;; vector res.
  resw = True
  resw@gsnFrame     = False
  resw@gsnDraw      = False
  resw@gsnMaximize  = False
  resw@gsnAddCyclic = False

  resw@vcGlyphStyle            = "CurlyVector"
  resw@vcLineArrowThicknessF   = 1.5
  resw@vcRefMagnitudeF         = 3.0
  resw@vcRefLengthF            = 0.05
  resw@vcLineArrowColor        = "black"
  resw@vcRefAnnoOn             = False
  resw@vcRefAnnoString2On      = False
  resw@vcRefAnnoString1On      = False
  resw@vcRefAnnoPerimSpaceF    = 0.4
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoFontHeightF    = 0.015
  resw@vcRefAnnoBackgroundColor = -1
  resw@vcRefAnnoString1        = "$VMG$" +" ~F21~m/s"

  respoly = True
  respoly@gsLineThicknessF = 3.5
  respoly@gsLineColor      = "green4"

;; plot.

  wks = gsn_open_wks("x11", pdir+"/"+"Fig7.NAJS_olr_omega500_UV")

  ;rest@vpXF = 0.25
  ;rest@vpYF = 0.95
  ;rest@vpWidthF  = 0.57
  ;rest@vpHeightF = 0.35
  rest@gsnLeftString = "~F21~(a) omega"
  rest@gsnRightString = "~F21~ ave: phase 8 day -10 to -5"

  plott = gsn_csm_contour_map(wks, llolr(0,:,:), rest)
  ;dum1  = gsn_add_polyline(wks, plott, (/-70,-40,-40,-70,-70/), \
  ;                        (/20,20,28,28,20/), respoly)
;; draw uv.
  plotuv  = gsn_csm_vector(wks, llu(0,::2,::2), llv(0,::2,::2), resw)

  respoly@gsLineDashPattern = 1
  dum2  = gsn_add_polyline(wks, plott, (/-100,-80,-80,-100,-100/), \
                          (/33,33,45,45,33/), respoly)

  respoly@gsLineDashPattern = 0
  dum3  = gsn_add_polyline(wks, plott, (/-70,-40,-40,-70,-70/), \
                          (/20,20,28,28,20/), respoly)

  overlay(plott, plotuv)

  draw(plott)

  i = 1
  do while (i.eq.1)
    i = 1
  end do
