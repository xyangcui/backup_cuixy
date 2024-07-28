;-----------------------------------------------------------------
;; draw general budget terms. 
;; MEJS box (10-70E,10-40N)
;; EAJS box (110-180E,15-50N)
;; time: 2024-01-12
;; m-value(low:high,wme:vme:wea:vea,dates,level,lat)
;-----------------------------------------------------------------
;;---------------------------------------------------------------------
;; This function reads topographic data off a binary file as "short"
;; and converts it to float.  It also attaches lat/lon coordinate
;; variables.

  lats = -10
  latn =  50
  lonw =  70
  lone =  210

  ddir = "/home/sunming/data5/cuixy/global_jets/data/meri_circu"

  u2dir = "/home/sunming/data5/cuixy/global_jets/data/ERA5_daily_iso_uvert_1979-2022_r2.5.nc"
  v2dir = "/home/sunming/data5/cuixy/global_jets/data/ERA5_daily_iso_vvert_1979-2022_r2.5.nc"
  t2dir = "/home/sunming/data5/cuixy/global_jets/data/ERA5_daily_low_tvert_1979-2022_r2.5.nc"

  olrdir = "/home/sunming/data5/cuixy/global_jets/data/olr_daily_anom_1979-2022.nc"

  sigdir = "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/Fig13_uv_p5_sig.nc"
  sig2dir = "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/Fig13_olr_p5_sig.nc"

  pdir = "/home/sunming/data5/cuixy/global_jets/plot/Fig10_UV_OLR_EAJS"

  shptp = "/home/sunming/data5/cuixy/data/add_TP/DBATP_Line.shp"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
  fsig   = addfile(sigdir,"r")
  mvalue = fsig->m_value
  fsig2  = addfile(sig2dir,"r")
  molr   = fsig2->m_olr

  fu = addfile(u2dir,"r")
  u2  = fu->u(:,{lats:latn},{lonw:lone})

  fv = addfile(v2dir,"r")
  v2  = fv->v(:,{lats:latn},{lonw:lone})

;; vdtdy
  ft = addfile(t2dir,"r")
  t2 = ft->t(:,{lats:latn},{lonw:lone})

  folr = addfile(olrdir,"r")

  olr = folr->olr(:,{lats:latn},{lonw:lone})

;---------------------------------------------------------------=
;; MJO phases.
  time = folr->time
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
;----------------------------------------------------------------------
  llolra = new((/5,dimsizes(olr&lat),dimsizes(olr&lon)/),float)
  llolra = 0.

  llu2 = new((/5,dimsizes(u2&lat),dimsizes(u2&lon)/),float)
  llv2 = new((/5,dimsizes(u2&lat),dimsizes(u2&lon)/),float)
  llt2 = new((/5,dimsizes(t2&lat),dimsizes(t2&lon)/),float)

  llolrt = new((/5,dimsizes(olr&lat),dimsizes(olr&lon)/),float)
  llolrt = 0.

  llu3 = new((/5,dimsizes(u2&lat),dimsizes(u2&lon)/),float)
  llv3 = new((/5,dimsizes(u2&lat),dimsizes(u2&lon)/),float)
  llt3 = new((/5,dimsizes(t2&lat),dimsizes(t2&lon)/),float)

  nt = ind(phase.eq.5.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 4
    ntt   = -7 + j*1
    nindd = nind + ntt
    llolra(j,:,:) = dim_avg_n(olr(nindd,:,:), 0)
    llu2(j,:,:) = dim_avg_n(u2(nindd,:,:), 0)
    llv2(j,:,:) = dim_avg_n(v2(nindd,:,:), 0)
    llt2(j,:,:) = dim_avg_n(t2(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do

  llolr = dim_avg_n(llolra,0)
  llu = dim_avg_n(llu2,0)
  llv = dim_avg_n(llv2,0)
  llt = dim_avg_n(llt2,0)

  copy_VarCoords(olr(0,:,:), llolr)
  copy_VarCoords(u2(0,:,:), llu)
  copy_VarCoords(v2(0,:,:), llv)
  copy_VarCoords(t2(0,:,:), llt)

  do j = 0, 4
    ntt   = -12 + j*1
    nindd = nind + ntt
    llolrt(j,:,:) = dim_avg_n(olr(nindd,:,:), 0)
    llu3(j,:,:) = dim_avg_n(u2(nindd,:,:), 0)
    llv3(j,:,:) = dim_avg_n(v2(nindd,:,:), 0)
    llt3(j,:,:) = dim_avg_n(t2(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do

  llolr1 = dim_avg_n(llolrt,0)
  llu1 = dim_avg_n(llu3,0)
  llv1 = dim_avg_n(llv3,0)
  llt1 = dim_avg_n(llt3,0)

  copy_VarCoords(olr(0,:,:), llolr1)
  copy_VarCoords(u2(0,:,:), llu1)
  copy_VarCoords(v2(0,:,:), llv1)
  copy_VarCoords(t2(0,:,:), llt1)

  print("complete calc leadlag data. "+ systemfunc("date"))

;; res for meri circulation.
  cmap = read_colormap_file("MPL_RdBu")

;; vector res.
  resw = True
  resw@gsnFrame     = False
  resw@gsnDraw      = False
  resw@gsnMaximize  = False
  resw@gsnAddCyclic = False
  resw@vcGlyphStyle            = "CurlyVector"
  resw@vcLineArrowThicknessF   = 2.0
  resw@vcRefMagnitudeF         = 5
  resw@vcRefLengthF            = 0.05
  resw@vcLineArrowColor        = "gray"
  resw@vcRefAnnoOn             = False
  resw@vcRefAnnoString2On      = False
  resw@vcRefAnnoPerimSpaceF    = 0.4
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoFontHeightF    = 0.015
  resw@vcRefAnnoBackgroundColor = -1

  reswt = True
  reswt@gsnFrame     = False
  reswt@gsnDraw      = False
  reswt@gsnMaximize  = False
  reswt@gsnAddCyclic = False

  reswt@vcGlyphStyle            = "CurlyVector"
  reswt@vcLineArrowThicknessF   = 2.5
  reswt@vcRefMagnitudeF         = 5
  reswt@vcRefLengthF            = 0.05
  reswt@vcLineArrowColor        = "black"
  reswt@vcRefAnnoOn             = True
  reswt@vcRefAnnoString2On      = False
  reswt@vcRefAnnoString1On      = True
  reswt@vcRefAnnoPerimSpaceF    = 0.4
  reswt@vcRefAnnoPerimOn        = False
  reswt@vcRefAnnoFontHeightF    = 0.015
  reswt@vcRefAnnoBackgroundColor = -1
  reswt@vcRefAnnoString1        = "$VMG$" +" ~F21~m/s"
  reswt@vcRefAnnoOrthogonalPosF =  0.12
  reswt@vcRefAnnoParallelPosF   =  1.0

;; res for olr.
  cmap1 = read_colormap_file("BlueWhiteOrangeRed")

  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = False
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 30

  rest@gsnLeftStringFontHeightF  = 0.015 
  rest@gsnRightStringFontHeightF = 0.015
  rest@gsnLeftStringOrthogonalPosF = 0.01
  rest@gsnRightStringOrthogonalPosF = rest@gsnLeftStringOrthogonalPosF

  rest@tmXTOn = False
  rest@tmYROn = False
  rest@tmYLLabelFontHeightF = 0.015
  rest@tmXBLabelFontHeightF = 0.015

  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF   = lats
  rest@mpMaxLatF   = latn
  rest@mpMinLonF   = lonw
  rest@mpMaxLonF   = lone
  rest@mpLandFillColor = "white"

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap1;(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -15
  rest@cnMaxLevelValF         =  15
  rest@cnLevelSpacingF        =  3
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = False
  rest@pmLabelBarOrthogonalPosF =  0.1
  rest@pmLabelBarParallelPosF   =  0.45

; contour phi plot for dtdt.
  resphi = True
  resphi@gsnDraw      = False
  resphi@gsnMaximize  = False
  resphi@gsnFrame     = False
  resphi@gsnAddCyclic = False
  resphi@gsnContourNegLineDashPattern = 1
  resphi@cnMonoLineColor    = True
  resphi@cnLinesOn          = True
  resphi@cnLineLabelsOn     = True
  resphi@cnLineColor        = "mediumorchid4"
  resphi@cnLineThicknessF   = 3.0
  resphi@cnInfoLabelPerimOn = False
  resphi@cnInfoLabelOn      = False
  resphi@cnLineLabelPlacementMode = "computed"

  resP = True
  resP@gsnMaximize      = False
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelRowSpec  = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02
  resP@gsnPanelYWhiteSpacePercent = 3

  resP@pmLabelBarOrthogonalPosF = 0.02
  ;resP@pmLabelBarParallelPosF   = 0.8
  resP@pmLabelBarWidthF         = 0.6
  resP@pmLabelBarHeightF        = 0.1

  resP@lbBoxLinesOn      = True
  resP@lbLabelFontHeightF = 0.012

  sres   = True 
  sres@cnLineDrawOrder    = "PostDraw"                    
  sres@gsnDraw            = False                          
  sres@gsnFrame           = False                        
  sres@cnLineLabelsOn     = False                              
  sres@cnLinesOn          = False                           
  sres@cnInfoLabelOn      = False
  sres@cnFillOn           = False                       
  sres@cnInfoLabelOn      = False
  sres@cnFillScaleF       = 0.5
  sres@cnFillDotSizeF     = 0.003
  sres@gsnLeftString      =" "
  sres@gsnRightString     =" "
  sres@cnMonoFillColor    = True
  sres@cnFillColor        = "gray"
  sres@gsnPresHgtHeightLabelOn = False 

  opt = True
  opt@gsnShadeFillType   = "pattern"
  opt@gsnShadeHigh      = 17
  opt@gsnShadeFillScaleF   = 1.5

  wks = gsn_open_wks("x11",pdir)

  rest@gsnLeftString  = "~F21~ (b) Phase 5 day 0"
  rest@gsnRightString = "~F21~ OLR_UV(ISO)_T(LOW)"

  x = llolr
  xt = x
  xt = 0.
  do i = 0, dimsizes(olr&lat)-1
    do j = 0, dimsizes(olr&lon)-1
      if(x(i,j).lt.molr(0,i,j).or.x(i,j).gt.molr(3,i,j)) then
        xt(i,j) = 2.
      else
        xt(i,j) = -2.
      end if
    end do
  end do

  plot  = gsn_csm_contour_map(wks, llolr, rest)
  dumolr = gsn_csm_contour(wks,xt,sres)
  dumolr = gsn_contour_shade(dumolr,1.9,-999.,opt)
  overlay(plot,dumolr)

;;plot for uv.
  undef("x")
  x = llu
  y = llv

  x = mask(x, ( llu.lt.mvalue(1,0,:,:).or.llu.gt.mvalue(2,0,:,:) ).or.( llv.lt.mvalue(1,1,:,:).or.llv.gt.mvalue(2,1,:,:) ), True)

  y = mask(y, ( llu.lt.mvalue(1,0,:,:).or.llu.gt.mvalue(2,0,:,:) ).or.( llv.lt.mvalue(1,1,:,:).or.llv.gt.mvalue(2,1,:,:) ), True)

  plotu = gsn_csm_vector(wks, llu(::2,::2),llv(::2,::2), resw)
  plotv = gsn_csm_vector(wks, x(::2,::2), y(::2,::2), reswt)

  plott = gsn_csm_contour(wks, llt, resphi)

  overlay(plot,plott)

  respoly = True
  respoly@gsLineThicknessF = 3.5
  respoly@gsLineColor      = "green3"
  respoly@gsLineDashPattern  = 1
  dum  = gsn_add_polyline(wks, plot, (/160,190,190,160,160/), \
                      (/30,30,40,40,30/), respoly)

  restext = True
  restext@txFontColor = "green3"
  restext@txFontHeightF = 0.02
  dumtext = gsn_add_text(wks, plot, "~F22~A", 108, 28, restext)

  dumtext2 = gsn_add_text(wks, plot, "~F22~C", 152, 38, restext)

  overlay(plot, plotu)
  overlay(plot, plotv)

;; draw day -5

  sigdir = "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/Fig13_uv_p5_sig_day-5.nc"
  sig2dir = "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/Fig13_olr_p5_sig_day-5.nc"
  fsig   = addfile(sigdir,"r")
  mvalue = fsig->m_value
  fsig2  = addfile(sig2dir,"r")
  molr   = fsig2->m_olr
  
  delete(x)
  reswt@vcRefAnnoOn             = False
  reswt@vcRefAnnoOrthogonalPosF = 0.1

  rest@gsnLeftString  = "~F21~ (a) Phase 5 day -10"
  rest@gsnRightString = "~F21~ OLR_UV(ISO)_T(LOW)"

  x = llolr1
  xt = x
  xt = 0.
  do i = 0, dimsizes(olr&lat)-1
    do j = 0, dimsizes(olr&lon)-1
      if(x(i,j).lt.molr(0,i,j).or.x(i,j).gt.molr(3,i,j)) then
        xt(i,j) = 2.
      else
        xt(i,j) = -2.
      end if
    end do
  end do

  plot1  = gsn_csm_contour_map(wks, llolr1, rest)
  dumolr1 = gsn_csm_contour(wks,xt,sres)
  dumolr1 = gsn_contour_shade(dumolr1,1.9,-999.,opt)
  overlay(plot1,dumolr1)

;;plot for uv.
  undef("x")
  x = llu1
  y = llv1

  x = mask(x, ( llu1.lt.mvalue(1,0,:,:).or.llu1.gt.mvalue(2,0,:,:) ).or.( llv1.lt.mvalue(1,1,:,:).or.llv1.gt.mvalue(2,1,:,:) ), True)

  y = mask(y, ( llu1.lt.mvalue(1,0,:,:).or.llu1.gt.mvalue(2,0,:,:) ).or.( llv1.lt.mvalue(1,1,:,:).or.llv1.gt.mvalue(2,1,:,:) ), True)

  plotu1 = gsn_csm_vector(wks, llu1(::2,::2),llv1(::2,::2), resw)
  plotv1 = gsn_csm_vector(wks, x(::2,::2), y(::2,::2), reswt)

  plott = gsn_csm_contour(wks, llt, resphi)

  overlay(plot1,plott)

  respoly = True
  respoly@gsLineThicknessF = 3.5
  respoly@gsLineColor      = "green3"
  respoly@gsLineDashPattern  = 1
  dum1  = gsn_add_polyline(wks, plot1, (/160,190,190,160,160/), \
                      (/25,25,35,35,25/), respoly)

  restext = True
  restext@txFontColor = "green3"
  restext@txFontHeightF = 0.02
  dumtext = gsn_add_text(wks, plot1, "~F22~A", 80, 28, restext)

  overlay(plot1, plotu1)
  overlay(plot1, plotv1)

  gsn_panel(wks, (/plot1,plot/), (/1,1/), resP)

  i = 0
  do while (i.eq.0)
    i = 0
  end do
