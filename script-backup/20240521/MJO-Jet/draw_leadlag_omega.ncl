;--------------------------------------------------------
; draw uwind anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour. with WAF.
; time: 2023-12-18
;--------------------------------------------------------
  latmin = -20
  latmax =  60
  lonmin =   0
  lonmax =  220

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  wpdir   = ddir+"/temp_budget/anom_"+"wdtdp"+"_r2.5.nc"
  wdir    = ddir+"/"+"ERA5_daily_anom_w500_1979-2022_r2.5.nc"
  olrdir  = ddir+"/"+"olr_daily_anom_1979-2022.nc"
  wsigdir = ddir+"/"+"monte_carlo/omega_sig_phase1.nc"
  ucdir   = ddir+"/"+"ERA5_daily_u250_1979-2022.nc"
;--------------------------------------------------------
; u250
;--------------------------------------------------------
  f    = addfile(ucdir,"r")

  time = f->time
  TIME = cd_calendar(time, 0)
  mon  = toint(TIME(:,1))
  nt   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  uc   = dim_avg_n_Wrap(f->u(nt,0,:,:),0)

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
  resc@cnLineThicknessF   = 3.0
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels      = (/30.,40./)

  delete([/time,f,TIME,mon,nt/])

  print("draw clim.u complete. "+ systemfunc("date"))

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  fw   = addfile(wdir,"r")
  fwp  = addfile(wpdir,"r")
  folr = addfile(olrdir,"r")

  olr = folr->olr(time|:,lat|:,lon|:)
  w   = fw->w(:,{latmin:latmax},{lonmin:lonmax})
  wp  = fwp->wdtdp(:,{latmin:latmax},{lonmin:lonmax})

  lat = w&lat
  lon = w&lon

  time = w&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  print("complete data selling "+ systemfunc("date"))
;; w sig.
  sigl = 0.05
  sigh = 0.95
  fsigw   = addfile(wsigdir,"r")
  m_value_w = fsigw->m_valueu

  m_valuew  = new((/2,8,dimsizes(w&lat),dimsizes(w&lon)/),float) 
  m_valuew  = 0.
  nl = toint(dimsizes(m_value_w(:,0,0,0))*sigl) - 1
  nh = toint(dimsizes(m_value_w(:,0,0,0))*sigh) - 1
  m_valuew(0,:,:,:) = m_value_w(nl,:,:,:)
  m_valuew(1,:,:,:) = m_value_w(nh,:,:,:)
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
  llw = new((/8,dimsizes(w&lat),dimsizes(w&lon)/),float)
  llw = 0.

  llwp = new((/8,dimsizes(wp&lat),dimsizes(wp&lon)/),float)
  llwp = 0.

  llolr = new((/8,dimsizes(olr&lat),dimsizes(olr&lon)/),float)
  llolr = 0.

  nt = ind(phase.eq.1.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 7
    ntt   = -15 + j*5
    nindd = nind + ntt
    llw(j,:,:) = dim_avg_n(w(nindd,:,:), 0)
    llwp(j,:,:)  = dim_avg_n(wp(nindd,:,:), 0)
    llolr(j,:,:) = dim_avg_n(olr(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do
  copy_VarCoords(w(0,:,:), llw(0,:,:))
  copy_VarCoords(wp(0,:,:), llwp(0,:,:))
  copy_VarCoords(olr(0,:,:), llolr(0,:,:))
  print("complete calc leadlag data. "+ systemfunc("date"))

;shading T plot
  cmap = read_colormap_file("MPL_RdBu")

  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = False
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 30
  rest@gsnRightString     = "~F21~500hPa_~F33~w"
  rest@gsnLeftStringFontHeightF = 0.022 
  rest@gsnRightStringFontHeightF = 0.02

  rest@tmXTOn = False
  rest@tmYROn = False
  rest@tmYLLabelFontHeightF = 0.015
  rest@tmXBLabelFontHeightF = 0.015

  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF   = -20
  rest@mpMaxLatF   = 60
  rest@mpMinLonF   = 0
  rest@mpCenterLonF= 120
  rest@mpMaxLonF   = 220
  rest@mpLandFillColor = "white"

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -0.03
  rest@cnMaxLevelValF         =  0.03
  rest@cnLevelSpacingF        =  0.0025
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = False

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
  resolr@cnLineThicknessF   = 3.0
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
  resP@gsnPanelYWhiteSpacePercent = 5

  resP@pmLabelBarOrthogonalPosF = -0.018
  resP@pmLabelBarWidthF         = 0.85

  resP@lbTitleString     = "~F21~ unit: Pa/s"
  resP@lbTitleFontHeightF= 0.013
  resP@lbTitlePosition   = "Bottom"
  ;resP@lbTitleDirection  = "Across"
  ;resP@lbTitleFontHeightF= 0.015
  ;resP@lbBoxLinesOn      = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012
  print("complete res. "+ systemfunc("date"))

  respoly = True
  respoly@gsLineThicknessF = 2.5
  respoly@gsLineColor      = "green4"

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
  sres@cnFillDotSizeF     = 0.0025
  sres@gsnLeftString      =" "
  sres@gsnRightString     =" "
  sres@cnMonoFillColor    = True
  sres@cnFillColor        = "black"
  ;sres@cnLevels    = (/-0.1,0,0.1/) 

  opt = True
  opt@gsnShadeFillType   = "pattern"
  opt@gsnShadeHigh      = 17
  ;opt@gsnShadeMid        = 10
  opt@gsnShadeFillScaleF   = 1.
;; plot.

  wks = gsn_open_wks("x11", pdir+"/"+"u-shading_phi_contour_olr_contour_waf_vector_phase")

  plott   = new(8,graphic)
  plotct  = plott
  plotolr = plott
  dum     = plott

  label = (/"(a) phase 1 lag -15","(e) phase 1 lag -10",\
                              "(b) phase 1 lag -5","(f) phase 1 lag 0",\
                              "(c) phase 1 lag 5","(g) phase 1 lag 10",\
                              "(d) phase 1 lag 15","(h) phase 1 lag 20"/)
  do n = 0,7
    rest@gsnLeftString = label(n)
    xt = llw(n,:,:)
    xt = 0.
    do i = 0, dimsizes(lat)-1
      do j = 0, dimsizes(lon)-1
        if(llw(n,i,j).lt.m_valuew(0,n,i,j).or.llw(n,i,j).gt.m_valuew(1,n,i,j)) then
          xt(i,j) = 2.
        else
          xt(i,j) = -2.
        end if
      end do
    end do
    plott(n)   = gsn_csm_contour_map(wks, llw(n,:,:), rest)
    plotolr(n) = gsn_csm_contour(wks, llolr(n,:,:), resolr)
    plotu = gsn_csm_contour(wks, uc, resc)
    dum1  = gsn_add_polyline(wks, plott(n), (/60,80,80,60,60/), \
                                (/20,20,30,30,20/), respoly)
    dum2  = gsn_add_polyline(wks, plott(n), (/170,190,190,170,170/), \
                                (/25,25,35,35,25/), respoly)
    dum(n)   = gsn_csm_contour(wks,xt,sres)
    dum(n)   = gsn_contour_shade(dum(n),1.9,-999.,opt)
    overlay(plott(n),dum(n))
    overlay(plott(n), dum1)
    overlay(plott(n), dum2)
    overlay(plott(n),plotu)
    ;overlay(plott(n),plotphi(n))
    overlay(plott(n),plotolr(n))
  end do

  resP@gsnPanelFigureStringsPerimOn = False
  gsn_panel(wks, (/plott(0),plott(4),plott(1),plott(5),\
                plott(2),plott(6),plott(3),plott(7)/), (/2,2,2,2/), resP)

  i = 1
  do while (i.eq.1)
    i = 1
  end do
