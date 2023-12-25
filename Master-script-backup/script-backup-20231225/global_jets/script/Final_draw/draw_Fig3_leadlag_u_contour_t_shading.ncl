;--------------------------------------------------------
; draw temperature anomaly in contour, MTG
; anomaly in shading and OLR anomaly in contour.
; time: 2023-11-22
;--------------------------------------------------------

;-----------------------------------------------------------

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  udir   = ddir+"/"+"ERA5_daily_anom_u250_1979-2022_2.5x2.5.nc"
  mtgdir = ddir+"/"+"ERA5_daily_anom_mtg_1979-2022_2.5x2.5.nc"
  tdir   = ddir+"/"+"ERA5_daily_anom_tvert_1979-2022_2.5x2.5.nc"
  ucdir   = ddir+"/"+"ERA5_daily_u250_1979-2022.nc"
  tsigdir = ddir+"/"+"monte_carlo/tvert_sig_phase1&5.nc"
  usigdir = ddir+"/"+"monte_carlo/u_sig_phase15.nc"
  latmin = 0
  latmax = 70
;--------------------------------------------------------
; u250
;--------------------------------------------------------
  f    = addfile(ucdir,"r")

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
  resc@cnLineThicknessF   = 3.0
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels      = (/30.,40./)

  delete([/time,f,TIME,mon,nt/])

print("draw clim.u complete. "+ systemfunc("date"))
;--------------------------------------------------------

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  fu = addfile(udir,"r")
  ft   = addfile(tdir,"r")

  mtg  = fu->u(:,{-30:90},:)
  t    = ft->t(:,{latmin:latmax},:)

print("complete data selling "+ systemfunc("date"))

  time = mtg&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  lat  = mtg&lat
  lon  = mtg&lon
;--------------------------------------------------------------
;; t sig.
;; sig
  sigl = 0.05
  sigh = 0.95
  fsigu   = addfile(usigdir,"r")
  fsigt   = addfile(tsigdir,"r")
  m_value_u  = fsigu->m_valueu
  m_value_t  = fsigt->m_valuet

  m_valuet  = new((/2,2,8,dimsizes(t&lat),dimsizes(t&lon)/),float) 
  m_valuet  = 0.
  nl  = toint(dimsizes(m_value_t(:,0,0,0,0))*sigl) - 1
  nh  = toint(dimsizes(m_value_t(:,0,0,0,0))*sigh) - 1
  m_valuet(0,:,:,:,:) = m_value_t(nl,:,:,:,:)
  m_valuet(1,:,:,:,:) = m_value_t(nh,:,:,:,:)

  m_valueu  = new((/2,2,8,dimsizes(mtg&lat),dimsizes(mtg&lon)/),float) 
  m_valueu  = 0.
  nl = toint(dimsizes(m_value_u(:,0,0,0,0))*sigl) - 1
  nh = toint(dimsizes(m_value_u(:,0,0,0,0))*sigh) - 1
  m_valueu(0,:,:,:,:) = m_value_u(nl,:,:,:,:)
  m_valueu(1,:,:,:,:) = m_value_u(nh,:,:,:,:)

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

;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
  llmtg = new((/8,dimsizes(mtg&lat),dimsizes(mtg&lon)/),float)
  copy_VarCoords(mtg(0,:,:), llmtg(0,:,:))
  llmtg = 0.

  llt = new((/8,dimsizes(t&lat),dimsizes(t&lon)/),float)
  copy_VarCoords(t(0,:,:), llt(0,:,:))
  llt = 0.

  i = 0
    ip = i + 1
    nt = ind(phase.eq.ip.and.amp.gt.1.)
    nind = ind_DJF(nt)
    do j = 0, 7
        ntt   = -15 + j*5
        nindd = nind + ntt
        llmtg(j,:,:) = dim_avg_n(mtg(nindd,:,:), 0)
        llt(j,:,:)   = dim_avg_n(t(nindd,:,:), 0)
        delete([/ntt,nindd/])
    end do
    delete([/ip,nt,nind/])

  ;llmtg = llmtg * 10^(6)
  print("complete calc leadlag data. "+ systemfunc("date"))

;; resources

;shading T plot
  cmap = read_colormap_file("MPL_RdBu")

  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = True
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 30
  rest@gsnLeftStringFontHeightF = 0.02

  rest@tmXBLabelFont = 21
  rest@tmYLLabelFont = 21
  rest@tmXBLabelFontHeightF = 0.015
  rest@tmYLLabelFontHeightF = 0.015

  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF   = 0
  rest@mpMaxLatF   = 70
  rest@mpMinLonF   = -10
  rest@mpMaxLonF   = 220
  rest@mpLandFillColor = "white"

  rest@vpHeightF   = 0.23

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

; contour phi plot.
  resu = True
  resu@gsnDraw      = False
  resu@gsnMaximize  = False
  resu@gsnFrame     = False
  resu@gsnAddCyclic = True
  resu@gsnContourNegLineDashPattern = 1
  resu@cnMonoLineColor    = True
  resu@cnLinesOn          = True
  resu@cnLineLabelsOn     = False
  resu@cnLineColor        = "black"
  resu@cnLineThicknessF   = 2.5
  resu@cnInfoLabelPerimOn = False
  resu@cnInfoLabelOn      = False
  resu@cnLevelSelectionMode = "Explicitlevels"
  resu@cnLevels      = (/-8.,-6.,-4.,-2.,2.,4.,6.,8./)
  ;resu@cnLevelSelectionMode   = "Manuallevels"
  ;resu@cnMinLevelValF         = -8.0
  ;resu@cnMaxLevelValF         =  8.0
  ;resu@cnLevelSpacingF        =  2.0

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
  resP@gsnMaximize      = False
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02

  resP@pmLabelBarOrthogonalPosF = -0.03
  resP@pmLabelBarWidthF         = 0.7

  resP@lbTitleString     = "~F21~Temp anom (K)"
  resP@lbTitlePosition   = "Bottom"
  resP@lbTitleDirection  = "Across"
  resP@lbTitleFontHeightF= 0.015
  resP@lbBoxLinesOn      = True
  resP@gsnPanelRowSpec    = True
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
  sres@cnFillScaleF       = 0.1
  sres@gsnLeftString      =" "
  sres@gsnRightString     =" "
  sres@cnMonoFillColor    = True
  sres@cnFillColor        = "black"
  ;sres@cnLevels    = (/-0.1,0,0.1/) 

  opt = True
  opt@gsnShadeFillType   = "pattern"
  opt@gsnShadeHigh      = 10
  ;opt@gsnShadeMid        = 10
  opt@gsnShadeFillScaleF   = 1.

  print("complete res. "+ systemfunc("date"))

;; plot.
  Label1 = (/"~F21~(a) phase 1 day -15","~F21~(b) phase 1 day -10",\
                              "~F21~(c) phase 1 day -5","~F21~(d) phase 1 day 0",\
                              "~F21~(e) phase 1 day 5","~F21~(f) phase 1 day 10",\
                              "~F21~(g) phase 1 day 15","~F21~(h) phase 1 day 20"/)
  wks = gsn_open_wks("eps", pdir+"/t_shading_u_contour")

  plott   = new(8,graphic)
  plotphi = plott
  dum1    = plott
  dum2    = plott
  dum     = plott

  respoly = True
  respoly@gsLineThicknessF = 4.0
  respoly@gsLineColor      = "black"

  do n = 0,7
    rest@gsnLeftString = "~F21~"+Label1(n)
    x = llt(n,:,:)
    x = mask(x,x.lt.m_valuet(0,0,n,:,:).or.x.gt.m_valuet(1,0,n,:,:),True)
    y  = llmtg(n,:,:)
    xt = llmtg(n,:,:)
    xt = 0.
    do i = 0, dimsizes(lat)-1
      do j = 0, dimsizes(lon)-1
        if(y(i,j).lt.m_valueu(0,0,n,i,j).or.y(i,j).gt.m_valueu(1,0,n,i,j)) then
          xt(i,j) = 2.
        else
          xt(i,j) = -2.
        end if
      end do
    end do
    plott(n)   = gsn_csm_contour_map(wks, x, rest)
    plotphi(n) = gsn_csm_contour(wks, llmtg(n,:,:), resu)
    dum1(n)  = gsn_add_polyline(wks, plott(n), (/60,80,80,60,60/), \
                                (/20,20,30,30,20/), respoly)
    dum2(n)  = gsn_add_polyline(wks, plott(n), (/160,190,190,160,160/), \
                                (/25,25,35,35,25/), respoly)
    plotu = gsn_csm_contour(wks, u, resc)
    dum(n)   = gsn_csm_contour(wks,xt,sres)
    dum(n)   = gsn_contour_shade(dum(n),1.9,-999.,opt)
    overlay(plott(n),dum(n))
    overlay(plott(n),plotu)
    overlay(plott(n),plotphi(n))
    delete(x)
  end do

  gsn_panel(wks, (/plott(0),plott(4),plott(1),plott(5),\
                plott(2),plott(6),plott(3),plott(7)/), (/2,2,2,2/), resP)
  delete([/plott,plotphi/])

  ;i = 1
  ;do while (i.eq.1)
  ;  i = 1
  ;end do
