;--------------------------------------------------------
; draw temperature anomaly in contour, MTG
; anomaly in shading and OLR anomaly in contour.
; time: 2023-11-22
;--------------------------------------------------------

;-----------------------------------------------------------

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  mtgdir = ddir+"/"+"ERA5_daily_anom_mtg_1979-2022_2.5x2.5.nc"
  tdir   = ddir+"/"+"ERA5_daily_anom_tvert_1979-2022_2.5x2.5.nc"
  udir   = ddir+"/"+"ERA5_daily_u250_1979-2022.nc"
  tsigdir = ddir+"/"+"monte_carlo/tvert_sig_phase1&5.nc"
  latmin = 0
  latmax = 70
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

  fmtg = addfile(mtgdir,"r")
  ft   = addfile(tdir,"r")

  mtg  = fmtg->mtg(:,{latmin:latmax},:)
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

  fsigt   = addfile(tsigdir,"r")
  m_value_t = fsigt->m_valuet

  m_valuet  = new((/2,2,8,dimsizes(t&lat),dimsizes(t&lon)/),float) 
  m_valuet  = 0.

  nl  = toint(dimsizes(m_value_t(:,0,0,0,0))*sigl) - 1
  nh  = toint(dimsizes(m_value_t(:,0,0,0,0))*sigh) - 1
  m_valuet(0,:,:,:,:) = m_value_t(nl,:,:,:,:)
  m_valuet(1,:,:,:,:) = m_value_t(nh,:,:,:,:)

  delete(m_value_t)

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

  llmtg = llmtg * 10^(6)
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
  rest@gsnLeftStringFontHeightF = 0.025

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

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -2.0
  rest@cnMaxLevelValF         =  2.0
  rest@cnLevelSpacingF        =  0.4
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = False

; contour phi plot.
  resphi = True
  resphi@gsnDraw      = False
  resphi@gsnMaximize  = False
  resphi@gsnFrame     = False
  resphi@gsnAddCyclic = True
  resphi@gsnContourNegLineDashPattern = 1
  resphi@cnMonoLineColor    = True
  resphi@cnLinesOn          = True
  resphi@cnLineLabelsOn     = False
  resphi@cnLineColor        = "black"
  resphi@cnLineThicknessF   = 2.5
  resphi@cnInfoLabelPerimOn = False
  resphi@cnInfoLabelOn      = False
  resphi@cnLevelSelectionMode = "Explicitlevels"
  resphi@cnLevels      = (/-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,\
                        4.0,3.5,3.0,2.5,2.0,1.5,1.0/)

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

  resP@pmLabelBarOrthogonalPosF = -0.03
  resP@pmLabelBarWidthF         = 0.7

  resP@lbTitleString     = "~F21~MTG anom (K/m)"
  resP@lbTitlePosition   = "Bottom"
  resP@lbTitleDirection  = "Across"
  resP@lbTitleFontHeightF= 0.015
  resP@lbBoxLinesOn      = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012
  print("complete res. "+ systemfunc("date"))

;; plot.
  Label1 = (/"~F21~(a) phase 1 day -15","~F21~(b) phase 1 day -10",\
                              "~F21~(c) phase 1 day -5","~F21~(d) phase 1 day 0",\
                              "~F21~(e) phase 1 day 5","~F21~(f) phase 1 day 10",\
                              "~F21~(g) phase 1 day 15","~F21~(h) phase 1 day 20"/)
  wks = gsn_open_wks("eps", pdir+"/mtg_shading_t_contour")

  plott   = new(8,graphic)
  plotphi = plott

  do n = 0,7
    rest@gsnLeftString = "~F21~"+Label1(n)
    x = llt(n,:,:)
    x = mask(x,x.lt.m_valuet(0,0,n,:,:).or.x.gt.m_valuet(1,0,n,:,:),True)
    plott(n)   = gsn_csm_contour_map(wks, llmtg(n,:,:), rest)
    plotphi(n) = gsn_csm_contour(wks, x, resphi)
    plotu = gsn_csm_contour(wks, u, resc)
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
