load "/home/sunming/data5/cuixy/script/ncl_func/cal_WFsf.ncl"
;--------------------------------------------------------
; draw uwind anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour. 
; waf vector
; This is for test sig. 95% by monte carlo method.
; time: 2023-12-18
; 2.5x2.5.
;--------------------------------------------------------

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  phidir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022_2.5x2.5.nc"
  udir   = ddir+"/"+"ERA5_daily_anom_u250_1979-2022_2.5x2.5.nc"
  olrdir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
  ucdir  = ddir+"/"+"ERA5_daily_u250_1979-2022.nc"
  ;sigdir = ddir+"/"+"monte_carlo/phiu_sig_phase15.nc"
  usigdir = ddir+"/"+"monte_carlo/u_sig_phase15.nc"
  phisigdir = ddir+"/"+"monte_carlo/phi_sig_phase15.nc"

  latmin = -30
  latmax = 90
;--------------------------------------------------------

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  fphi = addfile(phidir,"r")
  fu   = addfile(udir,"r")
  folr = addfile(olrdir,"r")

  olr  = folr->olr(time|:,lat|:,lon|:)
  phi  = fphi->phi(:,{latmax:latmin},:)
  u    = fu->u(:,{latmax:latmin},:)
  time = u&time
  
  TIME = cd_calendar(time,0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
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
;; sig
  sigl = 0.05
  sigh = 0.95

  fsigphi = addfile(phisigdir,"r")
  fsigu   = addfile(usigdir,"r")
  m_value_phi = fsigphi->m_valuephi(:,:,:,::-1,:)
  m_value_u   = fsigu->m_valueu(:,:,:,::-1,:)

  m_valuephi  = new((/2,2,8,dimsizes(phi&lat),dimsizes(phi&lon)/),float) 
  m_valuephi  = 0.
  m_valueu  = new((/2,2,8,dimsizes(u&lat),dimsizes(u&lon)/),float) 
  m_valueu  = 0.

  nl  = toint(dimsizes(m_value_phi(:,0,0,0,0))*sigl) - 1
  nh  = toint(dimsizes(m_value_phi(:,0,0,0,0))*sigh) - 1
  m_valuephi(0,:,:,:,:) = m_value_phi(nl,:,:,:,:)
  m_valuephi(1,:,:,:,:) = m_value_phi(nh,:,:,:,:)

  m_valuephi = m_valuephi*10^(-6)

  nl = toint(dimsizes(m_value_u(:,0,0,0,0))*sigl) - 1
  nh = toint(dimsizes(m_value_u(:,0,0,0,0))*sigh) - 1
  m_valueu(0,:,:,:,:) = m_value_u(nl,:,:,:,:)
  m_valueu(1,:,:,:,:) = m_value_u(nh,:,:,:,:)

  delete([/m_value_phi,m_value_u/])
;----------------------------------------------------------------------
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)

  llu = new((/8,dimsizes(u&lat),dimsizes(u&lon)/),float)
  copy_VarCoords(u(0,:,:), llu(0,:,:))
  llu = 0.

    ip = 1
    nt = ind(phase.eq.ip.and.amp.gt.1.)
    nind = ind_DJF(nt)
    do j = 0, 7
      ntt   = -15 + j*5
      nindd = nind + ntt
      llu(j,:,:)   = dim_avg_n(u(nindd,:,:), 0)
      delete([/ntt,nindd/])
    end do
    delete([/ip,nt,nind/])

  print("complete calc leadlag data. "+ systemfunc("date"))
;----------------------------------------------------------------------
;; resources
  cmap = read_colormap_file("MPL_RdBu")
;shading u plot
  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = True
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 30
  rest@gsnRightStringFontHeightF = 0.025 
  rest@gsnLeftStringFontHeightF = 0.025
  rest@gsnRightString = "~F21~psi: 10~S~6~N~m~S~2~N~/s"
  rest@tmXTOn = False
  rest@tmYROn = False
  rest@tmXBLabelFont = 21
  rest@tmYLLabelFont = 21
  rest@tmXBLabelFontHeightF = 0.015
  rest@tmYLLabelFontHeightF = 0.015
  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
;resolr@mpProjection= "Robinson"
  rest@mpMinLatF   = -20
  rest@mpMaxLatF   = 80
  ;rest@mpMinLonF   = 0
  ;rest@mpCenterLonF= 180
  ;rest@mpMaxLonF   = 360
  rest@mpLandFillColor = "white"
  rest@mpGridAndLimbDrawOrder = "PostDraw"
  rest@mpGridLineColor        = "Gray60"
;rest@mpGridLineThicknessF   = 1.0

  rest@cnFillOn               = True
;rest@cnFillMode             = "CellFill"
  rest@cnFillPalette          = cmap(::-1,:)
;rest@cnFillOpacityF         = 0.6
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -8.
  rest@cnMaxLevelValF         =  8.
  rest@cnLevelSpacingF        =  1.
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = False

  resP = True
  resP@gsnMaximize      = True
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02

  resP@pmLabelBarOrthogonalPosF = -0.02
  resP@pmLabelBarWidthF         = 0.75

  resP@lbLabelFont       = 21
  resP@lbTitleString     = "~F21~u anomaly (m/s)"
  resP@lbTitlePosition   = "Bottom"
  resP@lbTitleDirection  = "Across"
  resP@lbTitleFontHeightF= 0.012
  resP@lbBoxLinesOn      = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012
  print("complete res. "+ systemfunc("date"))

  sres   = True 
  sres@cnLineDrawOrder    = "PostDraw" 
  sres@gsnAddCyclic       = True                   
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
  ;sres@cnLevelSelectionMode   = "Manuallevels"
  ;sres@cnMinLevelValF         = -8.
  ;sres@cnMaxLevelValF         =  8.
  ;sres@cnLevelSpacingF        =  1.

  opt = True
  opt@gsnShadeFillType   = "pattern"
  ;opt@gsnShadeLow       = 10
  opt@gsnShadeHigh      = 10
  ;opt@gsnShadeMid        = 10
  opt@gsnShadeFillScaleF   = 1.
  ;opt@gsnShadeFillDotSizeF = 0.003
;----------------------------------------------------------------------
;; plot.
  wks = gsn_open_wks("x11", pdir+"/"+"u-shading_phi_contour_olr_contour")

  plott   = new(8,graphic)
  dum     = plott
  Label1 = (/"~F21~(a) phase 1 day -15","~F21~(b) phase 1 day -10",\
                              "~F21~(c) phase 1 day -5","~F21~(d) phase 1 day 0",\
                              "~F21~(e) phase 1 day 5","~F21~(f) phase 1 day 10",\
                              "~F21~(g) phase 1 day 15","~F21~(h) phase 1 day 20"/)

  lat = u&lat
  lon = u&lon
  do n = 0,7
    rest@mpMinLonF   = -10
    rest@mpMaxLonF   = 220
    rest@gsnLeftString = "~F21~"+Label1(n)
    x  = llu(n,:,:)
    xt = llu(n,:,:)
    xt = 0.
    do i = 0, dimsizes(lat)-1
      do j = 0, dimsizes(lon)-1
        if(x(i,j).lt.m_valueu(0,0,n,i,j).or.x(i,j).gt.m_valueu(1,0,n,i,j)) then
          xt(i,j) = 2.
        else
          xt(i,j) = -2.
        end if
      end do
    end do
    ;xt = mask(x,x.lt.m_valueu(0,0,n,:,:).or.x.gt.m_valueu(1,0,n,:,:),True)
    ;copy_VarCoords(x,xt)
    plott(n) = gsn_csm_contour_map(wks, x, rest)
    dum(n)   = gsn_csm_contour(wks,xt,sres)
    dum(n)   = gsn_contour_shade(dum(n),1.9,-999.,opt)
    overlay(plott(n),dum(n))
  end do
  gsn_panel(wks, (/plott(0),plott(4),plott(1),plott(5),\
                plott(2),plott(6),plott(3),plott(7)/), (/2,2,2,2/), resP)
  i = 0
  do while (i .eq. 0)
     i = 0
  end do
