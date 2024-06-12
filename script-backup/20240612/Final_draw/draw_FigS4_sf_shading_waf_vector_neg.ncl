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
  ucdir  = ddir+"/"+"ERA5_daily_u250_1979-2022_r2.5.nc.nc"
  sigdir = ddir+"/"+"monte_carlo/Fig9_sf_sig.nc"

  lats = -20
  latn =  90
  lonw = 100
  lone = 360
;--------------------------------------------------------
; u250
;--------------------------------------------------------
  f    = addfile(ucdir,"r")

  ut = f->u(:,0,:,:)
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)

  mon := mon(n1:n2)
  ut  := ut(n1:n2,:,:)
  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  uc = dim_avg_n_Wrap(ut(nt,:,:),0)

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

  delete([/time,f,TIME,mon,nt,ut,n1,n2,year,day/])

  print("draw clim.u complete. "+ systemfunc("date"))
;--------------------------------------------------------
;; for calc WAF.
  u1dir = ddir+"/"+"ERA5_daily_u250_1979-2022_r2.5.nc"
  vdir = ddir+"/"+"ERA5_daily_v250_1979-2022_r2.5.nc"

  fu = addfile(u1dir,"r")
  fv = addfile(vdir,"r")

  u0 = fu->u(:,{250},{lats:latn},{lonw:lone})
  v0 = fv->v(:,{250},{lats:latn},{lonw:lone})

  time = u0&time
  lat  = u0&lat
  lon  = u0&lon

  TIME = cd_calendar(time,0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  wgtl = filwgts_lanczos(201, 0, 1./100., -999., 1.)
  ;u1   = wgt_runave_n_Wrap(u0, wgtl, 0, 0)
  ;v1   = wgt_runave_n_Wrap(v0, wgtl, 0, 0)

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  fphi = addfile(phidir,"r")
  fu   = addfile(udir,"r")
  fv   = addfile(vdir,"r")
  folr = addfile(olrdir,"r")

  olr  = folr->olr(:,{lats:latn},{lonw:lone})
  phi  = fphi->phi(:,{lats:latn},{lonw:lone})
  u    = fu->u(:,{lats:latn},{lonw:lone})
  printVarSummary(phi)
  print("complete data selling "+ systemfunc("date"))

;; sig
  fsig = addfile(sigdir,"r")
  mvalue = fsig->mvalue
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
  u1 = u0(n1:n2,:,:)
  v1 = v0(n1:n2,:,:)

  nt     = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  um     = dim_avg_n_Wrap(u1(nt,:,:),0)
  vm     = dim_avg_n_Wrap(v1(nt,:,:),0)
  amp   := amp(nt)
  phase := phase(nt)

  delete(nt)
  ind_DJF = amp&timeind   ; indice of DJF between 1979-2022 DJF.
;----------------------------------------------------------------------
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
  llphi = new((/4,dimsizes(phi&lat),dimsizes(phi&lon)/),float)
  copy_VarCoords(phi(0,:,:), llphi(0,:,:))
  llphi = 0.

  llu = new((/4,dimsizes(u&lat),dimsizes(u&lon)/),float)
  copy_VarCoords(u(0,:,:), llu(0,:,:))
  llu = 0.

  llolr = new((/4,dimsizes(olr&lat),dimsizes(olr&lon)/),float)
  copy_VarCoords(olr(0,:,:), llolr(0,:,:))
  llolr = 0.

  lluw = new((/4,dimsizes(u&lat),dimsizes(u&lon)/),float)
  copy_VarCoords(u(0,:,:), lluw(0,:,:))
  lluw = 0.

  llvw = new((/4,dimsizes(u&lat),dimsizes(u&lon)/),float)
  copy_VarCoords(u(0,:,:), llvw(0,:,:))
  llvw = 0.

  prem = 25000./100000.
  nlati = ind(lat.ge.-20.and.lat.le.10)

  ip = 3
  nt = ind(phase.eq.8.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 3
    ntt   = -10 + j*5
    nindd = nind + ntt
    llphi(j,:,:) = dim_avg_n(phi(nindd,:,:), 0)
    llu(j,:,:)   = dim_avg_n(u(nindd,:,:), 0)
    llolr(j,:,:) = dim_avg_n(olr(nindd,:,:), 0)
    ;um = dim_avg_n(u0(nindd,:,:), 0)
    ;vm = dim_avg_n(v0(nindd,:,:), 0)
    WAFuv = cal_WFsf(llphi(j,:,:),um,vm,phi&lat,phi&lon)
    WAFuv(:,nlati,:) = 0.
    lluw(j,:,:) = WAFuv(0,:,:)*prem
    llvw(j,:,:) = WAFuv(1,:,:)*prem
    delete([/ntt,nindd,WAFuv/])
  end do
  delete([/ip,nt,nind/])

  llvalue = sqrt(lluw*lluw+llvw*llvw)

  lluw = mask(lluw,llvalue.gt.0.1,True)
  llvw = mask(llvw,llvalue.gt.0.1,True)

  llphi = llphi*10^(-6)

  llolr = smth9_Wrap(llolr, 0.5, 0.25, False)

  print("complete calc leadlag data. "+ systemfunc("date"))
;----------------------------------------------------------------------
;; resources
  cmap = read_colormap_file("BlueWhiteOrangeRed")
;shading u plot
  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = False
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 30
  rest@gsnRightStringFontHeightF = 0.02 
  rest@gsnLeftStringFontHeightF  = 0.02
  rest@gsnRightStringOrthogonalPosF = 0.01
  rest@gsnLeftStringOrthogonalPosF  = 0.01
  rest@gsnRightString = "~F21~250hPa_SF_WAF"


  rest@tmXTOn = False
  rest@tmYROn = False
  rest@tmXBLabelFont = 21
  rest@tmYLLabelFont = 21
  rest@tmXBLabelFontHeightF = 0.015
  rest@tmYLLabelFontHeightF = 0.015
  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF   = lats
  rest@mpMaxLatF   = latn-5
  rest@mpMinLonF   = lonw+10
  rest@mpMaxLonF   = lone-10
  rest@mpLandFillColor = "white"

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap;(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  ;rest@cnLevels = (/-4,-3.5,-3,-2.5,-2,-1.5,-1,0,1,1.5,2,2.5,3.5,4/)
  rest@cnMinLevelValF         = -6.
  rest@cnMaxLevelValF         =  6.
  rest@cnLevelSpacingF        =  1.
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = False

; contour olr plot.
  resolr = True
  resolr@gsnDraw      = False
  resolr@gsnMaximize  = False
  resolr@gsnFrame     = False
  resolr@gsnAddCyclic = False
  resolr@gsnContourNegLineDashPattern = 1
  resolr@cnMonoLineColor    = True
  resolr@cnLinesOn          = True
  resolr@cnLineLabelsOn     = False
  resolr@cnLineColor        = "red3"
  resolr@cnLineThicknessF   = 3.0
  resolr@cnInfoLabelPerimOn = False
  resolr@cnInfoLabelOn      = False
  resolr@cnLevelSelectionMode = "Explicitlevels"
  resolr@cnLevels      = (/-25,-20,-15,-10,-5,5,10,15,20,25/)

;; vector res.
  resw = True
  resw@gsnFrame     = False
  resw@gsnDraw      = False
  resw@gsnMaximize  = False
  resw@gsnAddCyclic = False
  resw@vcGlyphStyle            = "LineArrow"
  resw@vcLineArrowThicknessF   = 3.0
;resw@vcMinDistanceF          = 0.01
  resw@vcRefMagnitudeF         = 4
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

  resP = True
  resP@gsnMaximize      = True
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02
  ;resP@gsnPanelYWhiteSpacePercent = 3

  resP@pmLabelBarOrthogonalPosF = 0;-0.005
  resP@pmLabelBarWidthF         =  0.4
  resP@pmLabelBarParallelPosF   = -0.03;-0.1

  resP@lbLabelFont       = 21
  ;resP@lbTitleString     = "~F21~u anomaly (m/s)"
  ;resP@lbTitlePosition   = "Bottom"
  ;resP@lbTitleDirection  = "Across"
  ;resP@lbTitleFontHeightF= 0.01
  resP@lbBoxLinesOn      = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012
  print("complete res. "+ systemfunc("date"))

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
  sres@cnFillDotSizeF     = 0.0025
  sres@gsnLeftString      =" "
  sres@gsnRightString     =" "
  sres@cnMonoFillColor    = True
  sres@cnFillColor        = "gray"

  opt = True
  opt@gsnShadeFillType   = "pattern"
  opt@gsnShadeHigh       = 17
  opt@gsnShadeFillScaleF = 1.2

  respoly = True
  respoly@gsLineThicknessF  = 5.0
  respoly@gsLineDashPattern = 2
  respoly@gsLineColor       = "purple"

;----------------------------------------------------------------------
;; plot.
  wks = gsn_open_wks("x11", pdir+"/"+"Fig12_u_shading_waf_vector")

  plott   = new(4,graphic)
  plotolr = plott
  plotvc  = plott
  dum     = plott
  dum2    = plott
  Label1 = (/"~F21~(a) day -10","~F21~(b) day -5",\
             "~F21~(c) day 0","~F21~(d) day 5"/)

  do n = 0,3
    rest@gsnLeftString = "~F21~"+Label1(n)
    if (n.eq.3)then
      resw@vcRefAnnoOn             = True
      resw@vcRefAnnoString1On      = True
      resw@vcRefAnnoString1        = "$VMG$" +" ~F21~m~S~2~N~/s~S~2~N~"
      resw@vcRefAnnoOrthogonalPosF = 0.18
      resw@vcRefAnnoParallelPosF   = 1.0 
    end if
    x  = llphi(n,:,:)*10^(6)
    xt = llphi(n,:,:)
    xt = 0.
    do i = 0, dimsizes(lat)-1
      do j = 0, dimsizes(lon)-1
        if(x(i,j).lt.mvalue(1,0,n,i,j).or.x(i,j).gt.mvalue(1,1,n,i,j)) then
          xt(i,j) = 2.
        else
          xt(i,j) = -2.
        end if
      end do
    end do
    plott(n)   = gsn_csm_contour_map(wks, llphi(n,:,:), rest)
    plotolr(n) = gsn_csm_contour(wks, llolr(n,:,:), resolr)
    plotvc(n)  = gsn_csm_vector(wks, lluw(n,::3,::3), llvw(n,::3,::3), resw)
    plotu  = gsn_csm_contour(wks, uc, resc)
    dum(n)   = gsn_csm_contour(wks,xt,sres)
    dum(n)   = gsn_contour_shade(dum(n),1.9,-999.,opt)
    dum2(n)  = gsn_add_polyline(wks, plott(n), (/-110,-40,-40,-110,-110/), \
                                (/-90,-90,90,90,-90/), respoly)
    overlay(plott(n),dum2(n))
    overlay(plott(n),dum(n))
    overlay(plott(n),plotvc(n))
    overlay(plott(n),plotu)
    overlay(plott(n),plotolr(n))
  end do

  resP@lbTitleOn = True
  resP@lbTitleString = "[10~S~-6~N~ m~S~2~N~/s]"
  resP@lbTitleFontHeightF = 0.013
  resP@lbTitlePosition = "Bottom"
  gsn_panel(wks,(/plott(0),plott(1),plott(2),plott(3)/),(/1,1,1,1/),resP)

  i = 1
  do while (i.eq.1)
    i = 1
  end do
