;--------------------------------------------------------
; draw uwind anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour. with WAF.
; time: 2023-12-18
;--------------------------------------------------------
  lats = -10
  latn =  60
  lonw = -180
  lone =  20

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  olrdir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
  wdir = ddir+"/NAJS_temp_budget/ERA5_daily_anom_w500_1979-2022_r2.5.nc"

  sigdir = "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/NAJS_olr_uv_sig.nc"
  sigodir = "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/NAJS_omega_sig.nc"

  udir = ddir+"/NAJS_temp_budget/"+"ERA5_daily_iso_uvert_1979-2022_r2.5.nc"
  vdir = ddir+"/NAJS_temp_budget/"+"ERA5_daily_iso_vvert_1979-2022_r2.5.nc"
  tdir = ddir+"/NAJS_temp_budget/"+"ERA5_daily_low_tvert_1979-2022_r2.5.nc"

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  folr = addfile(olrdir,"r")

  olr = lonFlip(folr->olr)

  olr := olr(:,{lats:latn},{lonw:lone})

  lat = olr&lat
  lon = olr&lon

  time = olr&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  fw = addfile(wdir,"r")

  w  = fw->w(:,{lats:latn},{lonw:lone})

;; vdtdy
  ft = addfile(tdir,"r")
  t = ft->t(:,{lats:latn},{lonw:lone})

  fu = addfile(udir,"r")
  u  = fu->u(:,{lats:latn},{lonw:lone})

  fv = addfile(vdir,"r")
  v  = fv->v(:,{lats:latn},{lonw:lone})

  print("complete data selling "+ systemfunc("date"))

  printVarSummary(w)
  printVarSummary(olr)

;; sig
  fsig = addfile(sigdir,"r")
  mvalue = fsig->m_value
  fsigo  = addfile(sigodir,"r")
  momega = fsigo->m_value
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
  llolr1 = new((/5,dimsizes(olr&lat),dimsizes(olr&lon)/),float)
  llolr1 = 0.
  llw1 = new((/5,dimsizes(w&lat),dimsizes(w&lon)/),float)

  llt1 = new((/5,dimsizes(t&lat),dimsizes(t&lon)/),float)
  llu1 = new((/5,dimsizes(u&lat),dimsizes(u&lon)/),float)
  llv1 = new((/5,dimsizes(v&lat),dimsizes(v&lon)/),float)

  nt = ind(phase.eq.8.and.amp.gt.1.)
  nind = ind_DJF(nt)

  printVarSummary(w)
  printVarSummary(llw1)
  do j = 0, 4
    ntt   = 0 + j*1
    nindd = nind + ntt
    llolr1(j,:,:) = dim_avg_n(olr(nindd,:,:), 0)
    llw1(j,:,:) = dim_avg_n(w(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do
  copy_VarCoords(olr(0,:,:), llolr1(0,:,:))
  copy_VarCoords(w(0,:,:), llw1(0,:,:))

  llolr = dim_avg_n_Wrap(llolr1(0:5,:,:),0)
  ;llolr = smth9_Wrap(llolr, 0.5, 0.25, False)
  llw = dim_avg_n_Wrap(llw1(0:5,:,:),0)

  do j = 0, 4
    ntt   = -1 + j*1
    nindd = nind + ntt
    llt1(j,:,:) = dim_avg_n(t(nindd,:,:), 0)
    llu1(j,:,:) = dim_avg_n(u(nindd,:,:), 0)
    llv1(j,:,:) = dim_avg_n(v(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do

  llt = dim_avg_n(llt1(0:4,:,:),0)
  llu = dim_avg_n(llu1(0:4,:,:),0)
  llv = dim_avg_n(llv1(0:4,:,:),0)

  copy_VarCoords(t(0,:,:), llt)
  copy_VarCoords(u(0,:,:), llu)
  copy_VarCoords(v(0,:,:), llv)

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
  rest@mpMinLatF   =  lats
  rest@mpMaxLatF   =  latn
  rest@mpMinLonF   =  lonw
  rest@mpMaxLonF   =  lone
  rest@mpLandFillColor = "white"

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap;(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -0.024;-5
  rest@cnMaxLevelValF         =  0.024;5
  rest@cnLevelSpacingF        =  0.004;1
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = True
  rest@pmLabelBarWidthF       = 0.5
  rest@pmLabelBarOrthogonalPosF = 0.1

; contour plot for T.
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

;; contour plot for olr.
  reso = True
  reso@gsnDraw      = False
  reso@gsnMaximize  = False
  reso@gsnFrame     = False
  reso@gsnAddCyclic = False
  reso@cnMonoLineColor    = False
  reso@cnMonoLineDashPattern = 0
  reso@cnLinesOn          = True
  reso@cnLineLabelsOn     = False
  reso@cnLineThicknessF   = 3.0
  reso@cnInfoLabelPerimOn = False
  reso@cnInfoLabelOn      = False
  reso@cnLevelSelectionMode   = "Manuallevels"
  reso@cnMinLevelValF         = -15
  reso@cnMaxLevelValF         =  15
  reso@cnLevelSpacingF        =  3

;; vector res.
  resw = True
  resw@gsnFrame     = False
  resw@gsnDraw      = False
  resw@gsnMaximize  = False
  resw@gsnAddCyclic = False
  resw@gsnMajorLatSpacing = 15
  resw@gsnMajorLonSpacing = 30
  resw@gsnLeftStringFontHeightF  = 0.015 
  resw@gsnRightStringFontHeightF = 0.015
  resw@gsnLeftStringOrthogonalPosF  = 0.01
  resw@gsnRightStringOrthogonalPosF = 0.01

  resw@tmXTOn = False
  resw@tmYROn = False
  resw@tmYLLabelFontHeightF = 0.012
  resw@tmXBLabelFontHeightF = 0.012

  resw@mpFillOn        = False
  resw@mpFillDrawOrder = "Postdraw"
  resw@mpMinLatF   =  0
  resw@mpMaxLatF   =  60
  resw@mpMinLonF   = -110
  resw@mpMaxLonF   =  0
  resw@mpLandFillColor = "white"

  resw@vcGlyphStyle            = "CurlyVector"
  resw@vcLineArrowThicknessF   = 1.5
  resw@vcRefMagnitudeF         = 3.0
  resw@vcRefLengthF            = 0.05
  resw@vcLineArrowColor        = "gray64"
  resw@vcRefAnnoOn             = False
  resw@vcRefAnnoString2On      = False
  resw@vcRefAnnoString1On      = False
  resw@vcRefAnnoPerimSpaceF    = 0.4
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoFontHeightF    = 0.015
  resw@vcRefAnnoBackgroundColor = -1
  resw@vcRefAnnoString1        = "$VMG$" +" ~F21~m/s"


  reswt = True
  reswt@gsnFrame     = False
  reswt@gsnDraw      = False
  reswt@gsnMaximize  = False
  reswt@gsnAddCyclic = False

  reswt@vcGlyphStyle            = "CurlyVector"
  reswt@vcLineArrowThicknessF   = 2.5
  reswt@vcRefMagnitudeF         = 3.0
  reswt@vcRefLengthF            = 0.05
  reswt@vcLineArrowColor        = "black"
  reswt@vcRefAnnoOn             = True
  reswt@vcRefAnnoString2On      = False
  reswt@vcRefAnnoString1On      = True
  reswt@vcRefAnnoPerimSpaceF    = 0.4
  reswt@vcRefAnnoPerimOn        = False
  reswt@vcRefAnnoFontHeightF    = 0.015
  reswt@vcRefAnnoJust           = "BottomRight"
  reswt@vcRefAnnoBackgroundColor = -1
  reswt@vcRefAnnoString1        = "$VMG$" +" ~F21~m/s"
  reswt@vcRefAnnoOrthogonalPosF =  0.1
  reswt@vcRefAnnoParallelPosF   =  1.0

  respoly = True
  respoly@gsLineThicknessF = 3.5
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
  sres@cnFillDotSizeF     = 0.003
  sres@gsnLeftString      = " "
  sres@gsnRightString     = " "
  sres@cnMonoFillColor    = True
  sres@cnFillColor        = "gray" 

  opt = True
  opt@gsnShadeFillType   = "pattern"
  opt@gsnShadeHigh       = 17
  opt@gsnShadeFillScaleF = 1.2
;; plot.

  wks = gsn_open_wks("x11", pdir+"/"+"Fig7.NAJS_olr_omega500_UV")

  rest@vpXF = 0.25
  rest@vpYF = 0.95
  rest@vpWidthF  = 0.57
  rest@vpHeightF = 0.35
  rest@gsnLeftString = "~F21~(a) OLR_~F33~w~F21~500"
  rest@gsnRightString = "~F21~ ave: phase 8 day -10 to -5"
;; draw for omega.
  x = llw
  ;xt = x
  ;xt = 0.
  ;do i = 0, dimsizes(lat)-1
  ;  do j = 0, dimsizes(lon)-1
  ;    if(x(i,j).lt.momega(0,i,j).or.x(i,j).gt.momega(1,i,j)) then
   ;     xt(i,j) = 2.
   ;   else
   ;     xt(i,j) = -2.
   ;   end if
   ; end do
  ;end do
    
  ;dumw = gsn_csm_contour(wks,xt,sres)
  ;dumw = gsn_contour_shade(dumw,1.9,-999.,opt)

  plott = gsn_csm_contour_map(wks, llw, rest)
  ;overlay(plott,dumw)

;; draw for olr.
  delete(x)
  ploto = gsn_csm_contour(wks,llolr,reso)
  ploto = ColorNegDashZeroPosContour(ploto,"royalblue4","Transparent","red4")

  sres@cnFillColor = "black"
  opt@gsnShadeHigh = 4
  ;x = llolr
  ;xt = x
  ;xt = 0.
  ;do i = 0, dimsizes(lat)-1
  ;  do j = 0, dimsizes(lon)-1
  ;    if(x(i,j).lt.mvalue(0,0,i,j).or.x(i,j).gt.mvalue(1,0,i,j)) then
  ;      xt(i,j) = 2.
  ;    else
  ;      xt(i,j) = -2.
  ;    end if
   ; end do
  ;end do
    
  ;dumolr = gsn_csm_contour(wks,xt,sres)
  ;dumolr = gsn_contour_shade(dumolr,1.9,-999.,opt)
  ;overlay(ploto,dumolr)

  overlay(plott,ploto)
  dum1  = gsn_add_polyline(wks, plott, (/-70,-40,-40,-70,-70/), \
                          (/20,20,28,28,20/), respoly)
  draw(plott)

;; draw uv.
  ;delete(x)

  ;resw@vpXF = 0.25
  ;resw@vpYF = 0.48
  ;resw@vpWidthF  = 0.57
  ;resw@vpHeightF = 0.35
  ;resw@gsnLeftString = "~F21~(b) UV(ISO)_T(LOW)"
  ;resw@gsnRightString = "~F21~ ave: phase 8 day -1 to 4"

  ;x = llu
  ;y = llv

  ;x = mask(x,llu.lt.mvalue(0,1,:,:).or.llu.gt.mvalue(1,1,:,:).or.\
  ;        llv.lt.mvalue(0,2,:,:).or.llv.gt.mvalue(1,2,:,:),True)

  ;y = mask(y,llv.lt.mvalue(0,2,:,:).or.llv.gt.mvalue(1,2,:,:).or.\
  ;        llu.lt.mvalue(0,1,:,:).or.llu.gt.mvalue(1,1,:,:),True)


  ;plotuv  = gsn_csm_vector_map(wks, llu(::2,::2), llv(::2,::2), resw)
  ;plotuvt = gsn_csm_vector(wks, x(::2,::2), y(::2,::2), reswt)
  ;plotvt = gsn_csm_contour(wks, llt, resphi)

  ;overlay(plotuv,plotvt)
  ;overlay(plotuv,plotuvt)

  ;respoly@gsLineDashPattern = 1
  ;dum2  = gsn_add_polyline(wks, plotuv, (/-100,-70,-70,-100,-100/), \
  ;                        (/33,33,45,45,33/), respoly)

  ;restext = True
  ;restext@txFontColor = "green4"
  ;restext@txFontHeightF = 0.02
  ;dumtext = gsn_add_text(wks, plotuv, "~F22~C", -70, 45, restext)

  ;draw(plotuv)

  i = 1
  do while (i.eq.1)
    i = 1
  end do
