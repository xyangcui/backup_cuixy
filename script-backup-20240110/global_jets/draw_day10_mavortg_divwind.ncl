;--------------------------------------------------------
; draw temperature anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour.
; time: 2023-11-22
;--------------------------------------------------------
;-----------------------------------------------------------

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  phidir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022.nc"
  udir  = ddir+"/"+"ERA5_daily_divu_1979-2022_r2.5.nc"
  vdir  = ddir+"/"+"ERA5_daily_divv_1979-2022_r2.5.nc"
  vodir = ddir+"/"+"ERA5_daily_vor_1979-2022_r2.5.nc"
  olrdir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
  ucdir   = ddir+"/"+"ERA5_daily_u250_1979-2022.nc"
  idir  = ddir+"/"+"ERA5_daily_anom_rws_1979-2022.nc"

  latmin = -20
  latmax =  70
  lonmin =  0
  lonmax =  240
;--------------------------------------------------------
; u250
;--------------------------------------------------------
f    = addfile(ucdir,"r")

time = f->time
TIME = cd_calendar(time, 0)
mon  = toint(TIME(:,1))
nt   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)

uc    = dim_avg_n_Wrap(f->u(nt,0,:,:),0)

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

  print("Prepare data begin. "+systemfunc("date"))
  fu  = addfile(udir,"r")
  fv  = addfile(vdir,"r")
  fvo = addfile(vodir,"r")
  folr = addfile(olrdir,"r")

  olr  = folr->olr(time|:,lat|:,lon|:)
  vor= fvo->vor(:,{latmin:latmax},{lonmin:lonmax})

;; ISO filter.
  fi  = addfile(idir,"r")

  ui  = fi->u(:,{latmin:latmax},{lonmin:lonmax})
  vi  = fi->v(:,{latmin:latmax},{lonmin:lonmax})

  print("complete data selling "+ systemfunc("date"))

  time = ui&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

;-------------------------------------------------------
;--MJO date and DJF indices.
  print("sort MJO dates begin. "+systemfunc("date"))
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

  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.2.and.day.eq.28)

  amp   := amp(n1:n2)
  phase := phase(n1:n2)
  mon   := mon(n1:n2)
  time  := time(n1:n2)
  vor  := vor(n1:n2,:,:)

  nt     = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  amp   := amp(nt)
  phase := phase(nt)
  vorl  = dim_avg_n_Wrap(vor(nt,:,:),0)

  btemp = tofloat(beta_dfdy_rossby(vorl&lat,False))
  beta  = conform_dims(dimsizes(vorl), btemp, 0)
  copy_VarCoords(vorl, beta)

  grad_vorM = grad_latlon_cfd(vorl, vorl&lat,\
   vorl&lon, True, False)
  vorM_grad_lat = tofloat(grad_vorM[0])

  vorl := beta + vorM_grad_lat
  copy_VarCoords(ui(0,:,:), vorl)

  delete(nt)
  ind_DJF = amp&timeind   ; indice of DJF between 1979-2022 DJF.

  print("sort MJO dates done. "+systemfunc("date"))


;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
  llvor = new((/8,dimsizes(vor&lat),dimsizes(vor&lon)/),float)
  llvor = 0.

  llu = new((/8,dimsizes(ui&lat),dimsizes(ui&lon)/),float)
  llu = 0.

  llv = new((/8,dimsizes(vi&lat),dimsizes(vi&lon)/),float)
  llv = 0.

  llolr = new((/8,dimsizes(olr&lat),dimsizes(olr&lon)/),float)
  llolr = 0.

  ip = 0 + 1
  nt = ind(phase.eq.ip.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 7
    ntt   = -15 + j*5
    nindd = nind + ntt
    llvor(j,:,:) = vorl
    llu(j,:,:)   = dim_avg_n(ui(nindd,:,:), 0)
    llv(j,:,:)   = dim_avg_n(vi(nindd,:,:), 0)
    llolr(j,:,:) = dim_avg_n(olr(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do

  copy_VarCoords(vorl, llvor(0,:,:))
  copy_VarCoords(ui(0,:,:), llu(0,:,:))
  copy_VarCoords(vi(0,:,:), llv(0,:,:))
  copy_VarCoords(olr(0,:,:), llolr(0,:,:))
  ;llphi = llphi*10^(-6)

  print("complete calc leadlag data. "+ systemfunc("date"))

;; resources

;shading T plot
  cmap = read_colormap_file("MPL_RdBu")

  res = True
  res@gsnFrame           = False
  res@gsnDraw            = False
  res@gsnMaximize        = False
  res@gsnAddCyclic       = False
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30
  res@gsnLeftStringFontHeightF = 0.02 
  res@gsnRightStringFontHeightF = 0.02 
  res@gsnLeftString = "~F21~250hPa_div-wind_d~F33~z~F21~/dy"
  ;res@gsnRightString = "~F21~Phase1 day 10"

  res@mpFillOn        = False
  res@mpFillDrawOrder = "Postdraw"
  res@mpMinLatF   = latmin
  res@mpMaxLatF   = latmax
  res@mpMinLonF   = lonmin
  ;res@mpCenterLonF= 180
  res@mpMaxLonF   = lonmax
  res@mpLandFillColor = "white"

  res@cnFillOn               = True
  res@cnFillPalette          = "WhiteBlueGreenYellowRed"
  res@cnLinesOn              = False
  res@cnLineLabelsOn         = False
  res@cnLevelSelectionMode   = "Manuallevels"
  ;res@cnMinLevelValF         = -1.0
  ;res@cnMaxLevelValF         =  1.0
  ;res@cnLevelSpacingF        =  0.2
  res@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  res@lbLabelBarOn           = True
  res@lbTitleString          = "~F21~unit: (m s)~S~-1~N~"
  res@lbTitleFontHeightF     = 0.015
  res@lbTitlePosition    = "Bottom"

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
  resolr@cnLineColor        = "red4"
  resolr@cnLineThicknessF   = 2.5
  resolr@cnInfoLabelPerimOn = False
  resolr@cnInfoLabelOn      = False
  resolr@cnLevelSelectionMode = "Explicitlevels"
  resolr@cnLevels      = (/-20,-10,10,20/)

;; vector res.
  resw = True
  resw@gsnFrame     = False
  resw@gsnDraw      = False
  resw@gsnMaximize  = False
  resw@gsnAddCyclic = False
  resw@vcGlyphStyle            = "LineArrow"
  resw@vcLineArrowThicknessF   = 1.
;resw@vcMinDistanceF          = 0.01
  resw@vcRefMagnitudeF         = 2.5
  resw@vcRefLengthF            = 0.05
  resw@vcLineArrowColor        = "black"
  resw@vcRefAnnoString1On      = True
  resw@vcRefAnnoString2On      = True
  resw@vcRefAnnoString2        = ""
  resw@vcRefAnnoString1        = "$VMG$" +" m/s"
  resw@vcRefAnnoPerimSpaceF    = 0.3
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoFontHeightF    = 0.01
  resw@vcRefAnnoSide           = "Right"
  resw@vcRefAnnoJust           = "BottomRight"
  resw@vcRefAnnoOrthogonalPosF = -0.12
  resw@vcRefAnnoParallelPosF   = -0.28

  respoly = True
  respoly@gsLineThicknessF = 3.0
  respoly@gsLineColor      = "black"

  print("complete res. "+ systemfunc("date"))

;; plot.
  i  = 0
  ip = i + 1
  wks = gsn_open_wks("x11", pdir+"/"+"t-shading_phi_contour_olr_contour_phase"+ip)

  plott   = new(1,graphic)
  plotolr = plott
  plotvc  = plott

  plott   = gsn_csm_contour_map(wks, llvor(5,:,:), res)
  plotolr = gsn_csm_contour(wks, llolr(5,:,:), resolr)
  plotvc  = gsn_csm_vector(wks, llu(5,::2,::2), llv(5,::2,::2), resw)
  plotu = gsn_csm_contour(wks, uc, resc)
  dum1 = gsn_add_polyline(wks, plott,(/60,90,90,60,60/),\
    (/20,20,30,30,20/), respoly)
  dum2 = gsn_add_polyline(wks, plott,(/140,170,170,140,140/),\
    (/25,25,35,35,25/), respoly)
  overlay(plott,dum1)
  overlay(plott,dum2)
  overlay(plott,plotu)
  overlay(plott,plotolr)
  overlay(plott,plotvc)

  draw(plott)

  i = 1
  do while (i.eq.1)
    i = 1
  end do