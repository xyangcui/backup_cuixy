;--------------------------------------------------
; draw 20-70 days filtered OLR std.
; NE, EI and LA.
;--------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/ENSO-MJO-Jet/data"
  odir = ddir+"/olr_daily_anom_1979-2022.nc"
  ninodir = "/home/sunming/data5/cuixy/data/clim_indices/Nino_3.4_index.csv"

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------
  f = addfile(odir,"r")
  olr = f->olr

  time = olr&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
;----------------------------------------------------------
; prepare data for regression. (select DJF).
  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)
  year := year(n1:n2)
  mon  := mon(n1:n2)
  olr  := olr(n1:n2,:,:)

  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  year := year(nt)
  mon  := mon(nt)
  olr  := olr(nt,:,:)

  ni = ind(mon.eq.12)
  year(ni) = year(ni)+1
;; select for EI, LA or NE. Dec year + 1.
  delim  = ","
  ninof  = asciiread(ninodir,-1,"string")
  iyear  = str_get_field(ninof,1,delim)
  ievent = str_get_field(ninof,6,delim)

  iyear  := toint(iyear(1::))
  ievent := tofloat(ievent(1::))
;;EI
  nei = ind(ievent.eq.1)
  yei = iyear(nei)
;;LA
  nla = ind(ievent.eq.-1)
  yla = iyear(nla)
;;NE
  nne = ind(ievent.eq.0)
  yne = iyear(nne)

;; here for EI.
  ny = yne
  nc = ind(year.eq.ny(0))
  do i = 1, dimsizes(ny)-1
    ntemp = ind(year.eq.ny(i))
    nc := array_append_record(nc, ntemp, 0)
    delete(ntemp)
  end do
  nestd = dim_stddev_n_Wrap(olr(nc,:,:), 0)
  delete([/ny,nc/])
;; here for EI.
  ny = yei
  nc = ind(year.eq.ny(0))
  do i = 1, dimsizes(ny)-1
    ntemp = ind(year.eq.ny(i))
    nc := array_append_record(nc, ntemp, 0)
    delete(ntemp)
  end do
  eistd = dim_stddev_n_Wrap(olr(nc,:,:), 0)
  delete([/ny,nc/])
;; here for EI.
  ny = yla
  nc = ind(year.eq.ny(0))
  do i = 1, dimsizes(ny)-1
    ntemp = ind(year.eq.ny(i))
    nc := array_append_record(nc, ntemp, 0)
    delete(ntemp)
  end do
  lastd = dim_stddev_n_Wrap(olr(nc,:,:), 0)
  delete([/ny,nc/])

  estd = eistd-nestd
  lstd = lastd-nestd

  copy_VarCoords(olr(0,:,:), estd)
  copy_VarCoords(olr(0,:,:), lstd)
;;plot
  res = True
  res = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False
  res@gsnAddCyclic = True
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30

  res@gsnLeftStringOrthogonalPosF  = 0.012
  res@gsnRightStringOrthogonalPosF = 0.012
  res@gsnLeftStringFontHeightF  = 0.016
  res@gsnRightStringFontHeightF = 0.016
  res@tmXBLabelFontHeightF  = 0.012
  res@tmYLLabelFontHeightF  = 0.012

  res@mpFillOn = True
  res@mpLandFillColor = "white"
  ;res@mpCenterLonF    = 180.
  res@mpMinLatF       = -30
  res@mpMaxLatF       =  20
  res@mpMinLonF       =  30
  res@mpMaxLonF       =  180

  res@cnFillOn = True
  res@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
  res@cnFillPalette        = "MPL_coolwarm";"WhiteBlueGreenYellowRed"
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "Manuallevels"
  res@cnMinLevelValF       =  0.
  res@cnMaxLevelValF       =  3.5
  res@cnLevelSpacingF      =  0.25
  res@cnLevelSelectionMode = "Manuallevels"
  res@cnMinLevelValF       =  -5
  res@cnMaxLevelValF       =   5
  res@cnLevelSpacingF      =   1
  res@lbLabelBarOn = False
  res@pmLabelBarOrthogonalPosF = 0.25

  res@tmXBLabelFont = 21
  res@tmYRLabelFont = 21
  res@gsnRightString = ""  

  wks  = gsn_open_wks("x11","maponly")
;; NE
  res@gsnLeftString  = "~F21~(a) EI nino"
  res@vpYF = 0.95
  plot1  = gsn_csm_contour_map(wks, estd, res) 
;; EI
  res@lbLabelBarOn = True
  res@gsnLeftString  = "~F21~(b) LA nina"
  res@vpYF = 0.5
  plot2  = gsn_csm_contour_map(wks, lstd, res)

  draw(plot1)  
  draw(plot2) 

  i = 0
  do while (i.eq.0)
    i = 0
  end do