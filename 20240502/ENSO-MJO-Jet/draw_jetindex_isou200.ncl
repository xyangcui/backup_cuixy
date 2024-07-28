; ==============================================================

  yrStrt = 19790101
  yrLast = 20221231

  pi = atan(1.0)*4.

  ddir = "/home/sunming/data5/cuixy/ENSO-MJO-Jet/data"
;--------------------------------------------------------------
;; Main code.
;--------------------------------------------------------------
  f  = addfile (ddir+"/ERA5_u200_daily_anom_1979-2022_r1.nc", "r")
  
  YYYY   = cd_calendar(f->time,-2)   
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast) 
  u  = f->uwnd(iYYYY,:,:)

  TIME = cd_calendar(u&time,0)
  year = floattointeger(TIME(:,0))
  mon  = floattointeger(TIME(:,1))
  day  = floattointeger(TIME(:,2))  
  win = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2.or.mon.eq.3)
;;--------------------------------------------------------------
;; select MJO date.
;--------------------------------------------------------------
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
  RMMf  = rmmdir
  delim = ","
  rmmf  = asciiread(RMMf,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
  phase = str_get_field(rmmf,6,delim)
  amp   = str_get_field(rmmf,7,delim)

  iyear := toint(iyear(1::))
  imon  := toint(imon(1::))
  iday  := toint(iday(1::))
  phase := toint(phase(1::))
  amp   := tofloat(amp(1::))

  reamp   = new((/dimsizes(u&time)/),float)
  rephase = new((/dimsizes(u&time)/),integer)
  reamp   = 0.
  rephase = 0

  do i = 0, dimsizes(u&time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
    reamp(i)   = amp(nt)
    rephase(i) = phase(nt)
    delete(nt)
  end do
  delete(iyear)
;----------------------------------------------------------------------

  um = dim_avg_n_Wrap(u(win,:,:),0)

  year := year(win)
  mon  := mon(win)

  u := u(win,:,:)

  reamp   := reamp(win)
  rephase := rephase(win)

  year(ind(mon.eq.12)) = year(ind(mon.eq.12))+1

;; read in jet index and select years.

  delim = ","
  jif   = asciiread(ddir+"/interannual_jet_index.csv",-1,"string")
  iyear = str_get_field(jif,1,delim)
  iME   = str_get_field(jif,2,delim)
  iEA   = str_get_field(jif,3,delim)
  iNA   = str_get_field(jif,4,delim)

  iyear := toint(iyear(1::))
  iME   := tofloat(iME(1::))
  iEA   := tofloat(iEA(1::))
  iNA   := tofloat(iNA(1::))

;; MEJS
  ny = 10
  ip = dim_pqsort(iME, 1)
  ne = dimsizes(ip)-1
  nb = dimsizes(ip)-ny
  pme = iyear(ip(nb:ne))
  nme = iyear(ip(0:ny-1))
  delete([/ip,ne,nb/])
;; EAJS
  ip = dim_pqsort(iEA, 1)
  ne = dimsizes(ip)-1
  nb = dimsizes(ip)-ny
  pea = iyear(ip(nb:ne))
  nea = iyear(ip(0:ny-1))
  delete([/ip,ne,nb/])
;; NAJS
  ip = dim_pqsort(iNA, 1)
  ne = dimsizes(ip)-1
  nb = dimsizes(ip)-ny
  pna = iyear(ip(nb:ne))
  nna = iyear(ip(0:ny-1))
  delete([/ip,ne,nb/])

;; for store zonal wind.
  ampt = 1.
  nvalue = new((/21/),float)
  nvalue = 0.
;; MEJS.
print("MEJS")

do k = 0, 20
  nva = pme
  N = dimsizes(nva)*90.
  upme = u(0,:,:)
  upme = 0.
  do j = 1,8
    nt = ind(year.eq.nva(0).and.reamp.gt.ampt.and.rephase.eq.j.and.mon.ne.3)+k
    do i = 1,dimsizes(nva)-1
      ni = ind(year.eq.nva(i).and.reamp.gt.ampt.and.rephase.eq.j.and.mon.ne.3)+k
      nt := array_append_record(nt, ni, 0)
      delete(ni)
    end do
    igood = ind(.not.ismissing(nt))
    ngood = nt(igood)
    upme = dim_avg_n_Wrap(u(ngood,:,:), 0)*dimsizes(ngood) + upme
    delete([/nt,igood,ngood/])
  end do
  upme = upme/N
  copy_VarCoords(u(0,:,:),upme)
  delete(nva)

  nva = nme
  N = dimsizes(nva)*90.
  unme = u(0,:,:)
  unme = 0.
  do j = 1,8
    nt = ind(year.eq.nva(0).and.reamp.gt.ampt.and.rephase.eq.j.and.mon.ne.3)+k
    do i = 1,dimsizes(nva)-1
      ni = ind(year.eq.nva(i).and.reamp.gt.ampt.and.rephase.eq.j.and.mon.ne.3)+k
      nt := array_append_record(nt, ni, 0)
      delete(ni)
    end do
    igood = ind(.not.ismissing(nt))
    ngood = nt(igood)
    unme = dim_avg_n_Wrap(u(ngood,:,:), 0)*dimsizes(ngood) + unme
    delete([/nt,igood,ngood/])
  end do

  unme = unme/N
  copy_VarCoords(u(0,:,:),unme)
  delete(nva)

  medf = upme - unme
  copy_VarCoords(upme, medf)

  u0 = medf({25:35},{45:75})
  clat = cos(u0&lat*pi/180.)
  umean = wgt_areaave(u0, clat, 1.0, 1)
  nvalue(k) = umean
  delete([/u0,clat,umean/])
end do
print(nvalue)

;; EAJS.
print("EAJS")
do k = 0,20
  nva = pea
  N = dimsizes(nva)*90.
  upea = u(0,:,:)
  upea = 0.
  do j = 1,8
    nt = ind(year.eq.nva(0).and.reamp.gt.ampt.and.rephase.eq.j.and.mon.ne.3)+k
    do i = 1,dimsizes(nva)-1
      ni = ind(year.eq.nva(i).and.reamp.gt.ampt.and.rephase.eq.j.and.mon.ne.3)+k
      nt := array_append_record(nt, ni, 0)
      delete(ni)
    end do
    igood = ind(.not.ismissing(nt))
    ngood = nt(igood)
    upea = dim_avg_n_Wrap(u(ngood,:,:), 0)*dimsizes(ngood) + upea
    delete([/nt,igood,ngood/])
  end do
  upea = upea/N
  copy_VarCoords(u(0,:,:),upea)
  delete(nva)

  nva = nea
  N = dimsizes(nva)*90.
  unea = u(0,:,:)
  unea = 0.
  do j = 1,8
    nt = ind(year.eq.nva(0).and.reamp.gt.ampt.and.rephase.eq.j.and.mon.ne.3)+k
    do i = 1,dimsizes(nva)-1
      ni = ind(year.eq.nva(i).and.reamp.gt.ampt.and.rephase.eq.j.and.mon.ne.3)+k
      nt := array_append_record(nt, ni, 0)
      delete(ni)
    end do
    igood = ind(.not.ismissing(nt))
    ngood = nt(igood)
    unea = dim_avg_n_Wrap(u(ngood,:,:), 0)*dimsizes(ngood) + unea
    delete([/nt,igood,ngood/])
  end do
  unea = unea/N
  copy_VarCoords(u(0,:,:),unea)
  delete(nva)

  eadf = upea - unea
  copy_VarCoords(upea, eadf)

  u0 = eadf({25:40},{170:210})
  clat = cos(u0&lat*pi/180.)
  umean = wgt_areaave(u0, clat, 1.0, 1)
  nvalue(k) = umean
  delete([/u0,clat,umean/])
end do
print(nvalue)
;; NAJS.
print("NAJS")
do k = 0, 20
  nva = pna
  N = dimsizes(nva)*90.
  upna = u(0,:,:)
  upna = 0.
  do j = 1,8
    nt = ind(year.eq.nva(0).and.reamp.gt.ampt.and.rephase.eq.j.and.mon.ne.3)+k
    do i = 1,dimsizes(nva)-1
      ni = ind(year.eq.nva(i).and.reamp.gt.ampt.and.rephase.eq.j.and.mon.ne.3)+k
      nt := array_append_record(nt, ni, 0)
      delete(ni)
    end do
    igood = ind(.not.ismissing(nt))
    ngood = nt(igood)
    upna = dim_avg_n_Wrap(u(ngood,:,:), 0)*dimsizes(ngood) + upna
    delete([/nt,igood,ngood/])
  end do
  upna = upna/N
  copy_VarCoords(u(0,:,:),upna)
  delete(nva)

  nva = nna
  N = dimsizes(nva)*90.
  unna = u(0,:,:)
  unna = 0.
  do j = 1,8
    nt = ind(year.eq.nva(0).and.reamp.gt.ampt.and.rephase.eq.j.and.mon.ne.3)+k
    do i = 1,dimsizes(nva)-1
      ni = ind(year.eq.nva(i).and.reamp.gt.ampt.and.rephase.eq.j.and.mon.ne.3)+k
      nt := array_append_record(nt, ni, 0)
      delete(ni)
    end do
    igood = ind(.not.ismissing(nt))
    ngood = nt(igood)
    unna = dim_avg_n_Wrap(u(ngood,:,:), 0)*dimsizes(ngood) + unna
    delete([/nt,igood,ngood/])
  end do
  unna = unna/N
  copy_VarCoords(u(0,:,:),unna)
  delete(nva)

  nadf = upna - unna
  copy_VarCoords(upna, nadf)

  u0 = nadf({30:55},{300:335})
  clat = cos(u0&lat*pi/180.)
  umean = wgt_areaave(u0, clat, 1.0, 1)
  nvalue(k) = umean
  delete([/u0,clat,umean/])
end do
print(nvalue)                       
;--plot
  nlats = 10
  nlatn = 70
  nlonw = 240
  nlone = 360

  mlats = 10
  mlatn = 45
  mlonw = 0
  mlone = 90

  elats = 10
  elatn = 50
  elonw = 100
  elone = 250

  res = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False
  res@gsnAddCyclic = True
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringFontHeightF = 0.015
  res@gsnLeftStringOrthogonalPosF  = 0.01
  res@gsnRightStringOrthogonalPosF = 0.01

  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01

  res@mpFillOn = True
  res@mpLandFillColor = "white"

  res@cnFillOn = True
  res@cnFillPalette        = "BlueWhiteOrangeRed"
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "Manuallevels"
  res@cnMinLevelValF       =  -2
  res@cnMaxLevelValF       =   2
  res@cnLevelSpacingF      =   0.2
  res@lbLabelBarOn = False

  res@lbBoxLinesOn          = True
  res@lbBoxSeparatorLinesOn = False
  res@pmLabelBarOrthogonalPosF = 0.2

  res@tmXBLabelFont = 21
  res@tmYRLabelFont = 21

  res@gsnLeftString    = "~F21~(a) ~F21~200hPa U"
  res@gsnRightString   = "~F21~DJF"

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
  resc@cnLineColor        = "red"
  resc@cnLineThicknessF   = 2.0
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels      = (/30.,40.,50.,60.,70./)

  resp = True
  resp@gsLineColor = "black"
  resp@gsLineThicknessF = 3.0
  resp@gsLineDashPattern = 0

  resP = True
  resP@gsnMaximize      = True
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02

  resP@pmLabelBarOrthogonalPosF = -0.02
  resP@pmLabelBarParallelPosF   = 0.0
  resP@pmLabelBarWidthF         = 0.5

  resP@lbLabelFont       = 21
  resP@lbTitleString     = "~F21~u anomaly (m/s)"
  resP@lbTitlePosition   = "Bottom"
  resP@lbTitleDirection  = "Across"
  resP@lbTitleFontHeightF= 0.012
  resP@lbBoxLinesOn      = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012


  wks = gsn_open_wks("x11","eof") 

  plot = new(9,graphic)

  res@mpMinLatF   = mlats
  res@mpMaxLatF   = mlatn
  res@mpMinLonF   = mlonw
  res@mpMaxLonF   = mlone

  plot(0) = gsn_csm_contour_map(wks, upme, res)
  plotu = gsn_csm_contour(wks, um, resc)
  dum1  = gsn_add_polyline(wks, plot(0), (/45,45,75,75,45/), (/25,35,35,25,25/), resp)
  overlay(plot(0), plotu)

  plot(1) = gsn_csm_contour_map(wks, unme, res)
  plotu = gsn_csm_contour(wks, um, resc)
  dum2  = gsn_add_polyline(wks, plot(1), (/45,45,75,75,45/), (/25,35,35,25,25/), resp)
  overlay(plot(1), plotu)

  plot(2) = gsn_csm_contour_map(wks, medf, res)
  plotu = gsn_csm_contour(wks, um, resc)
  dum21  = gsn_add_polyline(wks, plot(2), (/45,45,75,75,45/), (/25,35,35,25,25/), resp)
  overlay(plot(2), plotu)

  res@mpMinLatF   = elats
  res@mpMaxLatF   = elatn
  res@mpMinLonF   = elonw
  res@mpMaxLonF   = elone

  plot(3) = gsn_csm_contour_map(wks, upea, res)
  plotu = gsn_csm_contour(wks, um, resc)
  dum3  = gsn_add_polyline(wks, plot(3), (/170,170,210,210,170/), (/25,40,40,25,25/), resp)
  overlay(plot(3), plotu)

  plot(4) = gsn_csm_contour_map(wks, unea, res)
  plotu = gsn_csm_contour(wks, um, resc)
  dum4  = gsn_add_polyline(wks, plot(4), (/170,170,210,210,170/), (/25,40,40,25,25/), resp)
  overlay(plot(4), plotu)

  plot(5) = gsn_csm_contour_map(wks, eadf, res)
  plotu = gsn_csm_contour(wks, um, resc)
  dum41  = gsn_add_polyline(wks, plot(5), (/170,170,210,210,170/), (/25,40,40,25,25/), resp)
  overlay(plot(5), plotu)

  res@mpMinLatF   = nlats
  res@mpMaxLatF   = nlatn
  res@mpMinLonF   = nlonw
  res@mpMaxLonF   = nlone

  plot(6) = gsn_csm_contour_map(wks, upna, res)
  plotu = gsn_csm_contour(wks, um, resc)
  dum5  = gsn_add_polyline(wks, plot(6), (/-60,-60,-25,-25,-60/), (/30,55,55,30,30/), resp)
  overlay(plot(6), plotu)

  plot(7) = gsn_csm_contour_map(wks, unna, res)
  plotu = gsn_csm_contour(wks, um, resc)
  dum6  = gsn_add_polyline(wks, plot(7), (/-60,-60,-25,-25,-60/), (/30,55,55,30,30/), resp)
  overlay(plot(7), plotu)

  plot(8) = gsn_csm_contour_map(wks, nadf, res)
  plotu = gsn_csm_contour(wks, um, resc)
  dum61  = gsn_add_polyline(wks, plot(8), (/-60,-60,-25,-25,-60/), (/30,55,55,30,30/), resp)
  overlay(plot(8), plotu)

  gsn_panel(wks,plot,(/3,3,3/),resP)

  i = 0
  do while (i.eq.0)
    i = 0
  end do