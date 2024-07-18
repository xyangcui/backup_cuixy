;-----------------------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  udir = ddir + "/ERA5_daily_anom_ulevel_1979-2022_r2.5.nc"
  zdir = ddir + "/ERA5_daily_anom_zlevel_1979-2022_r2.5.nc"

  lats =  15
  latn =  75
  lonw = -50
  lone =  20

  pi = atan(1.0)*4. 
  nsvd = 2
;------------------------------------------------------------
; Main code.
;------------------------------------------------------------
  fu = addfile(udir,"r")
  fz = addfile(zdir,"r")
  u = fu->u(:,{250:500},{lats:latn},{lonw:lone})
  z = fz->z(:,{500},{lats:latn},{lonw:lone})

  time = u&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.4.and.day.eq.1)

  mon:= mon(n1:n2)
  u  := u(n1:n2,:,:,:)
  z  := z(n1:n2,:,:)
  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  u := u(nt,:,:,:)
  z := z(nt,:,:)
  delete([/n1,n2,nt/])

  ntime  = dimsizes(u&time)
  nlevel = dimsizes(u&level)
  nlat   = dimsizes(u&lat)
  nlon   = dimsizes(u&lon)

  nu = dimsizes(u&level)*dimsizes(u&lat)*dimsizes(u&lon)
  ut = reshape(u, (/ntime,nu/))
  nz = dimsizes(z&lat)*dimsizes(z&lon)
  zt = reshape(z, (/ntime,nz/))

  ut!0 = "time"
  ut!1 = "col"
  zt!0 = "time"
  zt!1 = "col"

  homlft = new((/nsvd,nz/),float)
  hetlft = new((/nsvd,nz/),float)
  homrgt = new((/nsvd,nu/),float)
  hetrgt = new((/nsvd,nu/),float)

  print("data reshape done.")
;;SVD
  x = svdcov(zt(col|:,time|:),ut(col|:,time|:),\
             nsvd,homlft,hetlft,homrgt,hetrgt) 

  print("SVD calc done.")
;;return: time series.
  zk   = onedtond(x@ak,(/nsvd,ntime/))
  uk   = onedtond(x@bk,(/nsvd,ntime/))
  zk!0 = "sv"
  zk!1 = "time"
  uk!0 = "sv"
  uk!1 = "time"

;;return: Hect array.
  zh = reshape(hetlft, (/nsvd,nlat,nlon/))
  uh = reshape(hetrgt, (/nsvd,nlevel,nlat,nlon/))
  copy_VarCoords(z(0,:,:), zh(0,:,:))
  copy_VarCoords(u(0,:,:,:), uh(0,:,:,:))

  print("return value calc done.")
;============================================================
; PLOTS
;============================================================
  wks = gsn_open_wks("x11","svd_olr")         
  plot = new(nsvd,graphic)                
  plotc= plot   
  ploto= plot                                     
; EOF patterns
  res                      = True         
  res@gsnDraw              = False        
  res@gsnFrame             = False        
  res@gsnAddCyclic         = False       
 
  res@mpFillOn             = False        
  res@mpMinLatF            = lats         
  res@mpMaxLatF            = latn
  res@mpMinLonF            = lonw
  res@mpMaxLonF            = lone

  res@cnFillOn             = True         
  res@cnLinesOn            = False        
  res@cnFillPalette        = "BlWhRe"     
  res@cnLevelSelectionMode   = "Manuallevels"
  res@cnMinLevelValF         = -0.5
  res@cnMaxLevelValF         =  0.5
  res@cnLevelSpacingF        =  0.05 

  res@lbLabelBarOn         = False  

  cmap1 = read_colormap_file("BlueWhiteOrangeRed")
  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = False
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 30
  rest@gsnRightString     = "~F21~OLR"
  rest@gsnLeftStringFontHeightF  = 0.02
  rest@gsnRightStringFontHeightF = 0.02
  rest@gsnRightStringOrthogonalPosF = 0.01
  rest@gsnLeftStringOrthogonalPosF  = 0.01

  rest@tmXTOn = False
  rest@tmYROn = False
  rest@tmYLLabelFontHeightF = 0.015
  rest@tmXBLabelFontHeightF = 0.015

  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF     = lats         
  rest@mpMaxLatF     = latn
  rest@mpMinLonF     = lonw
  rest@mpMaxLonF     = lone
  rest@mpLandFillColor = "white"

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap1;(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -0.5
  rest@cnMaxLevelValF         =  0.5
  rest@cnLevelSpacingF        =  0.05
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = False

  label1 = (/"(a) u250_SVD1","(b) u250_SVD2"/)
  label2 = (/"(a) OLR_SVD1","(b) OLR_SVD2"/)
  res@gsnLeftStringFontHeightF  = 0.02
  res@gsnRightStringFontHeightF = 0.02
  res@gsnRightStringOrthogonalPosF = 0.01
  res@gsnLeftStringOrthogonalPosF  = 0.01
  res@tmYLLabelFontHeightF = 0.015
  res@tmXBLabelFontHeightF = 0.015

  do n=0,nsvd-1
    res@gsnLeftString  = label1(n)
    res@gsnRightString = sprintf("%5.1f", x(n)) +"%"
    plot(n) = gsn_csm_contour_map(wks,uh(n,{250},:,:),res)
  end do

  do n=0,nsvd-1
    rest@gsnLeftString  = label2(n)
    rest@gsnRightString = sprintf("%5.1f", x(n)) +"%"
    ploto(n) = gsn_csm_contour_map(wks,zh(n,:,:),rest)
  end do

  resP = True
  resP@gsnMaximize      = False
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True

  resP@pmLabelBarWidthF = 0.6

  resP@lbTitleFontHeightF = 0.012
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012

  resP@gsnPanelLabelBar = True
  resP@lbTitleOn = False
  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.013
  resP@lbTitleStringParallelPosF = -0.1

  resP@gsnPanelTop    = 0.98
  resP@gsnPanelBottom = 0.02
  
  gsn_panel(wks,(/plot(0),ploto(0),plot(1),ploto(1)/)\
            ,(/2,2/),resP) 

  i = 0
  do while (i.eq.0)
    i = 0
  end do