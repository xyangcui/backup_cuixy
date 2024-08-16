;----------------------------------------
;; select Jet events by PC1.
;; "candidate": consecutive 1 pentad 1.0.
;; Pre 12 days must grow and decay.
;; remove each year pre 12 days and late.
;----------------------------------------
  lats =  15
  latn =  80
  lonw = -90
  lone =  30

  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  indir = ddir+"/mhf_events.csv";"/centered_date.csv"
  udir = ddir+"/ERA5_daily_high_"+"u"+"250_1979-2022.nc"
  vdir = ddir+"/ERA5_daily_high_"+"v"+"250_1979-2022.nc"
  pdir = "/home/sunming/data5/cuixy/global_jets/plot/Fig10_E_vector"
; time table.
  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  fu   = addfile(udir,"r")
  fv   = addfile(vdir,"r")

  u = fu->u(:,{lats:latn},{lonw:lone})
  v = fv->v(:,{lats:latn},{lonw:lone})

  idir = ddir + "/ERA5_daily_anom_ulevel_1979-2022_r2.5.nc"
  f = addfile(idir,"r")
  ut = f->u(:,{250},{lats:latn},{lonw:lone})


  time = fu->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

;; selling events.
  delim  = ","
  rmmf   = asciiread(indir,-1,"string")
  iyear  = str_get_field(rmmf,1,delim)
  imon   = str_get_field(rmmf,2,delim)
  iday   = str_get_field(rmmf,3,delim)
  pct = str_get_field(rmmf,4,delim)

  iyear  := toint(iyear(1::))
  imon   := toint(imon(1::))
  iday   := toint(iday(1::))
  pct := toint(pct(1::))

  pev   = ind(pct.gt.0)
  pyear = iyear(pev)
  pmon  = imon(pev)
  pday  = iday(pev)

  pu = new((/dimsizes(pev),11,dimsizes(u&lat),dimsizes(u&lon)/),float)
  pu = 0.
  pv = pu
  put = pu

  do i = 0, dimsizes(pev)-1
    nt = ind(year.eq.pyear(i).and.mon.eq.pmon(i).and.day.eq.pday(i))
    pu(i,:,:,:) = u(nt-4:nt+6,:,:)
    pv(i,:,:,:) = v(nt-4:nt+6,:,:)
    put(i,:,:,:) = ut(nt-4:nt+6,:,:)
    delete(nt)
  end do

  nev   = ind(pct.lt.0)
  nyear = iyear(nev)
  nmon  = imon(nev)
  nday  = iday(nev)

  nu = new((/dimsizes(nev),11,dimsizes(u&lat),dimsizes(u&lon)/),float)
  nu = 0.
  nv = nu
  nut = nu

  do i = 0, dimsizes(nev)-1
    nt = ind(year.eq.nyear(i).and.mon.eq.nmon(i).and.day.eq.nday(i))
    nu(i,:,:,:) = u(nt-4:nt+6,:,:)
    nv(i,:,:,:) = v(nt-4:nt+6,:,:)
    nut(i,:,:,:) = ut(nt-4:nt+6,:,:)
    delete(nt)
  end do

  ;pu := dim_avg_n_Wrap(pu, 0)
  ;nu := dim_avg_n_Wrap(nu, 0)
  ;pv := dim_avg_n_Wrap(pv, 0)
  ;nv := dim_avg_n_Wrap(nv, 0)

  ;pu := dim_avg_n_Wrap(pu(6:8,:,:), 0)
  ;nu := dim_avg_n_Wrap(nu(6:8,:,:), 0)
  ;pv := dim_avg_n_Wrap(pv(6:8,:,:), 0)
  ;nv := dim_avg_n_Wrap(nv(6:8,:,:), 0)

  lat  = pu&lat
  rad  = atan(1.0)*4./180.
  cosl = cos(rad*lat)

  pEx = (pv*pv - pu*pu)*0.5;*conform(pv, cosl, 0)
  pEy = -1*pu*pv;*conform(pv, cosl, 0)

  nEx = (nv*nv - nu*nu)*0.5;*conform(pv, cosl, 0)
  nEy = -1*nu*nv;*conform(pv, cosl, 0)

  pEx := dim_avg_n_Wrap(pEx, 0)
  nEx := dim_avg_n_Wrap(nEx, 0)
  pEy := dim_avg_n_Wrap(pEy, 0)
  nEy := dim_avg_n_Wrap(nEy, 0)

  put := dim_avg_n_Wrap(put, 0)
  nut := dim_avg_n_Wrap(nut, 0)

  delete([/pu,pv,nu,nv/])
;;plot
  wks = gsn_open_wks("eps",pdir)         
  plot  = new(3,graphic)                
  plotc = plot  
  plotu = plot    
  plote = plot 

  cmap = read_colormap_file("MPL_RdBu") 
  res              = True         
  res@gsnDraw      = False        
  res@gsnFrame     = False       
  res@gsnAddCyclic = False  
  res@gsnRightStringFontHeightF  = 0.02 
  res@gsnLeftStringFontHeightF   = 0.02
  res@gsnCenterStringFontHeightF = 0.02
  res@gsnCenterStringOrthogonalPosF = 0.01
  res@gsnRightStringOrthogonalPosF  = 0.01
  res@gsnLeftStringOrthogonalPosF   = 0.01
  res@tmYLLabelFontHeightF = 0.015   
  res@tmXBLabelFontHeightF = 0.015   
 
  res@mpFillOn     = False        
  res@mpMinLatF    = lats        
  res@mpMaxLatF    = latn
  res@mpMinLonF    = lonw
  res@mpMaxLonF    = lone
  res@cnFillOn     = True                
  res@cnFillPalette = cmap(::-1,:) 
  res@cnLinesOn      = False
  res@cnLineLabelsOn = False
  ;res@cnLevelSelectionMode   = "Manuallevels"
  ;res@cnMinLevelValF         = -1.
  ;res@cnMaxLevelValF         =  1.
  ;res@cnLevelSpacingF        =  0.1                                       
  res@lbLabelBarOn = False

  resc              = True         
  resc@gsnDraw      = False        
  resc@gsnFrame     = False       
  resc@gsnAddCyclic = False 

  resc@cnLineLabelsOn = False
  resc@cnLineThicknessF = 3.5
  resc@cnLineColor      = "red"
  resc@cnInfoLabelOn    = False
  resc@gsnContourNegLineDashPattern = 1
  resc@gsnContourZeroLineThicknessF = 0
  ;res@cnLevelSelectionMode   = "Manuallevels"
  ;res@cnMinLevelValF         = -1.
  ;res@cnMaxLevelValF         =  1.
  ;res@cnLevelSpacingF        =  0.1 

  symMinMaxPlt(put, 5, False, resc)  

  resw = True
  resw@gsnFrame     = False
  resw@gsnDraw      = False
  resw@gsnMaximize  = False
  resw@gsnAddCyclic = False
  resw@vcGlyphStyle            = "LineArrow"
  resw@vcLineArrowThicknessF   = 2.5
  resw@vcRefMagnitudeF         = 100.
  resw@vcRefLengthF            = 0.05
  resw@vcLineArrowColor        = "mediumblue"
  resw@vcRefAnnoString2On      = False
  resw@vcRefAnnoString1        = "~F21~$VMG$" +" m s~S~-1~N~"
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoBackgroundColor = -1
  resw@vcRefAnnoOn = False

  resP                    = True         
  resP@gsnMaximize        = True         
  resP@gsnPanelLabelBar   = True 
  resP@gsnPanelYWhiteSpacePercent = 3

  resP@lbTitleString = "[~F21~10~S~-5~N~ m s~S~-2~N~]"   
  resP@lbTitlePosition  = "Right" 
  resP@lbTitleDirection = "Across" 
  resP@lbTitleFontHeightF = 0.015
  resP@pmLabelBarParallelPosF = 0.05
;*******************************************
; first plot
;*******************************************
  
  label = (/"~F21~(a) Day -4 to -2","(b) Day -1 to 1","(c) Day 2 to 4"/)
  i = 0
;; day -4 -2
  nb = 0
  ne = 2
  pEx1 = dim_avg_n_Wrap(pEx(nb:ne,:,:), 0)
  nEx1 = dim_avg_n_Wrap(nEx(nb:ne,:,:), 0)
  pEy1 = dim_avg_n_Wrap(pEy(nb:ne,:,:), 0)
  nEy1 = dim_avg_n_Wrap(nEy(nb:ne,:,:), 0)
  copy_VarCoords(u(0,:,:), pEx1(:,:))
  copy_VarCoords(u(0,:,:), pEy1(:,:))
  copy_VarCoords(u(0,:,:), nEx1(:,:))
  copy_VarCoords(u(0,:,:), nEy1(:,:))
  pdiv1 = uv2dv_cfd(pEx1, pEy1, u&lat, u&lon, 0)
  nEx1 = pEx1 - nEx1
  nEy1 = pEy1 - nEy1
  ndiv1 = uv2dv_cfd(nEx1, nEy1, u&lat, u&lon, 0)*10^(5)

  pE1 = sqrt(pEx1*pEx1+pEy1*pEy1)
  nE1 = sqrt(nEx1*nEx1+nEy1*nEy1)

  pu = dim_avg_n_Wrap(put(nb:ne,:,:), 0)
  nu = dim_avg_n_Wrap(nut(nb:ne,:,:), 0)
  pp = pu - nu
  pp = smth9_Wrap(pp, 0.5, 0.25, False)
  copy_VarCoords(pu,pp)

  ;pEx1 = mask(pEx1, pE1.gt.10., True)
  ;pEy1 = mask(pEx1, pE1.gt.10., True)
  ;nEx1 = mask(pEx1, nE1.gt.10., True)
  ;nEy1 = mask(pEx1, nE1.gt.10., True)
  copy_VarCoords(u(0,:,:),pdiv1(:,:))
  copy_VarCoords(u(0,:,:),ndiv1(:,:))

  symMinMaxPlt(ndiv1, 13, False, res) 
  res@gsnLeftString = label(i)
  res@gsnRightString = "~F21~High-Low"
  res@gsnCenterString = "~F22~E~F21~&div"
  plotu(i) = gsn_csm_vector(wks, nEx1(::2,::2), nEy1(::2,::2), resw)
  plot(i)  = gsn_csm_contour_map(wks,ndiv1(:,:),res)
  plote(i) = gsn_csm_contour(wks, pp, resc)
  overlay(plot(i), plotu(i))
  overlay(plot(i), plote(i))
  delete([/pEx1,pEy1,nEx1,nEy1,pE1,nE1,ndiv1,pdiv1/])

  i = 1
;; day -1 1
  nb = 3
  ne = 5
  pEx1 = dim_avg_n_Wrap(pEx(nb:ne,:,:), 0)
  nEx1 = dim_avg_n_Wrap(nEx(nb:ne,:,:), 0)
  pEy1 = dim_avg_n_Wrap(pEy(nb:ne,:,:), 0)
  nEy1 = dim_avg_n_Wrap(nEy(nb:ne,:,:), 0)
  copy_VarCoords(u(0,:,:), pEx1(:,:))
  copy_VarCoords(u(0,:,:), pEy1(:,:))
  copy_VarCoords(u(0,:,:), nEx1(:,:))
  copy_VarCoords(u(0,:,:), nEy1(:,:))
  pdiv1 = uv2dv_cfd(pEx1, pEy1, u&lat, u&lon, 0)
  nEx1 = pEx1 - nEx1
  nEy1 = pEy1 - nEy1
  ndiv1 = uv2dv_cfd(nEx1, nEy1, u&lat, u&lon, 0)*10^(5)

  pE1 = sqrt(pEx1*pEx1+pEy1*pEy1)
  nE1 = sqrt(nEx1*nEx1+nEy1*nEy1)

  pu = dim_avg_n_Wrap(put(nb:ne,:,:), 0)
  nu = dim_avg_n_Wrap(nut(nb:ne,:,:), 0)
  pp = pu - nu
  pp = smth9_Wrap(pp, 0.5, 0.25, False)
  copy_VarCoords(pu,pp)

  ;pEx1 = mask(pEx1, pE1.gt.10., True)
  ;pEy1 = mask(pEx1, pE1.gt.10., True)
  ;nEx1 = mask(pEx1, nE1.gt.10., True)
  ;nEy1 = mask(pEx1, nE1.gt.10., True)
  copy_VarCoords(u(0,:,:),pdiv1(:,:))
  copy_VarCoords(u(0,:,:),ndiv1(:,:))

  res@gsnLeftString = label(i)
  res@gsnRightString = "~F21~High-Low"
  res@gsnCenterString = "~F22~E~F21~&div"
  plotu(i) = gsn_csm_vector(wks, nEx1(::2,::2), nEy1(::2,::2), resw)
  plot(i)  = gsn_csm_contour_map(wks,ndiv1(:,:),res)
  plote(i) = gsn_csm_contour(wks, pp, resc)
  overlay(plot(i), plotu(i))
  overlay(plot(i), plote(i))
  delete([/pEx1,pEy1,nEx1,nEy1,pE1,nE1,ndiv1,pdiv1/])

  i = 2
;; day 2 4
  nb = 6
  ne = 8
  pEx1 = dim_avg_n_Wrap(pEx(nb:ne,:,:), 0)
  nEx1 = dim_avg_n_Wrap(nEx(nb:ne,:,:), 0)
  pEy1 = dim_avg_n_Wrap(pEy(nb:ne,:,:), 0)
  nEy1 = dim_avg_n_Wrap(nEy(nb:ne,:,:), 0)
  copy_VarCoords(u(0,:,:), pEx1(:,:))
  copy_VarCoords(u(0,:,:), pEy1(:,:))
  copy_VarCoords(u(0,:,:), nEx1(:,:))
  copy_VarCoords(u(0,:,:), nEy1(:,:))
  pdiv1 = uv2dv_cfd(pEx1, pEy1, u&lat, u&lon, 0)
  nEx1 = pEx1 - nEx1
  nEy1 = pEy1 - nEy1
  ndiv1 = uv2dv_cfd(nEx1, nEy1, u&lat, u&lon, 0)*10^(5)

  pE1 = sqrt(pEx1*pEx1+pEy1*pEy1)
  nE1 = sqrt(nEx1*nEx1+nEy1*nEy1)

  pu = dim_avg_n_Wrap(put(nb:ne,:,:), 0)
  nu = dim_avg_n_Wrap(nut(nb:ne,:,:), 0)
  pp = pu - nu
  pp = smth9_Wrap(pp, 0.5, 0.25, False)
  copy_VarCoords(pu,pp)

  ;pEx1 = mask(pEx1, pE1.gt.10., True)
  ;pEy1 = mask(pEx1, pE1.gt.10., True)
  ;nEx1 = mask(pEx1, nE1.gt.10., True)
  ;nEy1 = mask(pEx1, nE1.gt.10., True)
  copy_VarCoords(u(0,:,:),pdiv1(:,:))
  copy_VarCoords(u(0,:,:),ndiv1(:,:))

  resw@vcRefAnnoOn             = True
  resw@vcRefAnnoString1On      = True
  resw@vcRefAnnoString1        = "$VMG$" +" ~F21~m~S~2~N~ s~S~-2~N~"
  resw@vcRefAnnoOrthogonalPosF = -0.1
  resw@vcRefAnnoParallelPosF   = 1.2
  resw@vcRefAnnoFontHeightF    = 0.015

  res@gsnLeftString = label(i)
  res@gsnRightString = "~F21~High-Low"
  res@gsnCenterString = "~F22~E~F21~&div"
  plotu(i) = gsn_csm_vector(wks, nEx1(::2,::2), nEy1(::2,::2), resw)
  plot(i)  = gsn_csm_contour_map(wks,ndiv1(:,:),res)
  plote(i) = gsn_csm_contour(wks, pp, resc)
  overlay(plot(i), plotu(i))
  overlay(plot(i), plote(i))
  delete([/pEx1,pEy1,nEx1,nEy1,pE1,nE1,ndiv1,pdiv1/])

;; region
  lonw = -30
  lone =  10
  lats =  35
  latn =  75
  polyx = (/lonw,lone,lone,lonw,lonw/)
  polyy = (/lats,lats,latn,latn,lats/)
  respoly = True
  respoly@gsLineDashPattern = 0
  respoly@gsLineThicknessF  = 4.0
  ;dum3 = gsn_add_polyline(wks, plot(0), polyx, polyy, respoly)
  ;dum31 = gsn_add_polyline(wks, plot(0), (/-30,10/), (/55,55/), respoly)

  ;dum4 = gsn_add_polyline(wks, plot(1), polyx, polyy, respoly)
  ;dum41 = gsn_add_polyline(wks, plot(1), (/-30,10/), (/55,55/), respoly)

  ;dum5 = gsn_add_polyline(wks, plot(2), polyx, polyy, respoly)
  ;dum51 = gsn_add_polyline(wks, plot(2), (/-30,10/), (/55,55/), respoly)

  gsn_panel(wks,plot,(/3,1/),resP)    

  ;i = 0
  ;do while (i.eq.0)
  ;  i = 0
  ;end do