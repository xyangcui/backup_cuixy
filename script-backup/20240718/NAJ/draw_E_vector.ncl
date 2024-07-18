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
  pdir = ddir+"/Ulevel_PC"
; time table.
  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  fu   = addfile(udir,"r")
  fv   = addfile(vdir,"r")

  u = fu->u(:,{lats:latn},{lonw:lone})
  v = fv->v(:,{lats:latn},{lonw:lone})

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

  do i = 0, dimsizes(pev)-1
    nt = ind(year.eq.pyear(i).and.mon.eq.pmon(i).and.day.eq.pday(i))
    pu(i,:,:,:) = u(nt-4:nt+6,:,:)
    pv(i,:,:,:) = v(nt-4:nt+6,:,:)
    delete(nt)
  end do

  nev   = ind(pct.lt.0)
  nyear = iyear(nev)
  nmon  = imon(nev)
  nday  = iday(nev)

  nu = new((/dimsizes(nev),11,dimsizes(u&lat),dimsizes(u&lon)/),float)
  nu = 0.
  nv = nu

  do i = 0, dimsizes(nev)-1
    nt = ind(year.eq.nyear(i).and.mon.eq.nmon(i).and.day.eq.nday(i))
    nu(i,:,:,:) = u(nt-4:nt+6,:,:)
    nv(i,:,:,:) = v(nt-4:nt+6,:,:)
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

  nb = 0
  ne = 2
  pEx := dim_avg_n_Wrap(pEx(nb:ne,:,:), 0)
  nEx := dim_avg_n_Wrap(nEx(nb:ne,:,:), 0)
  pEy := dim_avg_n_Wrap(pEy(nb:ne,:,:), 0)
  nEy := dim_avg_n_Wrap(nEy(nb:ne,:,:), 0)

  copy_VarCoords(u(0,:,:), pEx(:,:))
  copy_VarCoords(u(0,:,:), pEy(:,:))

  copy_VarCoords(u(0,:,:), nEx(:,:))
  copy_VarCoords(u(0,:,:), nEy(:,:))

  pdiv = uv2dv_cfd(pEx, pEy, u&lat, u&lon, 0)

  nEx = pEx - nEx
  nEy = pEy - nEy
  ndiv = uv2dv_cfd(nEx, nEy, u&lat, u&lon, 0)

  pE = sqrt(pEx*pEx+pEy*pEy)
  nE = sqrt(nEx*nEx+nEy*nEy)

  pEx = mask(pEx, pE.gt.10., True)
  pEy = mask(pEx, pE.gt.10., True)
  nEx = mask(pEx, nE.gt.10., True)
  nEy = mask(pEx, nE.gt.10., True)
  copy_VarCoords(u(0,:,:),pdiv(:,:))
  copy_VarCoords(u(0,:,:),ndiv(:,:))
;;plot
  wks = gsn_open_wks("x11","eof")         
  plot = new(2,graphic)                
  plotc= plot         
  res              = True         
  res@gsnDraw      = False        
  res@gsnFrame     = False       
  res@gsnAddCyclic = False       
 
  res@mpFillOn     = False        
  res@mpMinLatF    = lats        
  res@mpMaxLatF    = latn
  res@mpMinLonF    = lonw
  res@mpMaxLonF    = lone
  res@cnFillOn      = True                
  res@cnFillPalette = "BlWhRe" 
  ;res@cnLevelSelectionMode   = "Manuallevels"
  ;res@cnMinLevelValF         = -1.
  ;res@cnMaxLevelValF         =  1.
  ;res@cnLevelSpacingF        =  0.1                                       
  symMinMaxPlt(pdiv, 11, False, res)   

  ;res = True
  ;res@gsnFrame     = False
  ;res@gsnDraw      = False
  ;res@gsnMaximize  = False
  ;res@gsnAddCyclic = False
  ;res@gsnContourNegLineDashPattern = 1
  ;res@gsnLeftString      = ""
  ;res@gsnRightString     = ""
  ;res@cnMonoLineColor    = True
  ;res@cnLinesOn          = True
  ;res@cnLineLabelsOn     = False
  ;res@cnLineColor        = "black"
  ;res@cnLineThicknessF   = 2.0
  ;res@cnInfoLabelPerimOn = False
  ;res@cnInfoLabelOn      = False
  ;res@cnLevelSelectionMode = "Explicitlevels"
  ;res@cnLevels   = (/20.,25.,30.,35./)

  resw = True
  resw@gsnFrame     = False
  resw@gsnDraw      = False
  resw@gsnMaximize  = False
  resw@gsnAddCyclic = False
  resw@vcGlyphStyle            = "LineArrow"
  resw@vcLineArrowThicknessF   = 2.5
  ;resw@vcRefMagnitudeF         = 4.
  resw@vcRefLengthF            = 0.05
  resw@vcLineArrowColor        = "mediumblue"
  resw@vcRefAnnoString2On      = False
  resw@vcRefAnnoString1        = "~F21~$VMG$" +" m s~S~-1~N~"
  resw@vcRefAnnoPerimOn        = False
  resw@vcRefAnnoBackgroundColor = -1

  resP                    = True         
  resP@gsnMaximize        = True         
  resP@gsnPanelLabelBar   = True        
;*******************************************
; first plot
;*******************************************
  ;ndiv = pdiv-ndiv
  ;copy_VarCoords(pdiv,ndiv)
  plotu  = gsn_csm_vector(wks, pEx(::2,::2), pEy(::2,::2), resw)
  plotu2 = gsn_csm_vector(wks, nEx(::2,::2), nEy(::2,::2), resw)
  plot(0)= gsn_csm_contour_map(wks,pdiv(:,:),res)
  plot(1)= gsn_csm_contour_map(wks,ndiv(:,:),res)
  overlay(plot(0), plotu)
  overlay(plot(1), plotu2)

  lonw = -30
  lone =  10
  lats =  35
  latn =  75
  polyx = (/lonw,lone,lone,lonw,lonw/)
  polyy = (/lats,lats,latn,latn,lats/)
  respoly = True
  respoly@gsLineDashPattern = 0
  respoly@gsLineThicknessF  = 4.0
  dum3 = gsn_add_polyline(wks, plot(1), polyx, polyy, respoly)
  dum31 = gsn_add_polyline(wks, plot(1), (/-30,10/), (/55,55/), respoly)
  dum4 = gsn_add_polyline(wks, plot(0), polyx, polyy, respoly)
  dum41 = gsn_add_polyline(wks, plot(0), (/-30,10/), (/55,55/), respoly)

  gsn_panel(wks,plot,(/2,1/),resP)    

  i = 0
  do while (i.eq.0)
    i = 0
  end do