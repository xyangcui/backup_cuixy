  
  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  idir = ddir+"/ERA5_daily_mean_mhf850_1979-2022_r2.5.nc"

  f = addfile(idir,"r")
  x = f->mhf

  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)

  mon := mon(n1:n2)
  x  := x(n1:n2,:,:)

  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  x := x(nt,:,:)

  x := dim_avg_n_Wrap(x, 0)

  res = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False
  res@gsnAddCyclic = False
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
  res@mpCenterLonF    = -160
  res@mpMinLatF       = 15
  res@mpMaxLatF       = 80
  res@mpMinLonF       = -120
  res@mpMaxLonF       = 0

  res@cnFillOn = True
  res@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
  res@cnFillPalette        = "MPL_Blues";
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "Manuallevels"
  res@cnMinLevelValF       =  0
  res@cnMaxLevelValF       =  20
  res@cnLevelSpacingF      =  2

  res@pmLabelBarOrthogonalPosF = 0.1
  res@lbTitleString = "[K m s~S~-1~N~]"
  res@lbTitlePosition  = "Right" 
  res@lbTitleDirection = "Across" 
  res@lbTitleFontHeightF = 0.015

  wks = gsn_open_wks("x11", "eof")

  res@vpWidthF  = 0.45
  res@vpHeightF = 0.3
  res@vpXF = 0.06
  res@vpYF = 0.9
  res@gsnLeftString = "~F21~(a) DJF Storm Track"
  plot = gsn_csm_contour_map(wks, x, res)
  lonw = -70
  lone = -40
  lats =  35
  latn =  50
  polyx = (/lonw,lone,lone,lonw,lonw/)
  polyy = (/lats,lats,latn,latn,lats/)
  respoly = True
  respoly@gsLineDashPattern = 0
  respoly@gsLineThicknessF  = 4.0
  dum3 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)
  draw(plot)

;; draw lead-lag correlation.
  jdir = ddir+"/Ulevel_PC.csv"
  mdir = ddir+"/meri_heatflux.csv"
; read in jet index. DJF.
  delim  = ","
  jf  = asciiread(jdir,-1,"string")
  ji  = str_get_field(jf,1,delim)
  ji := tofloat(ji)

  mf  = asciiread(mdir,-1,"string")
  mi  = str_get_field(mf,1,delim)
  mi := tofloat(mi)
; time table.
  idir = ddir + "/ERA5_daily_anom_ulevel_1979-2022_r2.5.nc"
  f = addfile(idir,"r")
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  undef("mon")
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

;;calc self corr.
  mxlag    = 2
  xym = esccr(mi,mi,mxlag)
  yxm = esccr(mi,mi,mxlag)    
  ccrm = new ( 2*mxlag+1, float)    
  ccrm(0:mxlag-1) = xym(1:mxlag:-1)  
  ccrm(mxlag:)    = xym(0:mxlag)    
  vm = -1*dimsizes(mi)*0.5*log(ccrm(1))
  xyj = esccr(ji,ji,mxlag)
  yxj = esccr(ji,ji,mxlag)    
  ccrj = new ( 2*mxlag+1, float)    
  ccrj(0:mxlag-1) = xyj(1:mxlag:-1)  
  ccrj(mxlag:)    = xyj(0:mxlag)    
  vj = -1*dimsizes(mi)*0.5*log(ccrj(1))
  vv = (vj+vm)/2
  print(vv)
;; calc corr.
  mxlag    = 12
  x_Lead_y = esccr(mi,ji,mxlag)
  y_Lead_x = esccr(ji,mi,mxlag)  

  ccr = new ( 2*mxlag+1, float)  
  ccr(0:mxlag-1) = x_Lead_y(1:mxlag:-1)  
  ccr(mxlag:)    = y_Lead_x(0:mxlag)   
  ;ccr(0:mxlag-1) = x_Lead_y(1:mxlag:-1)  
  ;ccr(mxlag:)    = y_Lead_x(0:mxlag)    

  undef("nt")
  nt = ind(ccr.gt.0.15.or.ccr.lt.-0.15)
  value = fspan(0., 2.4, 25)
;;0.154
  resx = True
  resx@gsnFrame = False
  resx@gsnDraw  = False
  resx@gsnMaximize = False
  resx@gsnLeftStringOrthogonalPosF  = 0.012
  resx@gsnRightStringOrthogonalPosF = 0.012
  resx@gsnLeftStringFontHeightF  = 0.016
  resx@gsnRightStringFontHeightF = 0.016
  resx@tmXBLabelFontHeightF  = 0.012
  resx@tmYLLabelFontHeightF  = 0.012
  resx@tiXAxisString = "lead                           lag"
  resx@xyDashPattern  = 0
  resx@gsnYRefLine = 0
  resx@xyLineThicknessF = 3                          
  resx@tiMainString   = ""   
  resx@tmXBMode   = "Explicit"
  resx@tmXBValues = (/ 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,\
                    1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4 /)
  resx@tmXBLabels = (/"12","","","","8","","","","4","","","",\
                     "0","","","","4","","","","8",\
                     "","","","12"/)  
  resx@trXMaxF = 2.4
  resx@vpWidthF  = 0.45
  resx@vpHeightF = 0.25
  resx@vpXF = 0.06
  resx@vpYF = 0.5
  resx@gsnLeftString = "~F21~(b) Lead-lag correlation   MHF lead Jet"
  plotxy  = gsn_csm_xy (wks,fspan(0., 2.4, 25),ccr,resx) 
;;dot
  resp = True
  resp@gsMarkerColor = "red"
  resp@gsMarkerIndex = 16
  resp@gsMarkerThicknessF = 0.1
  dum = new(dimsizes(nt),"graphic")
  do i = 0,dimsizes(dum)-1
    ni = nt(i)
    dum(i) = gsn_add_polymarker(wks, plotxy, value(ni), ccr(ni), resp)
  end do
  print("done.")
  draw(plotxy)

;; draw com event.
  lats =  10
  latn =  80
  lonw = -60
  lone =  30

  undef("indir")
  indir = ddir+"/mhf_events.csv"

;; EOF pattern
  edir = ddir+"/EOF_ulevel.nc"
  fe = addfile(edir,"r")
  eof1 = fe->u(0,{250},:,:)
; time table.
  idir = ddir + "/ERA5_daily_anom_ulevel_1979-2022_r2.5.nc"
  f = addfile(idir,"r")
  undef("u")
  u = f->u(:,{250},{lats:latn},{lonw:lone})
  printVarSummary(u)
  time = f->time
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

  pu = new((/dimsizes(pev),5,dimsizes(u&lat),dimsizes(u&lon)/),float)
  pu = 0.
  undef("nt")
  do i = 0, dimsizes(pev)-1
    nt = ind(year.eq.pyear(i).and.mon.eq.pmon(i).and.day.eq.pday(i))
    pu(i,:,:,:) = u(nt+1:nt+5,:,:)
    delete(nt)
  end do

  nev   = ind(pct.lt.0)
  nyear = iyear(nev)
  nmon  = imon(nev)
  nday  = iday(nev)

  nu = new((/dimsizes(nev),5,dimsizes(u&lat),dimsizes(u&lon)/),float)
  nu = 0.

  do i = 0, dimsizes(nev)-1
    nt = ind(year.eq.nyear(i).and.mon.eq.nmon(i).and.day.eq.nday(i))
    nu(i,:,:,:) = u(nt+1:nt+5,:,:)
    delete(nt)
  end do

  copy_VarCoords(u(0,:,:), pu(0,0,:,:))
  copy_VarCoords(u(0,:,:), nu(0,0,:,:))

  pu := dim_avg_n_Wrap(pu, 0)
  nu := dim_avg_n_Wrap(nu, 0)

  pu := dim_avg_n_Wrap(pu, 0)
  nu := dim_avg_n_Wrap(nu, 0)

  clat = cos(0.01745329*u&lat)
  pcc1 = pattern_cor(pu, eof1,clat, 0)
  pcc2 = pattern_cor(nu, eof1,clat, 0)
;;plot
  undef("res")         
  res              = True         
  res@gsnDraw      = False        
  res@gsnFrame     = False       
  res@gsnAddCyclic = False       
  res@gsnLeftStringOrthogonalPosF  = 0.012
  res@gsnRightStringOrthogonalPosF = 0.012
  res@gsnLeftStringFontHeightF  = 0.016
  res@gsnRightStringFontHeightF = 0.016
  res@tmXBLabelFontHeightF  = 0.012
  res@tmYLLabelFontHeightF  = 0.012
  res@gsnMajorLatSpacing = 15

  res@mpFillOn     = False        
  res@mpMinLatF    = lats        
  res@mpMaxLatF    = latn
  res@mpMinLonF    = lonw
  res@mpMaxLonF    = lone
                                        
  cmap = read_colormap_file("MPL_RdBu")
  res@cnFillOn = True
  res@cnMonoLineColor    = True
  res@cnLinesOn          = False
  res@cnLineLabelsOn     = False
  res@cnLineColor        = "black"
  res@cnLineThicknessF   = 2.0
  res@cnInfoLabelPerimOn = False
  res@cnInfoLabelOn      = False
  res@cnFillPalette = cmap(::-1,:)
  res@cnLevelSelectionMode   = "Manuallevels"

  ;res@cnMinLevelValF         = -1.
  ;res@cnMaxLevelValF         =  1.
  ;res@cnLevelSpacingF        =  0.1 
  symMinMaxPlt(pu, 11, False, res)

  resP                    = True         
  resP@gsnMaximize        = True         
  resP@gsnPanelLabelBar   = True        


  res@lbLabelBarOn = False
  res@vpXF = 0.52
  res@vpYF = 0.87
  res@vpWidthF  = 0.45
  res@vpHeightF = 0.25
  res@gsnLeftString  = "~F21~(c) u250_High MHF"
  res@gsnRightString = "";PCC: "+sprintf("%5.2f", pcc1)
  plot1 =gsn_csm_contour_map(wks,pu,res)

  res@lbLabelBarOn = True
  res@pmLabelBarOrthogonalPosF = 0.15
  res@lbTitleString = "[m s~S~-1~N~]"
  res@lbTitlePosition  = "Bottom" 
  res@lbTitleDirection = "Across" 
  res@lbTitleFontHeightF = 0.015

  res@vpYF = 0.53
  res@gsnLeftString = "~F21~(d) u250_Low MHF"
  res@gsnRightString = "";PCC: "+sprintf("%5.2f", pcc2)
  plot2 =gsn_csm_contour_map(wks,nu,res)

  draw(plot1)    
  draw(plot2) 

  i = 0
  do while (i.eq.0)
    i = 0
  end do