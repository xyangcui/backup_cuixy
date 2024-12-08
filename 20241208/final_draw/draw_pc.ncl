
  lats = 20
  latn = 70


  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  zdir = ddir+"/"+"ERA5_z500_daily_anom_1979-2022_r2.5.nc"
  ;idir = ddir+"/"+"pattern_consistency_omi_alldays_2phases.nc"
  idir = ddir+"/"+"pattern_consistency_omi_alldays.nc"
  sdir = ddir+"/"+"pattern_consistency_omi_all_alldays.nc"

  fz = addfile(zdir,"r")
  z  = fz->z(0,10:18,0)
  lat = z&lat

  fi = addfile(idir,"r")
  IO = fi->pc
  IO = IO * 100
  printVarSummary(IO)

  IOnP = new((/9,26/),float)
  IOnP = 0.
  IOiP = IOnP
  IOfP = IOnP
  IOsP = IOnP

  IOfP(0:7,:) = (/IO(0,:,:25)/)
  IOfP(8,:)  = (/IO(0,0,:25)/)
  copy_VarCoords(z, IOfP(:,0))
  IOfP!1 = "time"
  IOfP := IOfP(time|:,lat|:)

  IOsP(0:7,:) = (/IO(1,:,:25)/)
  IOsP(8,:)  = (/IO(1,0,:25)/)
  copy_VarCoords(z, IOsP(:,0))
  IOsP!1 = "time"
  IOsP := IOsP(time|:,lat|:)

  undef("fi")
  fi = IOsP - IOfP
  copy_VarCoords(IOsP,fi)

;smooth
  fi = smth9_Wrap(fi, 0.5, -0.25, False)

  IOfP = smth9_Wrap(IOfP, 0.5, -0.25, False)
  IOsP = smth9_Wrap(IOsP, 0.5, -0.25, False)
;plot
  cmap  = read_colormap_file("MPL_hot")
  cmap2 = read_colormap_file("ViBlGrWhYeOrRe")        

  res             = True  
  res@gsnDraw     = False
  res@gsnMaximize = False
  res@gsnFrame    = False
  res@gsnRightStringFontHeightF = 0.015
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringOrthogonalPosF = 0.0
  res@gsnLeftStringOrthogonalPosF  = 0.0

  res@tmYLLabelFontHeightF = 0.012
  res@tmXBLabelFontHeightF = 0.012

  res@tmYROn = False
  res@tmXTOn = False

  res@tiXAxisString = "~F21~RMM phase"
  res@tiYAxisString = "~F21~Lag days"
  res@tiYAxisFontHeightF = 0.015
  res@tiXAxisFontHeightF = 0.015

  res@tmXBMode = "Explicit" 
  res@tmXBValues = z&lat
  res@tmXBLabels = (/"1","2","3","4","5","6","7","8","1" /)

  res@cnFillOn  = True              
  res@cnLinesOn = False    
  res@cnLineLabelsOn = False
  res@cnInfoLabelOn  = False

  res@lbOrientation  = "Vertical" 
  res@lbLabelStride  = 1
     
  res@vpWidthF  = 0.22
  res@vpHeightF = 0.35

  resd = res

  resp = True
  resp@gsLineColor = "black"
  resp@gsLineThicknessF = 3.
  resp@gsLineDashPattern = 0


  wks  = gsn_open_wks ("x11", "hov" )  

  llat = lat(4)
  polyx = (/llat-1,llat+1,llat+1,llat-1,llat-1/)
  polyy = (/15,15,21,21,15/)

  label = (/"(a) IOsP","(b) IOfP","(c) IOsP - IOfP"/)
            
  res@vpXF = 0.08
  res@vpYF = 0.8           
  res@cnFillPalette  = cmap(::-1,:)
  res@cnLevelSelectionMode = "ManualLevels"     
  res@cnMinLevelValF       = 13      
  res@cnMaxLevelValF       = 20        
  res@cnLevelSpacingF      = 0.5
  res@lbLabelStride  = 2
  res@lbLabelBarOn  = False
  res@gsnLeftString = "~F21~"+label(0)
  plotf = gsn_csm_time_lat(wks, IOsP, res) 
  ;dumf  = gsn_add_polyline(wks, plotf, polyx, polyy, resp)

  res@vpXF = 0.36
  res@vpYF = 0.8
  res@lbLabelBarOn  = True         
  res@gsnLeftString = "~F21~"+label(1)
  res@tiYAxisString = ""
  plots = gsn_csm_time_lat(wks, IOfP, res) 
  ;dums  = gsn_add_polyline(wks, plots, polyx, polyy, resp)


  resd@cnFillOn       = True              
  resd@cnFillPalette  = cmap2 
  resd@cnLevelSelectionMode = "ManualLevels"     
  resd@cnMinLevelValF       = -6        
  resd@cnMaxLevelValF       =  6              
  resd@cnLevelSpacingF      =  1

  resd@vpXF = 0.7
  resd@vpYF = 0.8
  resd@gsnLeftString = "~F21~"+label(2) 
  resd@tiYAxisString = ""
  plotd = gsn_csm_time_lat(wks, fi, resd) 
  ;dumd  = gsn_add_polyline(wks, plotd, polyx, polyy, resp)

  draw(plotd)
  draw(plots)
  draw(plotf)

  i = 0
  do while(i.eq.0)
    i = 0
  end do
