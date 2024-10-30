
  lats = 20
  latn = 70


  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  zdir = ddir+"/"+"ERA5_z500_daily_anom_1979-2022_r2.5.nc"
  idir = ddir+"/"+"pattern_consistency_omi_alldays.nc"
  sdir = ddir+"/"+"pattern_consistency_omi_all_alldays.nc"

  fz = addfile(zdir,"r")
  z  = fz->z(0,10:18,0)

  fi = addfile(idir,"r")
  IO = fi->pc
  IO = IO * 100

  fs = addfile(sdir,"r")
  IOa = fs->pc
  IOa = IOa * 100
  print(dimsizes(IOa))

  IOnP = new((/9,31/),float)
  IOnP = 0.
  IOiP = IOnP
  IOfP = IOnP
  IOsP = IOnP

  IOnP(1:,:) = IOa
  IOnP(0,:)  = IOa(7,:)
  copy_VarCoords(z, IOnP(:,0))
  IOnP!1 = "time"
  IOnP := IOnP(time|:,lat|:)

  IOiP(1:,:) = IO(1,:,:)
  IOiP(0,:)  = IO(1,7,:)
  copy_VarCoords(z, IOiP(:,0))
  IOiP!1 = "time"
  IOiP := IOiP(time|:,lat|:)

  IOfP(1:,:) = IO(2,:,:)
  IOfP(0,:)  = IO(2,7,:)
  copy_VarCoords(z, IOfP(:,0))
  IOfP!1 = "time"
  IOfP := IOfP(time|:,lat|:)

  IOsP(1:,:) = IO(3,:,:)
  IOsP(0,:)  = IO(3,7,:)
  copy_VarCoords(z, IOsP(:,0))
  IOsP!1 = "time"
  IOsP := IOsP(time|:,lat|:)
;;diff
  undef("fi")
  fi = IOfP - IOsP
  si = IOsP - IOiP
  ni = si
  ni(:,1:5) = IOnP(:,1:5) - IOiP(:,1:5)
  ni(:,0)   = 0
  ni(:,6:8) = 0
  copy_VarCoords(IOnP,ni)
  copy_VarCoords(IOnP,fi)
  copy_VarCoords(IOnP,si)

  wks  = gsn_open_wks ("x11", "hov" )  

  cmap  = read_colormap_file("MPL_hot")
  cmap2 = read_colormap_file("BlueWhiteOrangeRed")        

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
  res@tmXBLabels = (/"8","1","2","3","4","5","6","7","8" /)

  res@cnFillOn       = True              
  ;res@cnFillPalette  = cmap(::-1,:)  
  res@cnLinesOn = False    
  res@cnLevelSelectionMode = "ManualLevels"     
  res@cnMinLevelValF       =  10;10         
  res@cnMaxLevelValF       =  15;25              
  res@cnLevelSpacingF      =  0.2;1

  res@lbOrientation  = "Vertical" 
  res@lbLabelStride  = 2
     
  res@vpWidthF  = 0.25
  res@vpHeightF = 0.35

  label = (/"(a) IOnP","(b) IOiP","(c) IOfP","(d) IOsP"/)
            
  res@cnFillPalette  = cmap(::-1,:)   
  res@cnLevelSelectionMode = "ManualLevels"     
  res@cnMinLevelValF       =  10         
  res@cnMaxLevelValF       =  15              
  res@cnLevelSpacingF      =  0.2

  res@vpXF = 0.18
  res@vpYF = 0.95
  res@gsnLeftString = "~F21~"+label(0)
  plot1 = gsn_csm_time_lat(wks, IOnP, res) 
  ;plot1 = gsn_csm_time_lat(wks, fi, res) 

  res@cnFillOn       = True              
  res@cnFillPalette  = cmap(::-1,:)  
  res@cnLevelSelectionMode = "ManualLevels"     
  res@cnMinLevelValF       =  10;10         
  res@cnMaxLevelValF       =  25;25              
  res@cnLevelSpacingF      =  1;1

  res@vpXF = 0.6
  res@vpYF = 0.95
  res@gsnLeftString = "~F21~"+label(1)
  plot2 = gsn_csm_time_lat(wks, IOiP, res) 

  res@vpXF = 0.18
  res@vpYF = 0.47
  res@gsnLeftString = "~F21~"+label(2)
  plot3 = gsn_csm_time_lat(wks, IOfP, res) 
  ;plot3 = gsn_csm_time_lat(wks, fi, res) 

  res@vpXF = 0.6
  res@vpYF = 0.47
  res@gsnLeftString = "~F21~"+label(3)
  plot4 = gsn_csm_time_lat(wks, IOsP, res) 
  ;plot4 = gsn_csm_time_lat(wks, si, res) 
  draw(plot1)
  draw(plot2)
  draw(plot3)
  draw(plot4)

  i = 0
  do while(i.eq.0)
    i = 0
  end do
