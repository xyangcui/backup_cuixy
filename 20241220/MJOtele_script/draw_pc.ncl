
  lats = 20
  latn = 70


  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  zdir = ddir+"/"+"ERA5_z500_daily_anom_1979-2022_r2.5.nc"
  ;idir = ddir+"/"+"pattern_consistency_omi_alldays_2phases.nc"
  idir = ddir+"/"+"pattern_consistency_omi_alldays.nc"
  sdir = ddir+"/"+"pattern_consistency_omi_all_alldays.nc"

  fz = addfile(zdir,"r")
  z  = fz->z(0,10:18,0)

  fi = addfile(idir,"r")
  IO = fi->pc
  IO = IO * 100


  IOnP = new((/9,31/),float)
  IOnP = 0.
  IOiP = IOnP
  IOfP = IOnP
  IOsP = IOnP

  IOfP(0:7,:) = IO(2,:,:)
  IOfP(8,:)  = IO(2,0,:)
  copy_VarCoords(z, IOfP(:,0))
  IOfP!1 = "time"
  IOfP := IOfP(time|:,lat|:)

  IOsP(0:7,:) = IO(3,:,:)
  IOsP(8,:)  = IO(3,0,:)
  copy_VarCoords(z, IOsP(:,0))
  IOsP!1 = "time"
  IOsP := IOsP(time|:,lat|:)
;;diff
  undef("fi")
  fi = IOsP - IOfP
  copy_VarCoords(IOsP,fi)

  wks  = gsn_open_wks ("x11", "hov" )  

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

  res@lbOrientation  = "Vertical" 
  res@lbLabelStride  = 2
     
  res@vpWidthF  = 0.25
  res@vpHeightF = 0.35

  label = (/"(a) IOnP","(b) IOiP","(c) IOfP","(d) IOsP"/)
            

  res@cnFillOn       = True              
  res@cnFillPalette  = cmap(::-1,:)  
  res@cnLevelSelectionMode = "ManualLevels"     
  res@cnMinLevelValF       =  10;-20;10         
  res@cnMaxLevelValF       =  25 ;20;25              
  res@cnLevelSpacingF      =  1; 2;1

  res@vpXF = 0.18
  res@vpYF = 0.47
  res@gsnLeftString = "~F21~"+label(2)
  ;plot3 = gsn_csm_time_lat(wks, IOfP, res) 
  plot3 = gsn_csm_time_lat(wks, IOfP, res) 

  res@vpXF = 0.6
  res@vpYF = 0.47
  res@gsnLeftString = "~F21~"+label(3)
  plot4 = gsn_csm_time_lat(wks, IOsP, res) 
  ;plot4 = gsn_csm_time_lat(wks, si, res) 
  draw(plot3)
  draw(plot4)

  i = 0
  do while(i.eq.0)
    i = 0
  end do
