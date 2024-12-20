

model_name = "ecmwf"

ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/hindcast/"
sdir = ddir+model_name+"/"+"acc_IOsP"+"_p"+i+"_a1.nc"
fdir = ddir+model_name+"/"+"acc_IOfP"+"_p"+i+"_a1.nc"
dfdir = ddir+model_name+"/"+"acc_diff_"+"p"+i+"_a1.nc" 

lats = 10
latn = 80
lonw = 140
lone = 270


  ff = addfile(ddir+model_name+"/"+"acc_IOfP"+"_p"+i+"_a1.nc" ,"r")  
  fd = addfile(ddir+model_name+"/"+"acc_diff"+"_p"+i+"_a1.nc" ,"r")  
  fs = addfile(ddir+model_name+"/"+"acc_IOsP"+"_p"+i+"_a1.nc" ,"r")  

  acc_s = fs->acc
  acc_f = ff->acc
  acc_d = fd->acc

  wks = gsn_open_wks("x11","acc")
  cmap = read_colormap_file("MPL_RdBu") 
  cmap2 = read_colormap_file("GMT_gray")

  res = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False
  res@gsnAddCyclic = False
  res@gsnMajorLatSpacing  = 20
  res@gsnMajorLonSpacing  = 30
  res@gsnLeftString = ""
  res@gsnRightString = ""

  res@gsnLeftStringOrthogonalPosF  = 0.012
  res@gsnRightStringOrthogonalPosF = 0.012
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringFontHeightF = 0.015
  res@tmXBLabelFontHeightF  = 0.012
  res@tmYLLabelFontHeightF  = 0.012

  res@mpFillOn = True
  res@mpLandFillColor = "white"
  res@mpMinLatF       = lats
  res@mpMaxLatF       = latn
  res@mpMinLonF       = lonw
  res@mpMaxLonF       = lone

  res@cnFillOn       = True              
 
  res@cnLinesOn = False   
  res@cnLevelSelectionMode = "ManualLevels"     
  res@cnInfoLabelPerimOn = False
  res@cnInfoLabelOn      = False
  res@cnLineLabelsOn     = False

  res@lbLabelBarOn = False
  res@pmLabelBarOrthogonalPosF = 0.23
  res@lbTitleOn = True
  res@lbTitleString = "[m s~S~-1~N~]"
  res@lbTitlePosition  = "Right"
  res@lbTitleDirection = "Across"
  res@lbTitleFontHeightF = 0.015
 
  res@tmXBLabelFont = 21
  res@tmYRLabelFont = 21 

; for diff
  resc = res 
  resc@cnFillPalette  = dcmap(::-1,:) 
  resc@cnMinLevelValF       = -1.0 
  resc@cnMaxLevelValF       = 1.0        
  resc@cnLevelSpacingF      = 0.1

  res@cnFillPalette  = dcmap(::-1,:) 
  res@cnMinLevelValF       = -1.0 
  res@cnMaxLevelValF       =  1.0        
  res@cnLevelSpacingF      =  0.1

  label = (/"~F21~(a) IOsP","~F21~(b) IOfP","~F21~(c) IOsP-IOfP"/)
  plots = new(4,graphic)
  plotf = new(4,graphic)
  plotd = new(4,graphic)

  do j = 0, 3

    plots(j) = gsn_csm_contour_map(wks, acc_s(j,:,:), res)
    plotf(j) = gsn_csm_contour_map(wks, acc_f(j,:,:), res)
    plotd(j) = gsn_csm_contour_map(wks, acc_d(j,:,:), resc)

  end do

  resP = True    
  resP@gsnMaximize        = True         
  resP@gsnPanelLabelBar   = False
  resP@gsnPanelRowSpec    = True

  gsn_panel(wks,(/plots(0),plotf(0),plotd(0),\
                  plots(1),plotf(1),plotd(1),\
                  plots(2),plotf(2),plotd(2),\
                  plots(3),plotf(3),plotd(3)/),(/3,3,3,3/),resP)  

  i = 0
  do while i.eq.0
    i = 0
  end do