

function return_acc(dir,model)
local f,acc,step,nt
begin
  lats = 10
  latn = 80
  lonw = 100
  lone = 240
  f = addfile(dir,"r")
  acc = f->acc(:,{lats:latn},{lonw:lone})
  printVarSummary(acc)
  step= acc&lead_time
  if model.eq."LIM"
    step = step*24
  end if
  nt = new(4,integer)
  do i = 0,3
    nt(i) = ind(step.eq.24*7+24*7*i)
  end do
  return acc(nt,:,:)
end

ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"
lats = 10
latn = 80
lonw = 100
lone = 240


model = (/"cma","cnrm","eccc","ecmwf",\
          "hmcr","isac","kma","ncep","ukmo"/)

  ldir = ddir+"/"+"LIM_theskill"+".nc"
  f = addfile(ldir,"r")
  x = f->acc(:,{lats:latn},{lonw:lone})
  dims = dimsizes(x)

  ;acc_mod = return_acc(ddir+"/"+"ecmwf_theskill"+".nc","ec")
  acc_mod = return_acc(ddir+"/"+"ecmwf_SNRskill"+".nc","ec")
  copy_VarCoords(x(0,:,:), acc_mod(0,:,:))

  acc_LIM = return_acc(ddir+"/"+"ecmwf_theskill"+".nc","ec")
  copy_VarCoords(x(0,:,:), acc_LIM(0,:,:))
  ;acc_LIM = acc_LIM - acc_mod
  printVarSummary(acc_LIM)
  printVarSummary(acc_mod)
  ;draw        
  cmap = read_colormap_file("WhiteBlueGreenYellowRed")
  res                       = True     
  res@gsnMaximize           = False             
  res@gsnDraw               = False
  res@gsnAddCyclic          = False
  res@gsnFrame              = False
  res@gsnMajorLatSpacing    = 10
  res@gsnLeftStringFontHeightF = 0.03
  res@gsnLeftStringOrthogonalPosF = 0.01

  res@tmXTOn = False
  res@tmYROn = False

  res@cnFillOn              = True             
  res@cnLinesOn             = False            
  res@cnLineLabelsOn        = False            

  res@tmYLLabelFontHeightF  = 0.02
  res@tmXBLabelFontHeightF  = 0.02

  res@mpFillOn              = False 
  res@mpMinLatF = lats
  res@mpMaxLatF = latn
  res@mpMinLonF = lonw
  res@mpMaxLonF = lone
  res@mpCenterLonF=180

  res@cnFillOn               = True
  res@cnFillPalette          = "BlueWhiteOrangeRed";cmap;(::-1,:)
  res@cnLinesOn              = False
  res@cnLineLabelsOn         = False
  res@cnLevelSelectionMode   = "Manuallevels"
  res@cnMinLevelValF         = 0.
  res@cnMaxLevelValF         = 1.
  res@cnLevelSpacingF        = 0.05
  res@lbLabelBarOn           = False

  wks   = gsn_open_wks("x11","regress")
  plot1 = new(4,graphic)
  plot2 = new(4,graphic)

  do i = 0, 3
    res@gsnLeftString = "LIM_"+"day"+(i*2+7)
    plot1(i) = gsn_csm_contour_map(wks, acc_LIM(i,:,:), res)
    res@gsnLeftString = "IFS_"+"day"+(i*7+7)
    plot2(i) = gsn_csm_contour_map(wks, acc_mod(i,:,:), res)
  end do


  resP = True
  resP@gsnMaximize      = True
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02

  resP@pmLabelBarOrthogonalPosF = 0.0
  resP@pmLabelBarWidthF         = 0.65

  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012
  resP@gsnPanelFigureStringsPerimOn = False
  resP@gsnPanelYWhiteSpacePercent = 3

  gsn_panel(wks, (/plot1(0),plot2(0),\
                   plot1(1),plot2(1),\
                   plot1(2),plot2(2),\
                   plot1(3),plot2(3)/), (/2,2,2,2/), resP)

  i=0
  do while(i.eq.0)
    i=0
  end do