

latn = 80
lats = 10
lonw = 100
lone = 240

ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"
;----------------------------------------------------------------
;main code
;----------------------------------------------------------------
  f = addfile(ddir+"/ERA5_daily_anom_u200_1982-2022.nc","r")
  fu= addfile(ddir+"/ERA5_daily_mean_u200_1979-2024.nc","r")

  TIME=cd_calendar(f->time,0)
  timi=cd_calendar(f->time,-2)
  mon = TIME(:,1)
  nt = ind((timi.ge.19821201.and.timi.le.20220331)\
           .and.(mon.eq.12.or.mon.eq.1.or.mon.eq.2.or.mon.eq.3) )
  print(dimsizes(nt))

  x = f->u(nt,{lats:latn},{lonw:lone})
 
  um = fu->u(nt,{lats:latn},{lonw:lone})
  um := dim_avg_n_Wrap(um, 0)

  delete([/TIME,timi,mon,nt/])

  ;pcs
  DIR = ddir+"/pcs_u200.nc";"/NPJ_pcs.nc"
  fout = addfile(DIR,"r")
  pc  = fout->x

  r1 = regCoef_n( pc(0,:),x, 0, 0)
  r2 = regCoef_n( pc(1,:),x, 0, 0)
  copy_VarCoords(x(0,:,:), r1)
  copy_VarCoords(x(0,:,:), r2)

  wks  = gsn_open_wks("x11","regress")         ; send graphics to PNG file
   cmap = read_colormap_file("BlueWhiteOrangeRed")
  res                       = True     
  res@gsnMaximize           = True             ; make large
  res@gsnDraw  = False
  res@gsnAddCyclic  = False
  res@gsnFrame = False
  res@gsnMajorLatSpacing = 15

  res@cnFillOn              = True             ; turn on color
  res@cnLinesOn             = False            ; turn off contour lines
  res@cnLineLabelsOn        = False            ; turn off contour line labels

  res@mpFillOn              = False            ; turn off default background gray
  res@mpMinLatF = lats
  res@mpMaxLatF = latn
  res@mpMinLonF = lonw
  res@mpMaxLonF = lone

  res@cnFillOn               = True
  res@cnFillPalette          = cmap;(::-1,:)
  res@cnLinesOn              = False
  res@cnLineLabelsOn         = False
  res@cnLevelSelectionMode   = "Manuallevels"
  res@cnMinLevelValF         = -12
  res@cnMaxLevelValF         =  12
  res@cnLevelSpacingF        =  2
  ;res@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  res@lbLabelBarOn           = False

  resc = True
  resc@gsnFrame     = False
  resc@gsnDraw      = False
  resc@gsnMaximize  = False
  resc@gsnAddCyclic = False
  resc@gsnContourNegLineDashPattern = 1
  resc@gsnLeftString      = ""
  resc@gsnRightString     = ""
  resc@cnMonoLineColor    = True
  resc@cnLinesOn          = True
  resc@cnLineLabelsOn     = False
  resc@cnLineColor        = "black";magenta
  resc@cnLineThicknessF   = 3.0
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels  = (/30.,40.,50.,60./)  

  res@gsnLeftString ="(a) NPJ EOF1 (14.56%~F18~S~F21~0.25%)"
  plot1 = gsn_csm_contour_map_ce(wks,r1,res)

  res@gsnLeftString ="(b) NPJ EOF2 (10.64%~F18~S~F21~0.19%)"
  plot2 = gsn_csm_contour_map_ce(wks,r2,res)

  plotu1 = gsn_csm_contour(wks,um,resc)
  plotu2 = gsn_csm_contour(wks,um,resc)
  overlay(plot1, plotu1)
  overlay(plot2, plotu2)

  resP = True
  resP@gsnMaximize      = True
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02

  resP@pmLabelBarOrthogonalPosF = 0.0
  resP@pmLabelBarWidthF         = 0.7

  resP@lbTitleString      = "u200 [m s~S~-1~N~]"
  resP@lbTitlePosition    = "Bottom"
  resP@lbTitleDirection   = "Across"
  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012

  resP@gsnPanelFigureStringsPerimOn = False
  gsn_panel(wks, (/plot1,plot2/),(/1,1/), resP)


  i=0
  do while(i.eq.0)
    i=0
  end do