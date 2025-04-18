

latn = 80
lats = 10
lonw = 100
lone = 240

ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"
;----------------------------------------------------------------
;main code
;----------------------------------------------------------------

  fu= addfile(ddir+"/ERA5_daily_mean_u200_1979-2024.nc","r")
  TIME=cd_calendar(fu->time,0)
  timi=cd_calendar(fu->time,-2)
  mon = TIME(:,1)
  nt = ind((timi.ge.19821201.and.timi.le.20220331)\
           .and.(mon.eq.12.or.mon.eq.1.or.mon.eq.2.or.mon.eq.3) )
  um = fu->u(nt,{lats:latn},{lonw:lone})
  um := dim_avg_n_Wrap(um, 0)
  delete([/TIME,timi,mon,nt/])

  f = addfile(ddir+"/PrCA_patterns.nc","r")
  x  = f->predictable_patterns
  var= tofloat(x@ratio*100)
;;read in apt
  f = addfile(ddir+"/PrCA_apt.nc","r")
  apt = tofloat(f->apt)

  cmap = read_colormap_file("BlueWhiteOrangeRed")
  res                = True     
  res@gsnMaximize    = True             ; make large
  res@gsnDraw        = False
  res@gsnAddCyclic   = False
  res@gsnFrame       = False
  res@gsnMajorLatSpacing = 15
  res@gsnRightStringFontHeightF = 0.03
  res@gsnLeftStringFontHeightF  = 0.03
  res@tmYLLabelFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.02

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


  label = (/"(a) MPM1","(b) MPM2","(c) MPM3","(d) MPM4",\
            "(e) MPM5","(f) MPM6","(g) MPM7","(h) MPM8"/)

  wks  = gsn_open_wks("eps","/home/sunming/data5/cuixy/Subpre_NPJ/"+"figure/PrCA_patterns")         ; send graphics to PNG file
  plot = new(8,graphic)
  plotu= new(8,graphic)

  do i = 0, dimsizes(plot)-1
    res@gsnLeftString  = label(i) +" (" + sprintf("%2.1f", var(i)) +"%)"
    res@gsnRightString = "APT: "+ sprintf("%3.1f", apt(i))+ " days"
    plot(i)  = gsn_csm_contour_map(wks,x(i,:,:),res)
    plotu(i) = gsn_csm_contour(wks,um,resc)
    overlay(plot(i), plotu(i))
  end do

  resP = True
  resP@gsnMaximize      = True
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02
  resP@gsnPanelYWhiteSpacePercent = 3

  resP@pmLabelBarOrthogonalPosF = -0.01
  resP@pmLabelBarWidthF         = 0.6

  resP@lbTitleString      = "u200 [m s~S~-1~N~]"
  resP@lbTitlePosition    = "Right"
  resP@lbTitleDirection   = "Across"
  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012

  resP@gsnPanelFigureStringsPerimOn = False
  gsn_panel(wks, plot,(/2,2,2,2/), resP)

  ;i=0
  ;do while(i.eq.0)
  ;  i=0
  ;end do