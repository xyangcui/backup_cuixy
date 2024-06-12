;-----------------------------------------------------------------
;; draw general budget terms. 
;; MEJS box (10-70E,10-40N)
;; EAJS box (110-180E,15-50N)
;; time: 2023-12-13
;-----------------------------------------------------------------
  rad   = 4.*atan(1.)/180. 
  tday  = 24.*60.*60.

  lats = 10
  latn = 50
  lonw = 100
  lone = 210

  ddir = "/home/sunming/data5/cuixy/global_jets/data/temp_budget"
  adir = ddir+"/anom_vdtdy_seperate_r2.5.nc"
  pdir = "/home/sunming/data5/cuixy/global_jets/plot/EAJS_budget"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------

  f  = addfile(adir,"r")

  ultl1 = f->vltl(:,{latmin:latmax},{lonmin:lonmax})
  printVarSummary(ultl1)
  lon  = ultl1&lon

  ultl = f->vltl(:,{latmin:latmax},{lonmin:lonmax})
  ulti = f->vlti(:,{latmin:latmax},{lonmin:lonmax})
  ulth = f->vlth(:,{latmin:latmax},{lonmin:lonmax})
  uitl = f->vitl(:,{latmin:latmax},{lonmin:lonmax})
  uiti = f->viti(:,{latmin:latmax},{lonmin:lonmax})
  uith = f->vith(:,{latmin:latmax},{lonmin:lonmax})
  uhtl = f->vhtl(:,{latmin:latmax},{lonmin:lonmax})
  uhti = f->vhti(:,{latmin:latmax},{lonmin:lonmax})
  uhth = f->vhth(:,{latmin:latmax},{lonmin:lonmax})


  ultl := dim_avg_n_Wrap(ultl, 1)
  ulti := dim_avg_n_Wrap(ulti, 1)
  ulth := dim_avg_n_Wrap(ulth, 1)
  uitl := dim_avg_n_Wrap(uitl, 1)
  uiti := dim_avg_n_Wrap(uiti, 1)
  uith := dim_avg_n_Wrap(uith, 1)
  uhtl := dim_avg_n_Wrap(uhtl, 1)
  uhti := dim_avg_n_Wrap(uhti, 1)
  uhth := dim_avg_n_Wrap(uhth, 1)

;---------------------------------------------------------------=
;; MJO phases.
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
  delim  = ","
  rmmf   = asciiread(rmmdir,-1,"string")
  iyear  = str_get_field(rmmf,1,delim)
  imon   = str_get_field(rmmf,2,delim)
  iday   = str_get_field(rmmf,3,delim)
  phaset = str_get_field(rmmf,6,delim)
  ampt   = str_get_field(rmmf,7,delim)

  iyear  := toint(iyear(1::))
  imon   := toint(imon(1::))
  iday   := toint(iday(1::))
  phaset := toint(phaset(1::))
  ampt   := tofloat(ampt(1::))

  amp   = new((/dimsizes(time)/),float)
  phase = new((/dimsizes(time)/),integer)
  amp   = 0.
  phase = 0

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
    amp(i)   = ampt(nt)
    phase(i) = phaset(nt)
    delete(nt)
  end do

  nn = ispan(0, dimsizes(time)-1, 1)
  amp!0   = "timeind"
  phase!0 = "timeind"
  amp&timeind   = nn
  phase&timeind = nn
  delete([/nn,ampt,phaset,iyear,imon,iday/])

  print("complete RMM selling "+ systemfunc("date"))

  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.2.and.day.eq.28)

  amp   := amp(n1:n2)
  phase := phase(n1:n2)
  mon   := mon(n1:n2)
  time  := time(n1:n2)

  nt     = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  amp   := amp(nt)
  phase := phase(nt)

  delete(nt)
  ind_DJF = amp&timeind   ; indice of DJF between 1979-2022 DJF.


;; store lead-lag values for plotting. 
  lltt = new((/10,5,dimsizes(lon)/),float)
  
  lltt = 0.

  nt = ind(phase.eq.1.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 4
    ntt   = -2 + j
    nindd = nind + ntt
    lltt(0,j,:) = dim_avg_n(ultl(nindd,:), 0)
    lltt(1,j,:) = dim_avg_n(ulti(nindd,:), 0)
    lltt(2,j,:) = dim_avg_n(ulth(nindd,:), 0)
    lltt(3,j,:) = dim_avg_n(uitl(nindd,:), 0)
    lltt(4,j,:) = dim_avg_n(uiti(nindd,:), 0)
    lltt(5,j,:) = dim_avg_n(uith(nindd,:), 0)
    lltt(6,j,:) = dim_avg_n(uhtl(nindd,:), 0)
    lltt(7,j,:) = dim_avg_n(uhti(nindd,:), 0)
    lltt(8,j,:) = dim_avg_n(uhth(nindd,:), 0)
    delete([/ntt,nindd/])
  end do

  llt = dim_avg_n(lltt,0)

  do i = 0,8
    llt(9,:,:) = llt(9,:,:) + llt(i,:,:)
  end do
  ;llt(9,:,:) = smth9(llt(9,:,:), 0.5, 0.0, False)
  llt = llt*tday*(-1.)
  print("calc MEJS data done. "+ systemfunc("date"))
  copy_VarCoords(ultl1(0,0,:), llt(0,0,:))
;------------------------------------------------------
;; plot format: Line chart. about 5 lines. RGBA color.
;; dtdt in (0,0,0,1.). udtdx in (112,160,205,1.0)
;; vdtdy in (0,52,102,1.) r in (178,178,178)


res = True
res@gsnDraw  = False
res@gsnFrame = False
;res@gsnRightStringFontHeightF  = 0.012
;res@gsnLeftStringFontHeightF   = 0.015

res@tmXTOn = False
res@tmYROn = False
;res@tmYRBorderOn = False
;res@tmXTBorderOn = False

res@tiXAxisString  = "Lag days"
res@tmXBMode   = "Explicit"
res@tmXBValues = (/ 0. , 5., 10., 15., 20., 25., 30., 35., 40. /)
res@tmXBLabels = (/"-15","-10","-5","0","5","10","15","20","25" /)

res@tiYAxisFontHeightF = 0.012
res@tiYAxisString = "Temperature variation (K/day)"
res@tmYLMode        = "Manual"
res@tmYLMinorPerMajor = 0
res@tmYLTickSpacingF  = 0.3
res@tmYLLabelFontHeightF = 0.015
res@tmYLMajorLengthF   = -0.01
res@trYMaxF  = 1.5
res@trYMinF  = -1.5
; t u v w r
color = (/"black","red","blue","green","purple"/)
xlabel = ispan(0,40,1)
;----------------------------------------------------------------
;shading T plot
  cmap = read_colormap_file("MPL_RdBu")

  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnMajorLonSpacing  = 30
  ;rest@gsnCenterStringFontHeightF = 0.025 

  rest@vpWidthF  = 0.45
  rest@vpHeightF = 0.2

  rest@gsnRightString = "~F21~ave: "+latmin+"-"+latmax+"N unit: K"
  rest@gsnRightStringFontHeightF  = 0.012
  rest@gsnLeftStringFontHeightF   = 0.015

  rest@tmXTOn = False
  rest@tmYROn = False

  rest@tiYAxisString  = "~F21~Lag days"

  rest@tmYLMode   = "Explicit"
  rest@tmYLValues = (/ 0. , 5., 10., 15., 20., 25., 30., 35., 40. /)
  rest@tmYLLabels = (/"-15","-10","-5","0","5","10","15","20","25" /)

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnInfoLabelOn          = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -1.8
  rest@cnMaxLevelValF         =  1.8
  rest@cnLevelSpacingF        =  0.3

  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = False
  ;rest@lbTitleString          = "~F21~unit: K"
  ;rest@lbTitlePosition        = "Right"
  ;rest@lbTitleDirection       = "Across"
  ;rest@lbTitleFontHeightF     = 0.015
  rest@lbBoxLinesOn           = True
  rest@lbBoxSeparatorLinesOn  = True
  rest@lbLabelFontHeightF     = 0.012

; contour phi plot.
  resc = True
  resc@gsnDraw      = False
  resc@gsnMaximize  = False
  resc@gsnFrame     = False
  resc@gsnAddCyclic = False
  resc@gsnContourNegLineDashPattern = 1
  resc@cnMonoLineColor    = True
  resc@cnLinesOn          = True
  resc@cnLineLabelsOn     = False
  resc@cnLineColor        = "black"
  resc@cnLineThicknessF   = 2.5
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  ;resc@cnLevelSelectionMode = "Explicitlevels"
  ;resc@cnLevels      = (/-10,-8,-6,-4,-2,2,4,6,8,10/)

  resP = True
  resP@gsnMaximize      = False
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02

  resP@gsnPanelYWhiteSpacePercent = 5
  resP@pmLabelBarWidthF         = 0.
  ;resP@pmLabelBarOrthogonalPosF = 0.05
  ;resP@pmLabelBarParallelPosF   = -0.1
  resP@pmLabelBarWidthF         = 0.6

  ;resP@lbTitleString     = "~F21~unit: K/day"
  resP@lbTitlePosition   = "Bottom"
  ;resP@lbTitleDirection  = "Across"
  resP@lbTitleFontHeightF= 0.015
  resP@lbBoxLinesOn      = True
  resP@lbBoxSeparatorLinesOn  = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012

  print("complete res. "+ systemfunc("date"))

;; plot.

  wks = gsn_open_wks("x11", pdir+"/"+"Fig10.hovmueller_temp_budget")

  label = (/"ultl","ulti","ulth","uitl","uiti","uith","uhtl","uhti","uhth","total"/)

  print("Prepare labels done. "+ systemfunc("date"))

  plott = new(10,graphic)
  plotc = plott

  do i = 0, 9
    if(i.eq.4)then
      resc@cnInfoLabelOn = True
      resc@cnInfoLabelString = "~F21~CONTOUR FROM -.25 TO .4 BY .05 K/day"
      resc@cnInfoLabelOrthogonalPosF = 0.1
      resc@cnInfoLabelFontHeightF    = 0.015
    end if
    rest@gsnLeftString = label(i)
    if (i.eq.0)then
      rest@gsnRightString = "~F21~ave: "+latmin+"-"+latmax+"N unit: K"
    else
      rest@gsnRightString = "~F21~ave: "+latmin+"-"+latmax+"N unit: K/day"
    end if
    
    plott(i)  = gsn_csm_hov(wks, llt(i,:,:), rest)

    resp = True
    ployx = (/70,90,90,70,70/)
    ployy = (/18,18,23,23,18/)
    resp@gsLineColor = "green4"
    resp@gsLineDashPattern = 0
    resp@gsLineThicknessF  = 3.5
    dum1 = gsn_add_polyline(wks, plott(i), ployx, ployy, resp)
    overlay(plott(i),dum1)
  
    ployx = (/170,190,190,170,170/)
    ployy = (/27,27,33,33,27/)
    resp@gsLineDashPattern = 1
    dum3  = gsn_add_polyline(wks, plott(i), ployx, ployy, resp)
    overlay(plott(i),dum3)
  
  end do

  gsn_panel(wks, plott, (/3,3,3/), resP)
  print("Draw done. "+ systemfunc("date"))

  i = 0
  do while (i.eq.0)
    i = 0
  end do
















