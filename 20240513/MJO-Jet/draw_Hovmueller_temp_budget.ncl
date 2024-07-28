;-----------------------------------------------------------
; draw Hovmueller. 
; zonal box 25-30N.
; shading: each terms.
; contour: dTdt.
; time: 2023-12-28
;-----------------------------------------------------------
;-----------------------------------------------------------

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  tdir   = ddir+"/temp_budget/anom_"+"dtdt"+"_r2.5.nc"
  udir   = ddir+"/temp_budget/anom_"+"udtdx"+"_r2.5.nc"
  vdir   = ddir+"/temp_budget/anom_"+"vdtdy"+"_r2.5.nc"
  wdir   = ddir+"/temp_budget/anom_"+"wdtdp"+"_r2.5.nc"
  Tdir   = ddir+"/ERA5_daily_anom_tvert_1979-2022_2.5x2.5.nc"

  latmin = 25
  latmax = 30
  lonmin = 0
  lonmax = 220
  tday = 24.*60.*60.
;-----------------------------------------------------------
; Main code.
;-----------------------------------------------------------

  ft = addfile(tdir,"r")
  fu = addfile(udir,"r")
  fv = addfile(vdir,"r")
  fw = addfile(wdir,"r")
  fT = addfile(Tdir,"r")

  t = ft->dtdt(:,{latmin:latmax},{lonmin:lonmax})
  u = fu->udtdx(:,{latmin:latmax},{lonmin:lonmax})*tday*(-1)
  v = fv->vdtdy(:,{latmin:latmax},{lonmin:lonmax})*tday*(-1)
  w = fw->wdtdp(:,{latmin:latmax},{lonmin:lonmax})*tday*(-1)
  T = fT->t(:,{latmin:latmax},{lonmin:lonmax})

  t := dim_avg_n_Wrap(t, 1)
  u := dim_avg_n_Wrap(u, 1)
  v := dim_avg_n_Wrap(v, 1)
  w := dim_avg_n_Wrap(w, 1)
  T := dim_avg_n_Wrap(T, 1)
  q  = t - w - u - v 
  copy_VarCoords(t, q)
  print("complete data selling "+ systemfunc("date"))

  time = t&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

;----------------------------------------------------------
;; selling MJO.

  delim = ","
  rmmf  = asciiread(rmmdir,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
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

;-------------------------------------------------------------

  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)

  amp   := amp(n1:n2)
  phase := phase(n1:n2)
  mon   := mon(n1:n2)
  time  := time(n1:n2)

  nt     = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  amp   := amp(nt)
  phase := phase(nt)

  delete(nt)
  ind_DJF = amp&timeind   ; indice of DJF between 1979-2022 DJF.
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
  llt = new((/6,41,dimsizes(t&lon)/),float)
  copy_VarCoords(t(0,:), llt(0,0,:))
  llt = 0.

  nt = ind(phase.eq.1.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 40
    ntt   = -15 + j
    nindd = nind + ntt
    llt(0,j,:) = dim_avg_n(T(nindd,:), 0)
    llt(1,j,:) = dim_avg_n(u(nindd,:), 0)
    llt(2,j,:) = dim_avg_n(v(nindd,:), 0)
    llt(3,j,:) = dim_avg_n(w(nindd,:), 0)
    llt(4,j,:) = dim_avg_n(q(nindd,:), 0)
    llt(5,j,:) = dim_avg_n(t(nindd,:), 0)
    delete([/ntt,nindd/])
  end do
  delete(nt)
  print("complete calc leadlag data. "+ systemfunc("date"))

;; resources

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

  rest@gsnRightString = "~F21~ave: "+latmin+"-"+latmax+"N"
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
  rest@cnMinLevelValF         = -2.1
  rest@cnMaxLevelValF         =  2.1
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
  resP@pmLabelBarOrthogonalPosF = 0.05
  resP@pmLabelBarParallelPosF   = -0.03
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

  wks = gsn_open_wks("eps", pdir+"/"+"Hovmueller_temp_budget")

  text1  = "~F21~(a) Temp"                 
  text2  = "~F21~(b) -(u~F34~6~F21~T~F21~/~F34~6~F21~x~F21~)'"
  text3  = "~F21~(d) -(v~F34~6~F21~T~F21~/~F34~6~F21~y~F21~)'"
  text4  = "~F21~(c) (~F33~w~F33~s~F21~)'"
  text5  = "~F21~(e) Q'/c~B~p~N~"
  label = (/text1,text2,text3,text4,text5/)

  print("Prepare labels done. "+ systemfunc("date"))

  plott = new(5,graphic)
  plotc = plott

  do i = 0, 4
    if(i.eq.4)then
      resc@cnInfoLabelOn = True
      resc@cnInfoLabelOrthogonalPosF = 0.1
      resc@cnInfoLabelFontHeightF    = 0.015
    end if
    rest@gsnLeftString = label(i)
    plott(i)  = gsn_csm_hov(wks, llt(i,:,:), rest)
    plotc(i) = gsn_csm_hov(wks, llt(5,:,:), resc)
    overlay(plott(i),plotc(i))
    ployy = (/0.,40./)

    resp = True
    ployx = (/60,60/)
    resp@gsLineColor = (/0.,79.,0.,1./)
    resp@gsLineDashPattern = 0
    resp@gsLineThicknessF  = 3.0
    dum1 = gsn_add_polyline(wks, plott(i), ployx, ployy, resp)
    overlay(plott(i),dum1)

    ployx = (/80,80/)
    dum2  = gsn_add_polyline(wks, plott(i), ployx, ployy, resp)
    overlay(plott(i),dum2)
  
    ployx = (/170,170/)
    resp@gsLineDashPattern = 1
    dum3  = gsn_add_polyline(wks, plott(i), ployx, ployy, resp)
    overlay(plott(i),dum3)
  
    ployx = (/190,190/)
    dum4  = gsn_add_polyline(wks, plott(i), ployx, ployy, resp)
    overlay(plott(i),dum4)
  end do

  gsn_panel(wks, (/plott(0),plott(1),plott(2),plott(3),plott(4)/), (/1,2,2/), resP)
  print("Draw done. "+ systemfunc("date"))

  ;i = 0
  ;do while (i.eq.0)
  ;  i = 0
  ;end do
