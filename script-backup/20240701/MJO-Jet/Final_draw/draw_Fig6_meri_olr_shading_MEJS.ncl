;-----------------------------------------------------------------
;; draw general budget terms. 
;; MEJS box (10-70E,10-40N)
;; EAJS box (110-180E,15-50N)
;; time: 2024-01-12
;; m-value(low:high,wme:vme:wea:vea,dates,level,lat)
;-----------------------------------------------------------------
;;---------------------------------------------------------------------
;; This function reads topographic data off a binary file as "short"
;; and converts it to float.  It also attaches lat/lon coordinate
;; variables.
  latn =  40
  lats = -20
  lonw =  100
  lone =  130

  ddir = "/home/sunming/data5/cuixy/global_jets/data/meri_circu"
  wdir = ddir + "/anom_"+"w"+"level_"+"EAJS"+"_r2.5.nc"
  vdir = ddir + "/anom_"+"v"+"level_"+"EAJS"+"_r2.5.nc"

  olrdir = "/home/sunming/data5/cuixy/global_jets/data/olr_daily_anom_1979-2022.nc"

  sdir =  "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/Fig7_meri_sig.nc"
  s1dir = "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/Fig7_olr_sig.nc"

  pdir = "/home/sunming/data5/cuixy/global_jets/plot/Fig6_meri_MEEAJS"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------

  fsig  = addfile(sdir,"r")
  mvalue = fsig->mvalue

  fsig1  = addfile(s1dir,"r")
  mvalue1 = fsig1->mvalue

  fw  = addfile(wdir,"r")
  fv  = addfile(vdir,"r")

  w = fw->w(:,:,{lats:latn},{lonw:lone})
  v = fv->v(:,:,{lats:latn},{lonw:lone})

  lat   = w&lat
  level = w&level

  w := dim_avg_n_Wrap(w, 3)
  v := dim_avg_n_Wrap(v, 3)

  folr = addfile(olrdir,"r")

  olr = folr->olr(:,{lats:latn},{0:180})
;---------------------------------------------------------------=
;; MJO phases.
  time = fw->time
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

  nt = ind(phase.eq.2.and.amp.gt.1.)
  nind = ind_DJF(nt)
;----------------------------------------------------------------------
  w_pvalue = new((/4,dimsizes(w&level),dimsizes(w&lat)/),float)
  w_pvalue = 0.
  v_pvalue = w_pvalue

  llolr = new((/4,dimsizes(olr&lat),dimsizes(olr&lon)/),float)
  llolr = 0.

  do j = 0, 3
    ntt   = -4 + j*2
    nindd = nind + ntt
    w_pvalue(j,:,:) = dim_avg_n(w(nindd,:,:), 0)
    v_pvalue(j,:,:) = dim_avg_n(v(nindd,:,:), 0)
    llolr(j,:,:) = dim_avg_n(olr(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do

  copy_VarCoords(w(0,:,:),w_pvalue(0,:,:))
  copy_VarCoords(v(0,:,:),v_pvalue(0,:,:))
  copy_VarCoords(olr(0,:,:), llolr(0,:,:))

  llolr = smth9_Wrap(llolr, 0.5, 0.25, False)
  print("complete calc leadlag data. "+ systemfunc("date"))
;; res for meri circulation.
  cmap = read_colormap_file("MPL_RdBu")

  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@gsnMaximize  = False
  res@gsnRightStringFontHeightF  = 0.015
  res@gsnLeftStringFontHeightF   = 0.015
  res@gsnRightStringOrthogonalPosF = 0.01
  res@gsnLeftStringOrthogonalPosF  = 0.01

  res@trXMinF = 0

  res@tmXTOn = False
  res@tmYROn = False
  res@gsnMajorLatSpacing      = 15
  res@gsnPresHgtHeightLabelOn = False
  res@tmYLLabelFont  = 21
  res@tmXBLabelFont  = 21
  res@tmYLLabelFontHeightF = 0.01
  res@tmXBLabelFontHeightF = 0.01
  res@tiYAxisFont = 21
  res@tiYAxisFontHeightF   = 0.012

  res@cnFillPalette        = cmap(::-1,:)
  res@cnFillOn             = True                
  res@cnLevelSelectionMode = "ManualLevels"      
  res@cnMinLevelValF       = -3              
  res@cnMaxLevelValF       =  3              
  res@cnLevelSpacingF      =  0.5  
  res@lbLabelBarOn         = False  
  res@lbOrientation        = "Vertical" 
  res@cnLinesOn = False    
  res@cnInfoLabelOn  = False  
  res@cnLineLabelsOn = False         
  res@cnLabelBarEndStyle     = "ExcludeOuterBoxes"

  res@vcGlyphStyle              = "CurlyVector"
  res@vcRefAnnoOrthogonalPosF   = -0.165                
  res@vcRefMagnitudeF           = 3       
  res@vcRefLengthF              = 0.1 
  res@vcLineArrowThicknessF     = 1.5       
  res@vcMapDirection            = False 
  res@vcRefAnnoOn               = False
  res@vcRefAnnoString2On        = False
  res@vcRefAnnoString1On        = False
  res@vcRefAnnoPerimOn          = False

  res@tmXBOn = True

  ress = res
  ress@cnFillOn = False
  ress@cnConstFEnableFill           = True
  ress@cnConstFLabelBackgroundColor = -1
  ress@cnConstFLabelOn  = False

  ress@vcRefAnnoOn             = True
  ress@vcRefAnnoString1On      = True
  ress@vcRefAnnoString1        = "$VMG$"
  ress@vcRefAnnoFontHeightF    =  0.015
  ress@vcRefAnnoOrthogonalPosF =  -0.3
  ress@vcRefAnnoParallelPosF   =  0.98

  ress@tmXBOn = True

  sres   = True 
  sres@cnLineDrawOrder    = "PostDraw"                    
  sres@gsnDraw            = False                          
  sres@gsnFrame           = False                        
  sres@cnLineLabelsOn     = False                              
  sres@cnLinesOn          = False                           
  sres@cnInfoLabelOn      = False
  sres@cnFillOn           = False                       
  sres@cnInfoLabelOn      = False
  sres@cnFillScaleF       = 0.5
  sres@cnFillDotSizeF     = 0.0025
  sres@gsnLeftString      =" "
  sres@gsnRightString     =" "
  sres@cnMonoFillColor    = True
  sres@cnFillColor        = "gray"
  sres@gsnPresHgtHeightLabelOn = False 

  opt = True
  opt@gsnShadeFillType   = "pattern"
  opt@gsnShadeHigh       = 17
  opt@gsnShadeFillScaleF = 1.

  respoly = True
  respoly@gsLineThicknessF = 3.5
  respoly@gsLineColor      = "green3"

  wks = gsn_open_wks("x11",pdir)

  plot  = new(4,graphic)
  plott = plot
  plotu = plot
  dum   = plot
  ploto = plot

  blankv = w_pvalue(0,:,:)
  blankv = 0.
;; m-value(low:high,wme:vme:wea:vea,dates,level,lat)
  polyx = (/5,15,15,5,5/)
  polyy = (/1100,1100,100,100,1100/)
  res@gsnRightString  = "~F21~ave:"+lonw+"-"+lone+"E"

  label = (/"(a) day-10","(b) day-7","(c) day-4","(d) day-1"/)

  res@vpWidthF  = 0.37
  res@vpHeightF = 0.23

  print("draw plot1 begin. " + systemfunc("date"))
  i = 0
  res@gsnLeftString   = label(i)
;; w,v sig.
  vp05 = mvalue(1,i,0,:,:)
  vp95 = mvalue(1,i,1,:,:)
  wp05 = mvalue(0,i,0,:,:)
  wp95 = mvalue(0,i,1,:,:)

;;for vector.
  x  = w_pvalue(i,:,:)
  y  = v_pvalue(i,:,:)
  wp = w_pvalue(i,:,:)
  vp = v_pvalue(i,:,:)

  y  = mask(y,vp.lt.vp05.or.vp.gt.vp95.or.wp.lt.wp05.or.wp.gt.wp95, True) ;\
  x  = mask(x,vp.lt.vp05.or.vp.gt.vp95.or.wp.lt.wp05.or.wp.gt.wp95, True)

  wp = wp*(-100)
  x  = x*(-100)
;; sig.
  wvalue = w_pvalue(i,:,:)

  xt = wvalue
  xt = 0.
  do k = 0, dimsizes(level)-1
    do j = 0, dimsizes(lat)-1
      if(wvalue(k,j).lt.wp05(k,j).or.wvalue(k,j).gt.wp95(k,j)) then
        xt(k,j) = 2.
      else
        xt(k,j) = -2.
      end if
    end do
  end do

  wvalue = wvalue*100

  dum(i) = gsn_csm_pres_hgt(wks, xt, sres )
  dum(i) = gsn_contour_shade(dum(i),1.9,-999.,opt)

  res@vcLineArrowColor   = "gray"
  plot(i)  = gsn_csm_pres_hgt_vector(wks,wvalue,vp,wp, res)

  ress@vcLineArrowColor   = "black"
  plotu(i) = gsn_csm_pres_hgt_vector(wks,blankv,y,x, ress)

  respoly@gsLineColor = "green3"
  dum1  = gsn_add_polyline(wks, plot(i), polyx, polyy, respoly)

  overlay(plot(i),plotu(i))
  overlay(plot(i),dum(i))

  delete([/vp05,vp95,wp05,wp95,wp,vp,wvalue/])
  print("draw plot1 done. " + systemfunc("date"))


  print("draw plot2 begin. " + systemfunc("date"))
  i = 1
  res@gsnLeftString   = label(i)
;; w,v sig.
  vp05 = mvalue(1,i,0,:,:)
  vp95 = mvalue(1,i,1,:,:)
  wp05 = mvalue(0,i,0,:,:)
  wp95 = mvalue(0,i,1,:,:)

  x  = w_pvalue(i,:,:)
  y  = v_pvalue(i,:,:)
  wp = w_pvalue(i,:,:)
  vp = v_pvalue(i,:,:)

  y  = mask(y,vp.lt.vp05.or.vp.gt.vp95.or.wp.lt.wp05.or.wp.gt.wp95, True) 
  x  = mask(x,vp.lt.vp05.or.vp.gt.vp95.or.wp.lt.wp05.or.wp.gt.wp95, True)

  wp = wp*(-100)
  x  = x*(-100)
;; sig.
  wvalue = w_pvalue(i,:,:)

  xt = wvalue
  xt = 0.
  do k = 0, dimsizes(level)-1
    do j = 0, dimsizes(lat)-1
      if(wvalue(k,j).lt.wp05(k,j).or.wvalue(k,j).gt.wp95(k,j)) then
        xt(k,j) = 2.
      else
        xt(k,j) = -2.
      end if
    end do
  end do

  wvalue = wvalue*100

  dum(i) = gsn_csm_pres_hgt(wks, xt, sres )
  dum(i) = gsn_contour_shade(dum(i),1.9,-999.,opt)

  res@vcLineArrowColor   = "gray"
  plot(i)  = gsn_csm_pres_hgt_vector(wks,wvalue,vp,wp, res)

  ress@vcLineArrowColor   = "black"
  plotu(i) = gsn_csm_pres_hgt_vector(wks,blankv,y,x, ress)

  respoly@gsLineColor = "green3"
  dum2  = gsn_add_polyline(wks, plot(i), polyx, polyy, respoly)

  overlay(plot(i),plotu(i))
  overlay(plot(i),dum(i))

  delete([/vp05,vp95,wp05,wp95,wp,vp,wvalue/])
  print("draw plot2 done. " + systemfunc("date"))

  print("draw plot3 begin. " + systemfunc("date"))
  i = 2
  res@gsnLeftString   = label(i)
;; w,v sig.
  vp05 = mvalue(1,i,0,:,:)
  vp95 = mvalue(1,i,1,:,:)
  wp05 = mvalue(0,i,0,:,:)
  wp95 = mvalue(0,i,1,:,:)

  x  = w_pvalue(i,:,:)
  y  = v_pvalue(i,:,:)
  wp = w_pvalue(i,:,:)
  vp = v_pvalue(i,:,:)

  y  = mask(y,vp.lt.vp05.or.vp.gt.vp95.or.wp.lt.wp05.or.wp.gt.wp95, True) ;\
  x  = mask(x,vp.lt.vp05.or.vp.gt.vp95.or.wp.lt.wp05.or.wp.gt.wp95, True)

  wp = wp*(-100)
  x  = x*(-100)
;; sig.
  wvalue = w_pvalue(i,:,:)

  xt = wvalue
  xt = 0.
  do k = 0, dimsizes(level)-1
    do j = 0, dimsizes(lat)-1
      if(wvalue(k,j).lt.wp05(k,j).or.wvalue(k,j).gt.wp95(k,j)) then
        xt(k,j) = 2.
      else
        xt(k,j) = -2.
      end if
    end do
  end do

  wvalue = wvalue*100

  dum(i) = gsn_csm_pres_hgt(wks, xt, sres )
  dum(i) = gsn_contour_shade(dum(i),1.9,-999.,opt)

  res@vcLineArrowColor   = "gray"
  plot(i)  = gsn_csm_pres_hgt_vector(wks,wvalue,vp,wp, res)

  ress@vcLineArrowColor   = "black"
  plotu(i) = gsn_csm_pres_hgt_vector(wks,blankv,y,x, ress)

  respoly@gsLineColor = "green3"
  dum3  = gsn_add_polyline(wks, plot(i), polyx, polyy, respoly)

  overlay(plot(i),plotu(i))
  overlay(plot(i),dum(i))

  delete([/vp05,vp95,wp05,wp95,wp,vp,wvalue/])
  print("draw plot3 done. " + systemfunc("date"))

  print("draw plot4 begin. " + systemfunc("date"))
  i = 3
  res@gsnLeftString   = label(i)
;; w,v sig.
  vp05 = mvalue(1,i,0,:,:)
  vp95 = mvalue(1,i,1,:,:)
  wp05 = mvalue(0,i,0,:,:)
  wp95 = mvalue(0,i,1,:,:)

  x  = w_pvalue(i,:,:)
  y  = v_pvalue(i,:,:)
  wp = w_pvalue(i,:,:)
  vp = v_pvalue(i,:,:)

  y  = mask(y,vp.lt.vp05.or.vp.gt.vp95.or.wp.lt.wp05.or.wp.gt.wp95, True) ;\
  x  = mask(x,vp.lt.vp05.or.vp.gt.vp95.or.wp.lt.wp05.or.wp.gt.wp95, True)

  wp = wp*(-100)
  x  = x*(-100)
;; sig.
  wvalue = w_pvalue(i,:,:)

  xt = wvalue
  xt = 0.
  do k = 0, dimsizes(level)-1
    do j = 0, dimsizes(lat)-1
      if(wvalue(k,j).lt.wp05(k,j).or.wvalue(k,j).gt.wp95(k,j)) then
        xt(k,j) = 2.
      else
        xt(k,j) = -2.
      end if
    end do
  end do

  wvalue = wvalue*100

  dum(i) = gsn_csm_pres_hgt(wks, xt, sres )
  dum(i) = gsn_contour_shade(dum(i),1.9,-999.,opt)

  res@vcLineArrowColor   = "gray"
  plot(i)  = gsn_csm_pres_hgt_vector(wks,wvalue,vp,wp, res)

  ress@vcLineArrowColor   = "black"
  plotu(i) = gsn_csm_pres_hgt_vector(wks,blankv,y,x, ress)

  respoly@gsLineColor = "green3"
  dum4  = gsn_add_polyline(wks, plot(i), polyx, polyy, respoly)

  overlay(plot(i),plotu(i))
  overlay(plot(i),dum(i))

  delete([/vp05,vp95,wp05,wp95,wp,vp,wvalue/])
  print("draw plot4 done. " + systemfunc("date"))

  resP = True
  resP@gsnMaximize      = False
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02
  resP@gsnPanelRight    = 0.55

  resP@pmLabelBarWidthF = 0.35

  resP@lbTitleFontHeightF = 0.012
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012

  resP@gsnPanelLabelBar = True
  resP@lbTitleOn = True
  resP@lbTitleString = "[10~S~-2~N~ Pa s~S~-1~N~]"
  resP@lbTitlePosition  = "Right"
  resP@lbTitleDirection = "Across"
  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.013
  resP@lbTitleStringParallelPosF = -0.1

  resP@gsnPanelXWhiteSpacePercent = 3

  gsn_panel(wks, plot, (/1,1,1,1/), resP)
;; res for olr.
  
  label1 = (/"(e) day-10","(f) day-7","(g) day-4","(h) day-1"/)
  cmap1 = read_colormap_file("BlueWhiteOrangeRed")

  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = False
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 30
  rest@gsnRightString     = "~F21~OLR"
  rest@gsnLeftStringFontHeightF  = 0.03
  rest@gsnRightStringFontHeightF = 0.03
  rest@gsnRightStringOrthogonalPosF = 0.01
  rest@gsnLeftStringOrthogonalPosF  = 0.01

  rest@tmXTOn = False
  rest@tmYROn = False
  rest@tmYLLabelFontHeightF = 0.02
  rest@tmXBLabelFontHeightF = 0.02

  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"
  rest@mpMinLatF   = -15
  rest@mpMaxLatF   =  40
  rest@mpMinLonF   =  0
  rest@mpMaxLonF   =  90
  rest@mpLandFillColor = "white"

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap1;(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -15
  rest@cnMaxLevelValF         =  15
  rest@cnLevelSpacingF        =  3
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = False

  delete(xt)
  do n = 0,3
    rest@gsnLeftString = "~F21~"+label1(n)
    ;xt = llolr(n,:,:)
    ;xt = 0.
    ;do i = 0, dimsizes(olr&lat)-1
    ;  do j = 0, dimsizes(olr&lon)-1
    ;    if(llolr(n,i,j).lt.mvalue1(n,0,i,j).or.llolr(n,i,j).gt.mvalue1(n,1,i,j)) then
    ;      xt(i,j) = 2.
    ;    else
    ;      xt(i,j) = -2.
    ;    end if
    ;  end do
    ;end do
    ploto(n) = gsn_csm_contour_map(wks, llolr(n,:,:), rest)
    dum1  = gsn_add_polyline(wks, ploto(n), (/lonw,lone,lone,lonw,lonw/), \
                                (/5,5,15,15,5/), respoly)
    ;dum(n) = gsn_csm_contour(wks,xt,sres)
    ;dum(n) = gsn_contour_shade(dum(n),1.9,-999.,opt)
    ;overlay(ploto(n),dum(n))
    overlay(ploto(n), dum1)
  end do

  delete(resP)
  resP = True
  resP@gsnMaximize      = False
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02
  resP@gsnPanelLeft     = 0.4

  resP@pmLabelBarWidthF         = 0.35

  resP@lbTitleFontHeightF = 0.012
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012

  resP@gsnPanelLabelBar = True
  resP@lbTitleOn = True
  resP@lbTitleString = "[W m~S~-2~N~]"
  resP@lbTitlePosition  = "Right"
  resP@lbTitleDirection = "Across"
  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.013
  resP@lbTitleStringParallelPosF = -0.1

  resP@gsnPanelXWhiteSpacePercent = 3

  gsn_panel(wks, ploto, (/1,1,1,1/), resP)

  i = 0
  do while (i.eq.0)
    i = 0
  end do
