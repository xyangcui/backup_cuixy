;--------------------------------------------------------
; draw temperature anomaly in contour, MTG
; anomaly in shading and OLR anomaly in contour.
; time: 2023-11-22
;--------------------------------------------------------

;-----------------------------------------------------------

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  udir   = ddir+"/"+"ERA5_daily_anom_u250_1979-2022_2.5x2.5.nc"
  mtgdir = ddir+"/"+"ERA5_daily_anom_mtg_1979-2022_2.5x2.5.nc"
  tdir   = ddir+"/"+"ERA5_daily_anom_tvert_1979-2022_2.5x2.5.nc"
  ucdir  = ddir+"/"+"ERA5_daily_u250_1979-2022.nc"

  sigdir  = "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/Fig6_sig.nc"
  signdir = "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/Fig6_NAJS_sig.nc"
  sigldir = "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/Fig6_left_sig.nc"

  lats = 10
  latn = 65
  lonw = 0
  lone = 360

  latsl = 10
  latnl = 60
  lonwl = 0
  lonel = 210

  latsn = 10
  latnn = 60
  lonwn = 230
  lonen = 350
;--------------------------------------------------------
; u250
;--------------------------------------------------------
  f    = addfile(ucdir,"r")

  time = f->time
  TIME = cd_calendar(time, 0)
  mon  = toint(TIME(:,1))
  nt   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  u    = dim_avg_n_Wrap(f->u(nt,0,:,:),0)

;; for max latitude of each Jet.
  jlats = 10
  jlatn = 45
  jlonw = 10
  jlone = 70

  ut  = u({jlats:jlatn},{jlonw:jlone})
  lat = ut&latitude

  umax = dim_avg_n(ut,1)
  mm   = maxind(umax)
  jlatm = lat(mm)
  delete([/ut,lat,umax,mm/])

  jlats = 15
  jlatn = 50
  jlonw = 110
  jlone = 180

  ut  = u({jlats:jlatn},{jlonw:jlone})
  lat = ut&latitude

  umax = dim_avg_n(ut,1)
  mm   = maxind(umax)
  jlate = lat(mm)
  delete([/ut,lat,umax,mm/])

  jlats = 15
  jlatn = 60
  jlonw = 250
  jlone = 320

  ut  = u({jlats:jlatn},{jlonw:jlone})
  lat = ut&latitude
  delete(jlatn)
  umax = dim_avg_n(ut,1)
  mm   = maxind(umax)
  jlatn = lat(mm)
  delete([/ut,lat,umax,mm/])

;contour plot
  resc = True
  resc@gsnFrame     = False
  resc@gsnDraw      = False
  resc@gsnMaximize  = False
  resc@gsnAddCyclic = True
  resc@gsnContourNegLineDashPattern = 1
  resc@gsnLeftString      = ""
  resc@gsnRightString     = ""
  resc@cnMonoLineColor    = True
  resc@cnLinesOn          = True
  resc@cnLineLabelsOn     = False
  resc@cnLineColor        = "magenta"
  resc@cnLineThicknessF   = 3.0
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels      = (/30.,40.,46.,60.,70./)

  delete([/time,f,TIME,mon,nt/])

print("draw clim.u complete. "+ systemfunc("date"))
;--------------------------------------------------------

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

  fu = addfile(udir,"r")
  ft = addfile(tdir,"r")

  mtg  = fu->u(:,{lats:latn},{lonw:lone})
  t    = ft->t(:,{lats:latn},{lonw:lone})

  ul  = fu->u(:,{latsl:latnl},{lonwl:lonel})
  tl  = ft->t(:,{latsl:latnl},{lonwl:lonel})

  mtgn  = fu->u(:,{latsn:latnn},{lonwn:lonen})
  tn    = ft->t(:,{latsn:latnn},{lonwn:lonen})

print("complete data selling "+ systemfunc("date"))

  time = mtg&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  lat  = mtg&lat
  lon  = mtg&lon
;--------------------------------------------------------------
;; t sig.
;; sig
  fsig = addfile(sigdir,"r")
  m_value = fsig->m_value

  fsign = addfile(signdir,"r")
  m_valuen = fsign->m_value

  fsigl = addfile(sigldir,"r")
  m_valuel = fsigl->m_value
;-------------------------------------------------------
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

;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
;; MEJS phase 3 day 5 and EAJS phase 6 day 5.
  llmtg = new((/6,dimsizes(mtg&lat),dimsizes(mtg&lon)/),float)
  copy_VarCoords(mtg(0,:,:), llmtg(0,:,:))
  llmtg = 0.

  llt = new((/6,dimsizes(t&lat),dimsizes(t&lon)/),float)
  copy_VarCoords(t(0,:,:), llt(0,:,:))
  llt = 0.

  nt = ind(phase.eq.1.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 5
    ntt   = 0 + j*5
    nindd = nind + ntt
    llmtg(j,:,:) = dim_avg_n(mtg(nindd,:,:), 0)
    llt(j,:,:)   = dim_avg_n(t(nindd,:,:), 0)
    delete([/ntt,nindd/])
  end do
  delete([/nt,nind/])

;; for NAJS.
  llmtgn = new((/dimsizes(mtgn&lat),dimsizes(mtgn&lon)/),float)
  copy_VarCoords(mtgn(0,:,:), llmtgn(:,:))
  llmtgn = 0.

  lltn = new((/dimsizes(tn&lat),dimsizes(tn&lon)/),float)
  copy_VarCoords(tn(0,:,:), lltn(:,:))
  lltn = 0.

  nt = ind(phase.eq.8.and.amp.gt.1.)
  nind = ind_DJF(nt)
  ntt   = 5
  nindd = nind + ntt
  llmtgn = dim_avg_n(mtgn(nindd,:,:), 0)
  lltn   = dim_avg_n(tn(nindd,:,:), 0)
  delete([/nt,ntt,nindd,nind/])

;; for MEJS and EAJS.

  llume = new((/2,dimsizes(ul&lat),dimsizes(ul&lon)/),float)
  copy_VarCoords(ul(0,:,:), llume(0,:,:))
  llume = 0.

  lltme = new((/2,dimsizes(tl&lat),dimsizes(tl&lon)/),float)
  copy_VarCoords(tl(0,:,:), lltme(0,:,:))
  lltme = 0.

  nt = ind(phase.eq.3.and.amp.gt.1.)
  nind = ind_DJF(nt)
  ntt   = 5
  nindd = nind + ntt
  llume(0,:,:) = dim_avg_n(ul(nindd,:,:), 0)
  lltme(0,:,:) = dim_avg_n(tl(nindd,:,:), 0)
  delete([/nt,ntt,nindd,nind/])

  nt = ind(phase.eq.5.and.amp.gt.1.)
  nind = ind_DJF(nt)
  ntt   = 10
  nindd = nind + ntt
  llume(1,:,:) = dim_avg_n(ul(nindd,:,:), 0)
  lltme(1,:,:) = dim_avg_n(tl(nindd,:,:), 0)
  delete([/nt,ntt,nindd,nind/])

  print("complete calc leadlag data. "+ systemfunc("date"))

;; resources

;shading uwind plot.
  cmap = read_colormap_file("MPL_RdBu")

  rest = True
  rest@gsnFrame           = False
  rest@gsnDraw            = False
  rest@gsnMaximize        = False
  rest@gsnAddCyclic       = False
  rest@gsnMajorLatSpacing = 15
  rest@gsnMajorLonSpacing = 15
  rest@gsnRightStringFontHeightF = 0.015
  rest@gsnLeftStringFontHeightF  = 0.015
  rest@gsnLeftStringOrthogonalPosF  = 0.01
  rest@gsnRightStringOrthogonalPosF = 0.01

  rest@gsnRightString = "~F21~Temp_250hPa_U"

  rest@tmXTOn = False
  rest@tmYROn = False
  rest@tmXBLabelFont = 21
  rest@tmYLLabelFont = 21
  rest@tmXBLabelFontHeightF = 0.012
  rest@tmYLLabelFontHeightF = 0.012

  rest@mpFillOn        = False
  rest@mpFillDrawOrder = "Postdraw"

  rest@mpMinLatF   = 15
  rest@mpMaxLatF   = 55
  rest@mpMinLonF   = lonw
  rest@mpMaxLonF   = lone
  rest@mpLandFillColor = "white"

  rest@cnFillOn               = True
  rest@cnFillPalette          = cmap(::-1,:)
  rest@cnLinesOn              = False
  rest@cnLineLabelsOn         = False
  rest@cnLevelSelectionMode   = "Manuallevels"
  rest@cnMinLevelValF         = -7
  rest@cnMaxLevelValF         =  7
  rest@cnLevelSpacingF        =  1
  rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  rest@lbLabelBarOn           = False

; contour temp plot.
  resu = True
  resu@gsnDraw      = False
  resu@gsnMaximize  = False
  resu@gsnFrame     = False
  resu@gsnAddCyclic = False
  resu@gsnContourNegLineDashPattern = 1
  resu@cnMonoLineColor    = True
  resu@cnLinesOn          = True
  resu@cnLineLabelsOn     = False
  resu@cnLineColor        = "black"
  resu@cnLineThicknessF   = 3.0
  resu@cnInfoLabelPerimOn = False
  resu@cnInfoLabelOn      = False
  resu@cnLevelSelectionMode = "Explicitlevels"
  resu@cnLevels      = (/0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0/) 

  resP = True
  resP@gsnMaximize      = False
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.95
  resP@gsnPanelBottom   = 0.05
  resP@gsnPanelYWhiteSpacePercent = 5

  resP@pmLabelBarWidthF   = 0.12
  resP@pmLabelBarHeightF  = 0.7

  resP@lbBoxLinesOn      = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012
  resP@lbOrientation      = "Vertical" 

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
  sres@cnFillDotSizeF     = 0.003
  sres@gsnLeftString      =" "
  sres@gsnRightString     =" "
  sres@cnMonoFillColor    = True
  sres@cnFillColor        = "gray"

  opt = True
  opt@gsnShadeFillType   = "pattern"
  opt@gsnShadeHigh       = 17
  opt@gsnShadeFillScaleF = 2.0

  print("complete res. "+ systemfunc("date"))

;; plot.
  Label1 = (/"~F21~(c) phase 8 day 5","~F21~(b) phase 1 day 5",\
            "~F21~(a) phase 1 day 10","~F21~(d) phase 1 day 15",\
            "~F21~(e) phase 1 day 20","~F21~(c) phase 1 day 25"/)

  wks = gsn_open_wks("x11", pdir+"/Fig6.t_shading_u_contour")

  plott   = new(8,graphic)
  plotphi = plott
  dum1    = plott
  dum2    = plott
  dum3    = plott
  dum4    = plott
  dum     = plott

  plott_me = new(2,graphic)
  plotu_me = new(2,graphic)

  respoly = True
  respoly@gsLineThicknessF = 4.0
  respoly@gsLineColor      = "green"

;; NAJS
  n = 0
  rest@vpHeightF  = 0.3
  rest@gsnLeftString = "~F21~"+Label1(n)
  rest@gsnRightString = "~F21~NAJS_Temp_250hPa_U"
  delete(resu@cnLevels)
  resu@cnLevels      = (/-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6/)
  rest@mpMinLonF   = -120
  rest@mpMaxLonF   = -25
  rest@mpMaxLatF   =  55

  x = llmtgn
  x = mask(x,x.lt.m_valuen(0,0,:,:).or.x.gt.m_valuen(1,0,:,:),True)
  y  = lltn

  xt = lltn
  xt = 0.
  do i = 0, dimsizes(tn&lat)-1
    do j = 0, dimsizes(tn&lon)-1
      if(y(i,j).lt.m_valuen(0,1,i,j).or.y(i,j).gt.m_valuen(1,1,i,j)) then
        xt(i,j) = 2.
      else
        xt(i,j) = -2.
      end if
    end do
  end do

  y = mask(y,y.lt.m_valuen(0,1,:,:).or.y.gt.m_valuen(1,1,:,:),True)
    
  plott(n)   = gsn_csm_contour_map(wks, x, rest)
  plotphi(n) = gsn_csm_contour(wks, y, resu)

  dum3(n)  = gsn_add_polyline(wks, plott(n), (/-80,-40,-40,-80,-80/), \
                            (/20,20,28,28,20/), respoly)

  respoly@gsLineDashPattern = 1
  dum4(n)  = gsn_add_polyline(wks, plott(n), (/-100,-70,-70,-100,-100/), \
                            (/33,33,45,45,33/), respoly)
  plotu = gsn_add_polyline(wks,plott(n),(/-110,-40/),(/jlatn,jlatn/),respoly);gsn_csm_contour(wks, u, resc)
  dum(n) = gsn_csm_contour(wks,xt,sres)
  dum(n) = gsn_contour_shade(dum(n),1.9,-999.,opt)
  ;overlay(plott(n),dum(n))
  overlay(plott(n),plotu)
  overlay(plott(n),plotphi(n))

;; phase 3 day 5 
  delete([/x,y,xt/])
  delete(resu@cnLevels)
  resu@cnLevels      = (/0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0/)   
  rest@gsnLeftString = "~F21~(a) phase 3 day 5"

  rest@mpMaxLatF   = 45
  rest@mpMinLonF   = 30
  rest@mpMaxLonF   = 115
  rest@gsnRightString = "~F21~MEJS_Temp_250hPa_U"

  rest@vpHeightF  = 0.25

;; sig level.
  x = llume(0,:,:)
  x = mask(x,x.lt.m_valuel(0,0,:,:).or.x.gt.m_valuel(1,0,:,:),True)
  y  = lltme(0,:,:)
  xt = lltme(0,:,:)
  xt = 0.
  do i = 0, dimsizes(ul&lat)-1
    do j = 0, dimsizes(ul&lon)-1
      if(y(i,j).lt.m_valuel(0,1,i,j).or.y(i,j).gt.m_valuel(1,1,i,j)) then
        xt(i,j) = 2.
      else
        xt(i,j) = -2.
      end if
    end do
  end do

  y = mask(y,y.lt.m_valuel(0,1,:,:).or.y.gt.m_valuel(1,1,:,:),True)

  plotu_me(0) = gsn_csm_contour_map(wks,x,rest)
  plott_me(0) = gsn_csm_contour(wks, y, resu)

  plotts1 = gsn_csm_contour(wks,xt,sres)
  plotts1 = gsn_contour_shade(plotts1,1.9,-999.,opt)

  plotu = gsn_add_polyline(wks,plotu_me(0),(/10,70/),(/jlatm,jlatm/),respoly);gsn_csm_contour(wks, u, resc)

  ;overlay(plotu_me(0),plotts1)
  overlay(plotu_me(0),plotu)
  overlay(plotu_me(0),plott_me(0))

;; phase 5 day 10
  delete([/x,y,xt/])
  rest@gsnLeftString = "~F21~(b) phase 5 day 10"
  rest@mpMinLonF   = 125
  rest@mpMaxLonF   = 210
  rest@gsnRightString = "~F21~EAJS_Temp_250hPa_U"

  rest@vpHeightF  = 0.25

;; sig level.
  x = llume(1,:,:)
  x = mask(x,x.lt.m_valuel(0,2,:,:).or.x.gt.m_valuel(1,2,:,:),True)
  y  = lltme(1,:,:)
  xt = lltme(1,:,:)
  xt = 0.
  do i = 0, dimsizes(ul&lat)-1
    do j = 0, dimsizes(ul&lon)-1
      if(y(i,j).lt.m_valuel(0,3,i,j).or.y(i,j).gt.m_valuel(1,3,i,j)) then
        xt(i,j) = 2.
      else
        xt(i,j) = -2.
      end if
    end do
  end do

  y  = mask(y,y.lt.m_valuel(0,3,:,:).or.y.gt.m_valuel(1,3,:,:),True)

  plotu_me(1) = gsn_csm_contour_map(wks,x,rest)
  plott_me(1) = gsn_csm_contour(wks, y, resu)

  plotts2 = gsn_csm_contour(wks,xt,sres)
  plotts2 = gsn_contour_shade(plotts2,1.9,-999.,opt)

  plotu = gsn_add_polyline(wks,plotu_me(1),(/110,180/),(/jlate,jlate/),respoly);gsn_csm_contour(wks, u, resc)

  ;overlay(plotu_me(1),plotts2)
  overlay(plotu_me(1),plotu)
  overlay(plotu_me(1),plott_me(1))

  dumme = gsn_add_polyline(wks, plotu_me(0), (/70,100,100,70,70/), \
                                (/25,25,30,30,25/), respoly)
  dumea = gsn_add_polyline(wks, plotu_me(1), (/170,200,200,170,170/), \
                                (/25,25,33,33,25/), respoly)

  gsn_panel(wks, (/plotu_me(0),plotu_me(1),plott(0)/), (/1,1,1/), resP)

  i = 1
  do while (i.eq.1)
    i = 1
  end do