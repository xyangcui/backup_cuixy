;-----------------------------------------------------------------
;; draw general budget terms. 
;; MEJS box (10-70E,10-40N)
;; EAJS box (110-180E,15-50N)
;; time: 2023-12-13
;-----------------------------------------------------------------
  rad   = 4.*atan(1.)/180. 
  tday  = 24.*60.*60.

  mlatn = 35
  mlats = 25
  mlonw = 150
  mlone = 190

  ddir = "/home/sunming/data5/cuixy/global_jets/data/temp_budget"
  tdir = ddir+"/"+"anom_dtdt_r2.5.nc"
  udir = ddir+"/"+"anom_udtdx_r2.5.nc"
  vdir = ddir+"/"+"anom_vdtdy_r2.5.nc"
  wdir = ddir+"/"+"anom_wdtdp_r2.5.nc"
  t1dir = "/home/sunming/data5/cuixy/global_jets/data/ERA5_daily_anom_tvert_1979-2022_2.5x2.5.nc"
  pdireajs = "/home/sunming/data5/cuixy/global_jets/plot/Fig8_bar_EAJS"

;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
  ft  = addfile(tdir,"r")
  fu  = addfile(udir,"r")
  fv  = addfile(vdir,"r")
  fw  = addfile(wdir,"r")

  ft1 = addfile(t1dir,"r")
  t1  = ft1->t

  t = ft->dtdt
  t = t
  u = fu->udtdx
  u = u*tday
  v = fv->vdtdy
  v = v*tday
  w = fw->wdtdp
  w = w*tday

;residual
  r = t + (u + v + w) ;- r1)
  copy_VarCoords(t, r)
  u = -1*u
  v = -1*v
  w = -1*w;+r1
  lat = t&lat

  adir = ddir+"/anom_udtdx_seperate_r2.5.nc"

  fa  = addfile(adir,"r")

  ultl = fa->ultl
  ulti = fa->ulti
  ulth = fa->ulth
  uitl = fa->uitl
  uiti = fa->uiti
  uith = fa->uith
  uhtl = fa->uhtl
  uhti = fa->uhti
  uhth = fa->uhth

;---------------------------------------------------------------=
;; MJO phases.
  time = ft->time
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
;;MEJS.
  llmejs = new((/8,41,6/),float)
  llmejs = 0.
  tm = t(:,{mlats:mlatn},{mlonw:mlone})
  clatjs = cos(tm&lat*rad)

  t1js = wgt_areaave_Wrap(t1(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  tjs = wgt_areaave_Wrap(t(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  ujs = wgt_areaave_Wrap(u(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  vjs = wgt_areaave_Wrap(v(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  wjs = wgt_areaave_Wrap(w(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  rjs = wgt_areaave_Wrap(r(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)

  do i = 0, 7
    ip = i + 1
    nt = ind(phase.eq.ip.and.amp.gt.1.)
    nind = ind_DJF(nt)
    do j = 0, 40
      ntt   = -15 + j
      nindd = nind + ntt
      llmejs(i,j,0) = dim_avg_n(tjs(nindd), 0)
      llmejs(i,j,1) = dim_avg_n(ujs(nindd), 0)
      llmejs(i,j,2) = dim_avg_n(vjs(nindd), 0)
      llmejs(i,j,3) = dim_avg_n(wjs(nindd), 0)
      llmejs(i,j,4) = dim_avg_n(rjs(nindd), 0)
      llmejs(i,j,5) = dim_avg_n(t1js(nindd), 0)
      delete([/ntt,nindd/])
    end do
    delete([/ip,nt,nind/])
  end do

  delete([/clatjs,tjs,ujs,vjs,wjs,rjs/])
  print("calc MEJS data done. "+ systemfunc("date"))

;; line for dtdt
  line = new((/2,41/),float)
  line(0,:) = llmejs(2,:,0)
  line(1,:) = llmejs(2,:,5)
;; bar for other value. EAJS:27-32 day: 12-17 MEJS: day 2-7 17-22 ;;
  bar = new((/5/),float)
  bar = dim_avg_n(llmejs(2,13:17,0:4),0)

;; seperate
;; store lead-lag values for plotting. 
;;MEJS.
  lsmejs = new((/10/),float)
  lsmejs = 0.

  tm = t(:,{mlats:mlatn},{mlonw:mlone})
  clatjs = cos(tm&lat*rad)
  
  ujs = wgt_areaave_Wrap(u(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  ultljs = wgt_areaave_Wrap(ultl(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  ultijs = wgt_areaave_Wrap(ulti(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  ulthjs = wgt_areaave_Wrap(ulth(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  uitljs = wgt_areaave_Wrap(uitl(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  uitijs = wgt_areaave_Wrap(uiti(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  uithjs = wgt_areaave_Wrap(uith(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  uhtljs = wgt_areaave_Wrap(uhtl(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  uhtijs = wgt_areaave_Wrap(uhti(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  uhthjs = wgt_areaave_Wrap(uhth(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)

  nt = ind(phase.eq.3.and.amp.gt.1.)
  nind  = ind_DJF(nt)
  ntt   = -2
  nin = nind + ntt
  do j = 1, 4
    ntt   = -2 + j
    nindd = nind + ntt
    nin  := array_append_record(nin, nindd, 0)
  end do
  lsmejs(0) = dim_avg_n(ujs(nin), 0)
  lsmejs(1) = dim_avg_n(ultljs(nin), 0)*tday*(-1)
  lsmejs(2) = dim_avg_n(ultijs(nin), 0)*tday*(-1)
  lsmejs(3) = dim_avg_n(ulthjs(nin), 0)*tday*(-1)
  lsmejs(4) = dim_avg_n(uitljs(nin), 0)*tday*(-1)
  lsmejs(5) = dim_avg_n(uitijs(nin), 0)*tday*(-1)
  lsmejs(6) = dim_avg_n(uithjs(nin), 0)*tday*(-1)
  lsmejs(7) = dim_avg_n(uhtljs(nin), 0)*tday*(-1)
  lsmejs(8) = dim_avg_n(uhtijs(nin), 0)*tday*(-1)
  lsmejs(9) = dim_avg_n(uhthjs(nin), 0)*tday*(-1)

  delete([/clatjs,ujs/])
  print("calc MEJS data done. "+ systemfunc("date"))

;------------------------------------------------------
;; plot format: Line chart. about 5 lines. RGBA color.

  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@gsnRightStringFontHeightF  = 0.015
  res@gsnLeftStringFontHeightF   = 0.015
  res@gsnLeftStringOrthogonalPosF = -0.02
  res@gsnRightStringOrthogonalPosF = -0.02

  res@tmXTOn = False
  res@tmYROn = False

  res@tiXAxisString  = "~F21~Lag days"
  res@tiYAxisString  = "~F21~Temp tendency (K/day)"
  res@tmXBMode   = "Explicit"
  res@tmXBValues = (/ 0. , 5., 10., 15., 20., 25., 30., 35., 40. /)
  res@tmXBLabels = (/"-15","-10","-5","0","5","10","15","20","25" /)
  res@tmXBLabelFontHeightF = 0.012
  res@tiXAxisFontHeightF = 0.012

  res@tiYAxisFontHeightF = 0.012

  res@tmYLMode        = "Manual"
  res@tmYLMinorPerMajor = 0
  res@tmYLTickSpacingF  = 0.3
  res@tmYLLabelFontHeightF = 0.012

  res@tmYLMajorLengthF   = -0.01
  res@trYMaxF  = 1.5
  res@trYMinF  = -1.5
  res@gsnYRefLine = 0.
  res@gsnYRefLineDashPattern = 1

  res@xyLineThicknessF = 4.0

  color = (/"blue3"/)  ;"royalblue4"

  xlabel = ispan(0,40,1)

;; bar chart.
  resb = True
  resb@gsnDraw     = False
  resb@gsnFrame    = False
  resb@gsnMaximize = False

  resb@gsnLeftStringFontHeightF   = 0.015
  resb@tiYAxisString  = "~F21~Temp tendency (K/day)"
  resb@tmXTOn = False
  resb@tmYROn = False

  resb@tiYAxisFontHeightF = 0.012
  resb@tmYLMode           = "Manual"
  resb@tmYLMinorPerMajor  = 0
  resb@tmYLTickSpacingF   = 0.1
  resb@tmYLLabelFontHeightF = 0.012

  resb@tmXBMode    = "Explicit"
  resb@tmXBValues  = (/ 0.3 , 0.9, 1.4, 1.9, 2.4/)
  text1  = "~F21~(~F34~6~F21~T/~F34~6~F21~t~F21~)'"                 
  text2  = "~F21~-(u~F34~6~F21~T~F21~/~F34~6~F21~x~F21~)'"
  text3  = "~F21~-(v~F34~6~F21~T~F21~/~F34~6~F21~y~F21~)'"
  text4  = "~F21~(~F33~w~F33~s~F21~)'"
  text5  = "~F21~Q'/c~B~p~N~"

  resb@tmXBLabels  = (/text1,text2,text3,text4,text5/)
  resb@tmXBLabelFontHeightF = 0.012

  resb@gsnXYBarChart = True
  resb@gsnXYBarChartBarWidth = 0.2

  resb@tmYLMajorLengthF   = -0.01
  resb@gsnYRefLine        = 0.
  resb@trXMinF  = 0.0
  resb@trXMaxF  = 2.7
;----------------------------------------------------------------
;; Line plot
  wks = gsn_open_wks("x11", pdireajs)
  ;drawNDCGrid(wks)
  res@vpWidthF  = 0.6
  res@vpHeightF = 0.25
  res@vpXF = 0.2
  res@vpYF = 0.95
  res@trYMaxF  =  0.5
  res@trYMinF  = -0.4
  res@trXMaxF  =  30
  res@trXMinF  =  5
  res@tmYLTickSpacingF  = 0.1
  res@gsnLeftString  = "~F21~(a) (~F34~6~F21~T/~F34~6~F21~t~F21~)' & T' "
  res@gsnRightString = "~F21~ave: "+mlats+"-"+mlatn+"N, "+mlonw+"-"+mlone+"E"
  res@xyLineColor = color(0)
  res@xyDashPattern = 1

  resR = True
  resR@trYMaxF  =  3.0
  resR@trYMinF  = -3.0
  resR@xyLineThicknessF = 4.0
  resR@xyLineColor = color(0)
  resR@xyDashPattern = 0
  resR@tiYAxisString  = "~F21~Temp anomaly (K)"
  dum = gsn_csm_xy2(wks, xlabel, line(0,:),line(1,:), res,resR)

;;draw four lines.
  resl = True
  resl@gsLineColor = "gray64"
  resl@gsLineDashPattern = 2
  resl@gsLineThicknessF  = 4.0  
  dum21 = gsn_add_polyline(wks, dum, (/13,13/), (/-10,10/), resl)
  dum22 = gsn_add_polyline(wks, dum, (/17,17/), (/-10,10/), resl)

  rest               = True            
  rest@txFontHeightF = 0.014    

  rest@txFontColor = color(0)

  gsn_text_ndc(wks, "(~F34~6~F21~T/~F34~6~F21~t~F21~)'", 0.45, 0.9, rest)
  gsn_text_ndc(wks, "T'", 0.6, 0.9, rest)

  draw(dum)

;;draw Bar chart.

  resb@trYMaxF  =  1.0
  resb@trYMinF  = -1.0

  resb@tmYLTickSpacingF   = 0.3
  resb@vpXF = 0.2
  resb@vpYF = 0.63
  resb@vpWidthF  = 0.6
  resb@vpHeightF = 0.25
  resb@tmXBLabelsOn = True
  resb@gsnLeftString = "~F21~(b) budget terms"

  resb@gsnXYBarChartColors  = color(0)
  dumbar  = gsn_csm_xy(wks,(/ 0.3 , 0.9, 1.4, 1.9, 2.4/),bar,resb)

;;polyline.
  respoly = True
  respoly@gsLineColor = "black"
  respoly@gsLineDashPattern = 2
  respoly@gsLineThicknessF  = 4.0
  respoly@tfPolyDrawOrder   = "PostDraw"
  dum1 = gsn_add_polyline(wks,dumbar,(/0.6,0.6/),(/-10,10/),respoly)
;; shade color for lseajs(0) and lleajs(2).
  resea = True
  resea@gsFillIndex = 6
  resea@gsFillColor = "black"
  ;dum4 = gsn_add_polygon(wks, dumbar, (/0.8,1.0,1.0,0.8,0.8/), (/0.,0.,bar(1),bar(1),0./), resea)

  rest               = True            
  rest@txFontHeightF = 0.014    

  rest@txFontColor = color(0)
  undef("t1")
  t1 = gsn_add_text(wks, dumbar, sprintf("%3.2f",bar(0)), 0.3, bar(0)+0.08, rest)
  t2 = gsn_add_text(wks, dumbar, sprintf("%3.2f",bar(1)), 0.9, bar(1)+0.08, rest)
  t3 = gsn_add_text(wks, dumbar, sprintf("%3.2f",bar(2)), 1.4, bar(2)-0.08, rest)
  t4 = gsn_add_text(wks, dumbar, sprintf("%3.2f",bar(3)), 1.9, bar(3)-0.08, rest)
  t5 = gsn_add_text(wks, dumbar, sprintf("%3.2f",bar(4)), 2.4, bar(4)+0.08, rest)

  draw(dumbar)

  undef("resb")

;------------------------------------------------------
;; plot format: Line chart. about 5 lines. RGBA color.
;; dtdt in (0,0,0,1.). udtdx in (112,160,205,1.0)
;; vdtdy in (0,52,102,1.) r in (178,178,178)
  text2  = "~F21~-(u~F34~6~F21~T~F21~/~F34~6~F21~x~F21~)'"
;; bar chart.
  resb = True
  resb@gsnDraw     = False
  resb@gsnFrame    = False
  resb@gsnMaximize = False

  resb@tmXTOn = False
  resb@tmYROn = False

  resb@tiYAxisFontHeightF = 0.012
  resb@tiYAxisString      = "Temp tendency (K/day)"
  resb@tmYLMode           = "Manual"
  resb@tmYLMinorPerMajor  = 0
  resb@tmYLTickSpacingF   = 0.3
  resb@tmYLLabelFontHeightF = 0.012

  resb@tmXBMode    = "Explicit"
  resb@tmXBValues  = (/ 0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5./)
  resb@tmXBLabels  = (/text2,"~F21~LL","~F21~LI",\
                    "~F21~LH","~F21~IL","~F21~II","~F21~IH","~F21~HL","~F21~HI","~F21~HH"/)
  resb@tmXBLabelFontHeightF = 0.012

  resb@gsnXYBarChart = True
  resb@gsnXYBarChartBarWidth = 0.2

  resb@tmYLMajorLengthF   = -0.01
  resb@gsnYRefLine        = 0.
  resb@trXMinF  = 0.2
  resb@trXMaxF  = 5.3

  resb@vpWidthF  = 0.6
  resb@vpHeightF = 0.25
  resb@vpXF = 0.2
  resb@vpYF = 0.3
;----------------------------------------------------------------
;; MEJS
  resb@gsnLeftStringFontHeightF    = 0.016
  resb@gsnRightStringFontHeightF   = 0.016            
;;draw Bar chart.
  resb@trYMaxF  =  1.0
  resb@trYMinF  = -0.5
  resb@gsnLeftString  = "~F21~(c) "+text2
  resb@gsnRightString = ""
  dummebar = gsn_csm_xy(wks,fspan(0.5,5.,10),lsmejs,resb)
;;polyline.
  respoly = True
  respoly@gsLineColor = "gray64"
  respoly@gsLineDashPattern = 2
  respoly@gsLineThicknessF  = 5.0
  undef("dum3")
  dum3 = gsn_add_polyline(wks,dummebar,(/0.7,0.7/),(/10,-10/),respoly)
;; shade color for lseajs(0) and lleajs(2).
  resme = True
  resme@gsFillColor = color(0)
  undef("dum4")
  dum4 = gsn_add_polygon(wks, dummebar, (/1.35,1.65,1.65,1.35,1.35/), (/0.,0.,lsmejs(2),lsmejs(2),0./), resme)
  draw(dummebar)

  i = 1
  do while (i.eq.1)
    i = 1
  end do

