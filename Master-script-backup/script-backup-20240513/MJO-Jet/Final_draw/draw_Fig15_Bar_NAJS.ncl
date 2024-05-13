;-----------------------------------------------------------------
;; draw general budget terms. 
;; time: 2023-12-13
;; + temp: (-70- -40 W, 20-28 N) day -5 - -10
;; - temp: (-100 - -70 W, 33-45N)
;-----------------------------------------------------------------
  rad   = 4.*atan(1.)/180. 
  tday  = 24.*60.*60.
  mlatn =  28
  mlats =  15
  mlonw = -70
  mlone = -40

  elatn =  50
  elats =  35
  elonw = -100
  elone = -70

  latn = mlatn
  lats = mlats
  lonw = mlonw
  lone = mlone

  color = (/"blue3","maroon1"/);(/"royalblue1","maroon1"/) 

  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJS_temp_budget"
  tdir = ddir+"/"+"anom_dtdt_r2.5.nc"
  udir = ddir+"/"+"anom_udtdx_r2.5.nc"
  vdir = ddir+"/"+"anom_vdtdy_r2.5.nc"
  wdir = ddir+"/"+"anom_wdtdp_r2.5.nc"
  t1dir = "/home/sunming/data5/cuixy/global_jets/data/ERA5_daily_anom_tvert_1979-2022_2.5x2.5.nc"
  pdir = "/home/sunming/data5/cuixy/global_jets/plot/Fig13_bar_NAJS"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
  ft  = addfile(tdir,"r")
  fu  = addfile(udir,"r")
  fv  = addfile(vdir,"r")
  fw  = addfile(wdir,"r")

  ft1 = addfile(t1dir,"r")
  t1  = ft1->t
  t1  = lonFlip(t1)

  t = ft->dtdt
  t = t
  u = fu->udtdx
  u = u*tday
  v = fv->vdtdy
  v = v*tday
  w = fw->wdtdp
  w = w*tday

  printVarSummary(t1)
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

  nt = ind(phase.eq.4.and.amp.gt.1.)
  nind = ind_DJF(nt)
;; store lead-lag values for plotting.  
;; + temp.
  llmejs = new((/41,6/),float)
  llmejs = 0.
  tm = t(:,{mlats:mlatn},{mlonw:mlone})
  clatjs = cos(tm&lat*rad)
  printVarSummary(tm)
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

  do j = 0, 40
    ntt   = -10 + j
    nindd = nind + ntt
    llmejs(j,0) = dim_avg_n(tjs(nindd), 0)
    llmejs(j,1) = dim_avg_n(ujs(nindd), 0)
    llmejs(j,2) = dim_avg_n(vjs(nindd), 0)
    llmejs(j,3) = dim_avg_n(wjs(nindd), 0)
    llmejs(j,4) = dim_avg_n(rjs(nindd), 0)
    llmejs(j,5) = dim_avg_n(t1js(nindd), 0)
    delete([/ntt,nindd/])
  end do
  delete([/clatjs,tjs,ujs,vjs,wjs,rjs/])
  print("calc MEJS data done. "+ systemfunc("date"))

;; line for dtdt
  line = new((/2,41/),float)
  line(0,:) = llmejs(:,0)
  line(1,:) = llmejs(:,5)
  print(line(0,:))
;; - temp.
  lleajs = new((/41,6/),float)
  lleajs = 0.
  te = t(:,{elats:elatn},{elonw:elone})
  clatjs = cos(te&lat*rad)

  t1js = wgt_areaave_Wrap(t1(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  tjs = wgt_areaave_Wrap(t(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  ujs = wgt_areaave_Wrap(u(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  vjs = wgt_areaave_Wrap(v(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  wjs = wgt_areaave_Wrap(w(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  rjs = wgt_areaave_Wrap(r(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  do j = 0, 40
    ntt   = -10 + j
    nindd = nind + ntt
    lleajs(j,0) = dim_avg_n(tjs(nindd), 0)
    lleajs(j,1) = dim_avg_n(ujs(nindd), 0)
    lleajs(j,2) = dim_avg_n(vjs(nindd), 0)
    lleajs(j,3) = dim_avg_n(wjs(nindd), 0)
    lleajs(j,4) = dim_avg_n(rjs(nindd), 0)
    lleajs(j,5) = dim_avg_n(t1js(nindd), 0)
    delete([/ntt,nindd/])
  end do
  delete([/tm,clatjs,tjs,ujs,vjs,wjs,rjs/])
  print("calc MEJS data done. "+ systemfunc("date"))
;; line for dtdt
  line1 = new((/2,41/),float)
  line1(0,:) = lleajs(:,0)
  line1(1,:) = lleajs(:,5)
;------------------------------------------------------
;; plot format: Line chart. about 5 lines. RGBA color.

  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@gsnRightStringFontHeightF  = 0.015
  res@gsnLeftStringFontHeightF   = 0.015
  res@gsnLeftStringOrthogonalPosF  = 0.01
  res@gsnRightStringOrthogonalPosF = 0.01

  res@tmXTOn = False
  res@tmYROn = False

  res@tiXAxisString  = "~F21~Lag days"
  res@tiYAxisString  = "~F21~Temp tendency (K/day)"
  res@tmXBMode   = "Explicit"
  res@tmXBValues = (/ 0. , 5., 10., 15., 20., 25., 30., 35., 40. /)
  res@tmXBLabels = (/"-10","-5","0","5","10","15","20","25","30" /)
  res@tmXBLabelFontHeightF = 0.012
  res@tiXAxisFontHeightF = 0.012

  res@tiYAxisFontHeightF = 0.012

  res@tmYLMode        = "Manual"
  res@tmYLMinorPerMajor = 0
  res@tmYLTickSpacingF  = 0.3
  res@tmYLLabelFontHeightF = 0.012

  res@trXMaxF = 25

  res@tmYLMajorLengthF   = -0.01
  res@trYMaxF  =  1.5
  res@trYMinF  = -1.5
  res@gsnYRefLine = 0.
  res@gsnYRefLineDashPattern = 0

  res@xyLineThicknessF = 4.0

  xlabel = ispan(0,40,1)
;----------------------------------------------------------------
;; Line plot
  wks = gsn_open_wks("x11", pdir)
  res@vpWidthF  = 0.6
  res@vpHeightF = 0.25
  res@vpXF = 0.2
  res@vpYF = 0.95
  res@trYMaxF  =  0.12
  res@trYMinF  = -0.12
  res@trXMinF  = 0
  res@trXMaxF  = 20
  res@tmYLTickSpacingF  = 0.03
  res@gsnLeftString  = "~F21~(a) (~F34~6~F21~T/~F34~6~F21~t~F21~)' & T' "
  res@xyLineColor = color(0)
  res@xyDashPattern = 1
  res@gsnRightString = ""

  resR = True
  resR@trYMaxF  =  0.8
  resR@trYMinF  = -0.8
  resR@tmYLTickSpacingF = 0.2
  resR@xyLineThicknessF = 4.0
  resR@xyLineColor = color(0)
  resR@xyDashPattern = 0
  resR@tiYAxisString  = "~F21~Temp anomaly (K)"
  dum  = gsn_csm_xy2(wks, xlabel, line(0,:),line(1,:), res,resR)

  delete([/res@gsnLeftString,res@gsnRightString/])

  res@xyLineColor  = color(1)
  resR@xyLineColor = color(1)

  dum2 = gsn_csm_xy2(wks, xlabel, line1(0,:),line1(1,:), res,resR)
  overlay(dum,dum2)

;;draw four lines.
  resl = True
  resl@gsLineColor = color(0)
  resl@gsLineDashPattern = 2
  resl@gsLineThicknessF  = 4.0  
  resl@gsLineColor = color(1)
  draw(dum)

;; text.
  rest               = True            
  rest@txFontHeightF = 0.012    

  rest@txFontColor = color(0)
  gsn_text_ndc(wks, "~F21~box1: "+mlats+"-"+mlatn+"N, "+mlonw+"-"+mlone+"E", 0.7, 0.935, rest)
  rest@txFontColor = color(1)
  gsn_text_ndc(wks, "~F21~box2: "+elats+"-"+elatn+"N, "+elonw+"-"+elone+"E", 0.702, 0.915, rest)
  draw(dum)

;; draw labels
  rest@txFontColor = "black"
  rest@txFontHeightF = 0.015

  resl = True
  resl@gsLineColor = "black"
  resl@gsLineThicknessF  = 4
  resl@gsLineDashPattern = 1

  gsn_polyline_ndc(wks, (/0.55,0.6/), (/0.96,0.96/), resl)
  gsn_text_ndc(wks,"~F21~(~F34~6~F21~T/~F34~6~F21~t~F21~)'", 0.64,0.96,rest)

  resl@gsLineDashPattern = 0
  gsn_polyline_ndc(wks, (/0.7,0.75/), (/0.96,0.96/), resl)
  gsn_text_ndc(wks,"~F21~T'", 0.77,0.96,rest)

  i = 1
  do while (i.eq.1)
    i = 1
  end do
