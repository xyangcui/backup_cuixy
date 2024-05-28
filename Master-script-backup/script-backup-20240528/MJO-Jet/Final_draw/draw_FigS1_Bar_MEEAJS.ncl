;-----------------------------------------------------------------
;; draw general budget terms. 
;; MEJS box (10-70E,10-40N)
;; EAJS box (110-180E,15-50N)
;; time: 2023-12-13
;-----------------------------------------------------------------
  rad   = 4.*atan(1.)/180. 
  tday  = 24.*60.*60.

  mlatn = 20
  mlats = 30
  mlonw = 55
  mlone = 95

  elatn = 25
  elats = 35
  elonw = 160
  elone = 200

  color = (/"blue3","red3"/) 

  ddir = "/home/sunming/data5/cuixy/global_jets/data/temp_budget"
  tdir = ddir+"/"+"anom_dtdt_r2.5.nc"
  udir = ddir+"/"+"anom_udtdx_r2.5.nc"
  vdir = ddir+"/"+"anom_vdtdy_r2.5.nc"
  wdir = ddir+"/"+"anom_wdtdp_r2.5.nc"
  pdir  = "/home/sunming/data5/cuixy/global_jets/plot/Fig8_bar_EAJS"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
  ft  = addfile(tdir,"r")
  fu  = addfile(udir,"r")
  fv  = addfile(vdir,"r")
  fw  = addfile(wdir,"r")

  t = ft->dtdt
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
  w = -1*w

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
  nt = ind(phase.eq.2.and.amp.gt.1.)
  nind  = ind_DJF(nt)
  ntt   = -10
  nin = nind + ntt
  do j = 1, 9
    ntt   = -10 + j
    nindd = nind + ntt
    nin  := array_append_record(nin, nindd, 0)
  end do 

;; bar for other value. EAJS:27-32 day: 12-17 MEJS: day 2-7 17-22 ;;
  bar1 = new((/5,2/),float)
  lin1 = new((/10,2/),float)
  lin1 = 0.
  tm = t(:,{mlats:mlatn},{mlonw:mlone})
  clatjs = cos(tm&lat*rad)

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

  bar1(0,0) = dim_avg_n(tjs(nindd), 0)
  bar1(1,0) = dim_avg_n(ujs(nindd), 0)
  bar1(2,0) = dim_avg_n(vjs(nindd), 0)
  bar1(3,0) = dim_avg_n(wjs(nindd), 0)
  bar1(4,0) = dim_avg_n(rjs(nindd), 0)
;; seperate
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

  lin1(0,0) = dim_avg_n(ujs(nin), 0)
  lin1(1,0) = dim_avg_n(ultljs(nin), 0)*tday*(-1)
  lin1(2,0) = dim_avg_n(ultijs(nin), 0)*tday*(-1)
  lin1(3,0) = dim_avg_n(ulthjs(nin), 0)*tday*(-1)
  lin1(4,0) = dim_avg_n(uitljs(nin), 0)*tday*(-1)
  lin1(5,0) = dim_avg_n(uitijs(nin), 0)*tday*(-1)
  lin1(6,0) = dim_avg_n(uithjs(nin), 0)*tday*(-1)
  lin1(7,0) = dim_avg_n(uhtljs(nin), 0)*tday*(-1)
  lin1(8,0) = dim_avg_n(uhtijs(nin), 0)*tday*(-1)
  lin1(9,0) = dim_avg_n(uhthjs(nin), 0)*tday*(-1)
  delete([/tm,clatjs,tjs,ujs,vjs,wjs,rjs/])

;; EAJS.
  tm = t(:,{elats:elatn},{elonw:elone})
  clatjs = cos(tm&lat*rad)

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

  bar1(0,1) = dim_avg_n(tjs(nindd), 0)
  bar1(1,1) = dim_avg_n(ujs(nindd), 0)
  bar1(2,1) = dim_avg_n(vjs(nindd), 0)
  bar1(3,1) = dim_avg_n(wjs(nindd), 0)
  bar1(4,1) = dim_avg_n(rjs(nindd), 0)

;; seperate values.
  ultljs = wgt_areaave_Wrap(ultl(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  ultijs = wgt_areaave_Wrap(ulti(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  ulthjs = wgt_areaave_Wrap(ulth(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  uitljs = wgt_areaave_Wrap(uitl(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  uitijs = wgt_areaave_Wrap(uiti(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  uithjs = wgt_areaave_Wrap(uith(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  uhtljs = wgt_areaave_Wrap(uhtl(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  uhtijs = wgt_areaave_Wrap(uhti(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  uhthjs = wgt_areaave_Wrap(uhth(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)

  lin1(0,1) = dim_avg_n(ujs(nin), 0)
  lin1(1,1) = dim_avg_n(ultljs(nin), 0)*tday*(-1)
  lin1(2,1) = dim_avg_n(ultijs(nin), 0)*tday*(-1)
  lin1(3,1) = dim_avg_n(ulthjs(nin), 0)*tday*(-1)
  lin1(4,1) = dim_avg_n(uitljs(nin), 0)*tday*(-1)
  lin1(5,1) = dim_avg_n(uitijs(nin), 0)*tday*(-1)
  lin1(6,1) = dim_avg_n(uithjs(nin), 0)*tday*(-1)
  lin1(7,1) = dim_avg_n(uhtljs(nin), 0)*tday*(-1)
  lin1(8,1) = dim_avg_n(uhtijs(nin), 0)*tday*(-1)
  lin1(9,1) = dim_avg_n(uhthjs(nin), 0)*tday*(-1)
  delete([/nt,nind,ntt,nin,nindd,tm,clatjs,tjs,ujs,vjs,wjs,rjs/])
  print("calc phase2 done. "+ systemfunc("date"))


;; store lead-lag values for plotting.
  nt = ind(phase.eq.6.and.amp.gt.1.)
  nind  = ind_DJF(nt)
  ntt   = -10
  nin = nind + ntt
  do j = 1, 9
    ntt   = -10 + j
    nindd = nind + ntt
    nin  := array_append_record(nin, nindd, 0)
  end do 
;; bar for other value. EAJS:27-32 day: 12-17 MEJS: day 2-7 17-22 ;;
  bar2 = new((/5,2/),float)
  lin2 = new((/10,2/),float)
  lin2 = 0.
  tm = t(:,{mlats:mlatn},{mlonw:mlone})
  clatjs = cos(tm&lat*rad)

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

  bar2(0,0) = dim_avg_n(tjs(nindd), 0)
  bar2(1,0) = dim_avg_n(ujs(nindd), 0)
  bar2(2,0) = dim_avg_n(vjs(nindd), 0)
  bar2(3,0) = dim_avg_n(wjs(nindd), 0)
  bar2(4,0) = dim_avg_n(rjs(nindd), 0)

;; seperate
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

  lin2(0,0) = dim_avg_n(ujs(nin), 0)
  lin2(1,0) = dim_avg_n(ultljs(nin), 0)*tday*(-1)
  lin2(2,0) = dim_avg_n(ultijs(nin), 0)*tday*(-1)
  lin2(3,0) = dim_avg_n(ulthjs(nin), 0)*tday*(-1)
  lin2(4,0) = dim_avg_n(uitljs(nin), 0)*tday*(-1)
  lin2(5,0) = dim_avg_n(uitijs(nin), 0)*tday*(-1)
  lin2(6,0) = dim_avg_n(uithjs(nin), 0)*tday*(-1)
  lin2(7,0) = dim_avg_n(uhtljs(nin), 0)*tday*(-1)
  lin2(8,0) = dim_avg_n(uhtijs(nin), 0)*tday*(-1)
  lin2(9,0) = dim_avg_n(uhthjs(nin), 0)*tday*(-1)

  delete([/tm,clatjs,tjs,ujs,vjs,wjs,rjs/])

;; EAJS.
  tm = t(:,{elats:elatn},{elonw:elone})
  clatjs = cos(tm&lat*rad)

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

  bar2(0,1) = dim_avg_n(tjs(nindd), 0)
  bar2(1,1) = dim_avg_n(ujs(nindd), 0)
  bar2(2,1) = dim_avg_n(vjs(nindd), 0)
  bar2(3,1) = dim_avg_n(wjs(nindd), 0)
  bar2(4,1) = dim_avg_n(rjs(nindd), 0)
;; seperate values.
  ultljs = wgt_areaave_Wrap(ultl(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  ultijs = wgt_areaave_Wrap(ulti(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  ulthjs = wgt_areaave_Wrap(ulth(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  uitljs = wgt_areaave_Wrap(uitl(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  uitijs = wgt_areaave_Wrap(uiti(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  uithjs = wgt_areaave_Wrap(uith(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  uhtljs = wgt_areaave_Wrap(uhtl(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  uhtijs = wgt_areaave_Wrap(uhti(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)
  uhthjs = wgt_areaave_Wrap(uhth(:,{elats:elatn},{elonw:elone})\
                         , clatjs, 1.0, 1)

  lin2(0,1) = dim_avg_n(ujs(nin), 0)
  lin2(1,1) = dim_avg_n(ultljs(nin), 0)*tday*(-1)
  lin2(2,1) = dim_avg_n(ultijs(nin), 0)*tday*(-1)
  lin2(3,1) = dim_avg_n(ulthjs(nin), 0)*tday*(-1)
  lin2(4,1) = dim_avg_n(uitljs(nin), 0)*tday*(-1)
  lin2(5,1) = dim_avg_n(uitijs(nin), 0)*tday*(-1)
  lin2(6,1) = dim_avg_n(uithjs(nin), 0)*tday*(-1)
  lin2(7,1) = dim_avg_n(uhtljs(nin), 0)*tday*(-1)
  lin2(8,1) = dim_avg_n(uhtijs(nin), 0)*tday*(-1)
  lin2(9,1) = dim_avg_n(uhthjs(nin), 0)*tday*(-1)
  delete([/tm,clatjs,tjs,ujs,vjs,wjs,rjs/])

  print("calc phase6 done. "+ systemfunc("date"))
;------------------------------------------------------
;; plot format: Line chart. about 5 lines. RGBA color. 
  text1  = "~F21~(~F34~6~F21~T/~F34~6~F21~t~F21~)'"                 
  text2  = "~F21~-(u~F34~6~F21~T~F21~/~F34~6~F21~x~F21~)'"
  text3  = "~F21~-(v~F34~6~F21~T~F21~/~F34~6~F21~y~F21~)'"
  text4  = "~F21~(~F33~w~F33~s~F21~)'"
  text5  = "~F21~Q'/c~B~p~N~"
;; bar chart.
  resb = True
  resb@gsnDraw     = False
  resb@gsnFrame    = False
  resb@gsnMaximize = False

  resb@gsnLeftStringFontHeightF  = 0.015
  resb@tiYAxisString  = "~F21~Temp tendency (K/day)"
  resb@tmXTOn = False
  resb@tmYROn = False

  resb@tiYAxisFontHeightF = 0.012
  resb@tmYLMode           = "Manual"
  resb@tmYLMinorPerMajor  = 0
  resb@tmYLTickSpacingF   = 0.1
  resb@tmYLLabelFontHeightF = 0.012

  resb@tmXBMode    = "Explicit"
  resb@tmXBValues  = (/  0.3,  0.8,  1.3,  1.8,  2.3/)
  resb@tmXBLabels  = (/text1,text2,text3,text4,text5/)
  resb@tmXBLabelFontHeightF = 0.012

  resb@gsnXYBarChart = True
  resb@gsnXYBarChartBarWidth = 0.2

  resb@tmYLMajorLengthF = -0.01
  resb@gsnYRefLine      = 0.
  resb@trXMinF  =  0.0
  resb@trXMaxF  =  2.7
  resb@trYMaxF  =  1.0
  resb@trYMinF  = -1.0
  resb@tmYLTickSpacingF = 0.2
  resb@vpWidthF  = 0.4
  resb@vpHeightF = 0.25
;----------------------------------------------------------------
  wks = gsn_open_wks("x11", pdir)
  ;drawNDCGrid(wks)
;;draw Bar chart for Phase 2.
  resb@vpXF = 0.1
  resb@vpYF = 0.9
  resb@tmXBLabelsOn = True
  resb@gsnLeftString = "~F21~(a) budget terms@P2"
  resb@gsnXYBarChartColors  = color(0)
  dumbar  = gsn_csm_xy(wks,(/ 0.2, 0.7, 1.2, 1.7, 2.2/),bar1(:,0),resb)
  resb@gsnXYBarChartColors  = color(1)
  dumbar2 = gsn_csm_xy(wks,(/ 0.4, 0.9, 1.4, 1.9, 2.4/),bar1(:,1),resb)
  overlay(dumbar,dumbar2)

;;polyline.
  respoly = True
  respoly@gsLineColor = "black"
  respoly@gsLineDashPattern = 2
  respoly@gsLineThicknessF  = 4.0
  respoly@tfPolyDrawOrder   = "PostDraw"
  dum1 = gsn_add_polyline(wks,dumbar,(/0.55,0.55/),(/-10,10/),respoly)
  draw(dumbar)
  undef("resb")
;; plot format: Line chart. about 5 lines. RGBA color.
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
  resb@tmXBValues  = (/ 1.5,4.5,7.5,10.5,13.5,15.5,18.5,21.5,24.5/)
  resb@tmXBLabels  = (/"~F21~LL","~F21~LI","~F21~LH",\
                       "~F21~IL","~F21~II","~F21~IH",\
                       "~F21~HL","~F21~HI","~F21~HH"/)
  resb@tmXBLabelFontHeightF = 0.012

  resb@gsnXYBarChart = True
  resb@gsnXYBarChartBarWidth = 1.0

  resb@tmYLMajorLengthF   = -0.01
  resb@gsnYRefLine        = 0.
  resb@trXMinF  = 0
  resb@trXMaxF  = 26

  resb@vpWidthF  = 0.4
  resb@vpHeightF = 0.25
  resb@vpXF = 0.1
  resb@vpYF = 0.55
;----------------------------------------------------------------
  resb@gsnLeftStringFontHeightF    = 0.016
  resb@gsnRightStringFontHeightF   = 0.016            
;;draw Bar chart.
  resb@trYMaxF  =  0.8
  resb@trYMinF  = -0.7
  resb@gsnLeftString  = "~F21~(c) "+text2+"@P2"
  resb@gsnRightString = ""
  resb@gsnXYBarChartColors  = color(0)
  dummebar  = gsn_csm_xy(wks,fspan(1,24,9),lin1(1::,0),resb)
  resb@gsnXYBarChartColors  = color(1)
  dummebar2 = gsn_csm_xy(wks,fspan(2,25,9),lin1(1::,1),resb)
  overlay(dummebar,dummebar2)
  draw(dummebar)

;;draw Bar chart for Phase 6.
  undef("resb")
;; bar chart.
  resb = True
  resb@gsnDraw     = False
  resb@gsnFrame    = False
  resb@gsnMaximize = False

  resb@gsnLeftStringFontHeightF  = 0.015
  resb@tiYAxisString  = ""
  resb@tmXTOn = False
  resb@tmYROn = False

  resb@tiYAxisFontHeightF = 0.012
  resb@tmYLMode           = "Manual"
  resb@tmYLMinorPerMajor  = 0
  resb@tmYLTickSpacingF   = 0.1
  resb@tmYLLabelFontHeightF = 0.012

  resb@tmXBMode    = "Explicit"
  resb@tmXBValues  = (/  0.3,  0.8,  1.3,  1.8,  2.3/)
  resb@tmXBLabels  = (/text1,text2,text3,text4,text5/)
  resb@tmXBLabelFontHeightF = 0.012

  resb@gsnXYBarChart = True
  resb@gsnXYBarChartBarWidth = 0.2

  resb@tmYLMajorLengthF = -0.01
  resb@gsnYRefLine      = 0.
  resb@trXMinF  =  0.0
  resb@trXMaxF  =  2.7
  resb@trYMaxF  =  0.6
  resb@trYMinF  = -0.6
  resb@tmYLTickSpacingF = 0.15
  resb@vpWidthF  = 0.4
  resb@vpHeightF = 0.25
;;------------------------------------------------------
  resb@vpXF = 0.57
  resb@vpYF = 0.9
  resb@tmXBLabelsOn = True
  resb@gsnLeftString = "~F21~(b) budget terms@P6"
  resb@gsnXYBarChartColors  = color(0)
  dumbart  = gsn_csm_xy(wks,(/ 0.2, 0.7, 1.2, 1.7, 2.2/),bar2(:,0),resb)
  resb@gsnXYBarChartColors  = color(1)
  dumbar2t = gsn_csm_xy(wks,(/ 0.4, 0.9, 1.4, 1.9, 2.4/),bar2(:,1),resb)
  overlay(dumbart,dumbar2t)

;;polyline.
  respoly = True
  respoly@gsLineColor = "black"
  respoly@gsLineDashPattern = 2
  respoly@gsLineThicknessF  = 4.0
  respoly@tfPolyDrawOrder   = "PostDraw"
  dum1 = gsn_add_polyline(wks,dumbart,(/0.55,0.55/),(/-10,10/),respoly)
  draw(dumbart)
  undef("resb")
;; plot format: Line chart. about 5 lines. RGBA color.
  text2  = "~F21~-(u~F34~6~F21~T~F21~/~F34~6~F21~x~F21~)'"
;; bar chart.
  resb = True
  resb@gsnDraw     = False
  resb@gsnFrame    = False
  resb@gsnMaximize = False

  resb@tmXTOn = False
  resb@tmYROn = False

  resb@tiYAxisFontHeightF = 0.012
  resb@tiYAxisString      = ""
  resb@tmYLMode           = "Manual"
  resb@tmYLMinorPerMajor  = 0
  resb@tmYLTickSpacingF   = 0.3
  resb@tmYLLabelFontHeightF = 0.012

  resb@tmXBMode    = "Explicit"
  resb@tmXBValues  = (/ 1.5,4.5,7.5,10.5,13.5,15.5,18.5,21.5,24.5/)
  resb@tmXBLabels  = (/"~F21~LL","~F21~LI","~F21~LH",\
                       "~F21~IL","~F21~II","~F21~IH",\
                       "~F21~HL","~F21~HI","~F21~HH"/)
  resb@tmXBLabelFontHeightF = 0.012

  resb@gsnXYBarChart = True
  resb@gsnXYBarChartBarWidth = 1.0

  resb@tmYLMajorLengthF   = -0.01
  resb@gsnYRefLine        = 0.
  resb@trXMinF  = 0
  resb@trXMaxF  = 26

  resb@vpWidthF  = 0.4
  resb@vpHeightF = 0.25
  resb@vpXF = 0.57
  resb@vpYF = 0.55
;----------------------------------------------------------------
  resb@gsnLeftStringFontHeightF    = 0.016
  resb@gsnRightStringFontHeightF   = 0.016            
;;draw Bar chart.
  resb@trYMaxF  =  0.8
  resb@trYMinF  = -0.7
  resb@gsnLeftString  = "~F21~(d) "+text2+"@P6"
  resb@gsnRightString = ""
  resb@gsnXYBarChartColors  = color(0)
  dummebart  = gsn_csm_xy(wks,fspan(1,24,9),lin2(1::,0),resb)
  resb@gsnXYBarChartColors  = color(1)
  dummebar2t = gsn_csm_xy(wks,fspan(2,25,9),lin2(1::,1),resb)
  overlay(dummebart,dummebar2t)
  draw(dummebart)

  i = 1
  do while (i.eq.1)
    i = 1
  end do