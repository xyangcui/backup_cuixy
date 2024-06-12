;-----------------------------------------------------------------
;; draw general budget terms. 
;; MEJS box (10-70E,10-40N)
;; EAJS box (110-180E,15-50N)
;; time: 2023-12-13
;-----------------------------------------------------------------
  rad   = 4.*atan(1.)/180. 
  tday  = 24.*60.*60.

  elatnv = 33
  elatsv = 25
  elonwv = 170
  elonev = 200

  color = (/"blue3","darkgoldenrod1"/) 

  ddir = "/home/sunming/data5/cuixy/global_jets/data/temp_budget"
  tdir = ddir+"/"+"anom_dtdt_r2.5.nc"
  udir = ddir+"/"+"anom_udtdx_r2.5.nc"
  vdir = ddir+"/"+"anom_vdtdy_r2.5.nc"
  adir = ddir+"/anom_udtdx_seperate_r2.5.nc"
  bdir = ddir+"/anom_vdtdy_seperate_r2.5.nc"
  pdir = "/home/sunming/data5/cuixy/global_jets/plot/FigS1.budget_seperate"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
  fvs  = addfile(bdir,"r")

  vltl = fvs->vltl
  vlti = fvs->vlti
  vlth = fvs->vlth
  vitl = fvs->vitl
  viti = fvs->viti
  vith = fvs->vith
  vhtl = fvs->vhtl
  vhti = fvs->vhti
  vhth = fvs->vhth

;-----------------------------------------------------------------
  fv  = addfile(vdir,"r")

  v = fv->vdtdy
  v = v*tday*(-1)
;---------------------------------------------------------------=
;; MJO phases.
  time = fv->time
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

;;;  EAJS for vdtdy.
  lseajsv = new((/10/),float)
  lseajsv = 0.
  ve = v(:,{elatsv:elatnv},{elonwv:elonev})
  clatjs = cos(ve&lat*rad)

  vjs = wgt_areaave_Wrap(v(:,{elatsv:elatnv},{elonwv:elonev})\
                         , clatjs, 1.0, 1)

;; seperate.
  vltljs = wgt_areaave_Wrap(vltl(:,{elatsv:elatnv},{elonwv:elonev})\
                         , clatjs, 1.0, 1)
  vltijs = wgt_areaave_Wrap(vlti(:,{elatsv:elatnv},{elonwv:elonev})\
                         , clatjs, 1.0, 1)
  vlthjs = wgt_areaave_Wrap(vlth(:,{elatsv:elatnv},{elonwv:elonev})\
                         , clatjs, 1.0, 1)

  vitljs = wgt_areaave_Wrap(vitl(:,{elatsv:elatnv},{elonwv:elonev})\
                         , clatjs, 1.0, 1)
  vitijs = wgt_areaave_Wrap(viti(:,{elatsv:elatnv},{elonwv:elonev})\
                         , clatjs, 1.0, 1)
  vithjs = wgt_areaave_Wrap(vith(:,{elatsv:elatnv},{elonwv:elonev})\
                         , clatjs, 1.0, 1)

  vhtljs = wgt_areaave_Wrap(vhtl(:,{elatsv:elatnv},{elonwv:elonev})\
                         , clatjs, 1.0, 1)
  vhtijs = wgt_areaave_Wrap(vhti(:,{elatsv:elatnv},{elonwv:elonev})\
                         , clatjs, 1.0, 1)
  vhthjs = wgt_areaave_Wrap(vhth(:,{elatsv:elatnv},{elonwv:elonev})\
                         , clatjs, 1.0, 1)

  nt = ind(phase.eq.1.and.amp.gt.1.)
  nind  = ind_DJF(nt)
  ntt   = -2
  nin   = nind + ntt
  do j = 1, 4
    ntt   = -2 + j
    nindd = nind + ntt
    nin  := array_append_record(nin, nindd, 0)
  end do
  lseajsv(0) = dim_avg_n(vjs(nin), 0)
  lseajsv(1) = dim_avg_n(vltljs(nin), 0)*tday*(-1)
  lseajsv(2) = dim_avg_n(vltijs(nin), 0)*tday*(-1)
  lseajsv(3) = dim_avg_n(vlthjs(nin), 0)*tday*(-1)
  lseajsv(4) = dim_avg_n(vitljs(nin), 0)*tday*(-1)
  lseajsv(5) = dim_avg_n(vitijs(nin), 0)*tday*(-1)
  lseajsv(6) = dim_avg_n(vithjs(nin), 0)*tday*(-1)
  lseajsv(7) = dim_avg_n(vhtljs(nin), 0)*tday*(-1)
  lseajsv(8) = dim_avg_n(vhtijs(nin), 0)*tday*(-1)
  lseajsv(9) = dim_avg_n(vhthjs(nin), 0)*tday*(-1)

  print("calc EAJS data done. "+ systemfunc("date"))

;------------------------------------------------------
;; plot format: Line chart. about 5 lines. RGBA color.
;; dtdt in (0,0,0,1.). udtdx in (112,160,205,1.0)
;; vdtdy in (0,52,102,1.) r in (178,178,178)
  text2  = "~F21~-(u~F34~6~F21~T~F21~/~F34~6~F21~x~F21~)'"
  text3  = "~F21~-(v~F34~6~F21~T~F21~/~F34~6~F21~y~F21~)'"
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
  resb@tmYLTickSpacingF   = 0.2
  resb@tmYLLabelFontHeightF = 0.012

  resb@tmXBMode    = "Explicit"
  resb@tmXBValues  = (/ 0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5./)
  resb@tmXBLabels  = (/text2,"~F21~LL","~F21~LI",\
                    "~F21~LH","~F21~IL","~F21~II","~F21~IH","~F21~HL","~F21~HI","~F21~HH"/)
  resb@tmXBLabelFontHeightF = 0.012

  resb@gsnXYBarChart = True
  resb@gsnXYBarChartBarWidth = 0.3

  resb@tmYLMajorLengthF   = -0.01
  resb@gsnYRefLine        = 0.
  resb@trXMinF  = 0.2
  resb@trXMaxF  = 5.3

  resb@vpWidthF  = 0.8
  resb@vpHeightF = 0.2

;----------------------------------------------------------------
;; MEJS
  wks = gsn_open_wks("x11", pdir)
;;eajs vdtdy
  resb@vpXF = 0.1
  resb@vpYF = 0.3

  resb@trYMinF  = -0.4
  resb@trYMaxF  =  0.4

  resb@gsnLeftString  = "~F21~(c) vadv"
  resb@gsnRightString = "~F21~("+elatsv+"-"+elatnv+"N"+", "+elonwv+"-"+elonev+"E"+")"
  resb@tmXBLabels  = (/text3,"~F21~LL","~F21~LI",\
                    "~F21~LH","~F21~IL","~F21~II","~F21~IH","~F21~HL","~F21~HI","~F21~HH"/)
  dumeabarv = gsn_csm_xy(wks,fspan(0.5,5.,10),lseajsv,resb)
;;polyline.
  respoly = True
  respoly@gsLineColor = "gray64"
  respoly@gsLineDashPattern = 2
  respoly@gsLineThicknessF  = 5.0
  dum3 = gsn_add_polyline(wks,dumeabarv,(/0.7,0.7/),(/10,-10/),respoly)
;; shade color for lseajs(0) and lleajs(2).
  resea = True
  resea@gsFillColor = color(1)
  dum4 = gsn_add_polygon(wks, dumeabarv, (/2.35,2.65,2.65,2.35,2.35/), (/0.,0.,lseajsv(4),lseajsv(4),0./), resea)
  draw(dumeabarv)

  i = 1
  do while (i.eq.1)
    i = 1
  end do













