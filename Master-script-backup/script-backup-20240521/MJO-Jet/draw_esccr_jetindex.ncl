;------------------------------------------------
;; draw esccr. 3 jet intensity index.
;; mxlag = 9
;; author: cui xiangyang  time: 2024-04-26
;------------------------------------------------
;-----------------------------------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/global_jets/data"
  udir = ddir+"/ERA5_daily_anom_u250_1979-2022.nc"
  tdir = ddir+"/"+"ERA5_daily_anom_tvert_1979-2022.nc"
;; jet indices.
  eadir = ddir+"/jet_indices_EAJS.nc"
  medir = ddir+"/jet_indices_MEJS.nc"
  nadir = ddir+"/jet_indices_NAJS_1.nc"

  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
  pdir = "/home/sunming/data5/cuixy/global_jets/plot/Fig2_uclim_jet_indices1"

  nwa = 100.
  nwb = 20.
;-----------------------------------------------------------------------------
; Function to filter.
;-----------------------------------------------------------------------------
undef("ISO_Filter")
function ISO_Filter(x,time,npt,nwa,nwb)
local time
begin
  TIME1 = cd_calendar(time, 0)
  year1 = toint(TIME1(:,0))
  mon1  = toint(TIME1(:,1))
  day1  = toint(TIME1(:,2))
  ddd   = day_of_year(year1, mon1, day1)
  yyyyddd = year1*1000 + ddd
  wts     = filwgts_lanczos (npt, 2, 1./nwa,1./nwb , 1.0)
  xtmp  = new((/dimsizes(x),2,2/),float)
  xtmp  = 0.0
  xtmp(:,0,0) = x
  xClim = clmDayTLL(xtmp, yyyyddd)
  xAnom = calcDayAnomTLL(xtmp, yyyyddd, xClim)
  xr    = xAnom(:,0,0)
  xr := wgt_runave_Wrap(xr, wts, 0)

  return(xr)
  delete([/TIME1,year1,mon1,day1,xtmp,xAnom,xr/])

end
;-----------------------------------------------------------------------------
; Main code.
;-----------------------------------------------------------------------------
;read in data and selected strong MJO days.
;----------------------------------------------------------------
  fea = addfile(eadir,"r")
  fme = addfile(medir,"r")
  fna = addfile(nadir,"r")

  time = fea->time
  TIME = cd_calendar(time, 0)
  year = toint( TIME(:,0) )
  mon  = toint( TIME(:,1) )
  day  = toint( TIME(:,2) )

  RMMf  = rmmdir
  delim = ","
  rmmf  = asciiread(RMMf,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
  phase = str_get_field(rmmf,6,delim)
  amp   = str_get_field(rmmf,7,delim)

  iyear := toint(iyear(1::))
  imon  := toint(imon(1::))
  iday  := toint(iday(1::))
  phase := toint(phase(1::))
  amp   := tofloat(amp(1::))

  reamp   = new((/dimsizes(time)/),float)
  rephase = new((/dimsizes(time)/),integer)
  reamp   = 0.
  rephase = 0

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
    reamp(i)   = amp(nt)
    rephase(i) = phase(nt)
    delete(nt)
  end do
;--------------------------------------------------------------
;; jet index.
  var1 = "jets"
  var2 = "jetl"
;; NAJS.
  jetsna1 = fna->$var1$ ; jet speed
  jetsna = ISO_Filter(jetsna1,time,201,nwa,nwb)
;; EAJS.
  jetsea1 = fea->$var1$ ; jet speed
  jetsea = ISO_Filter(jetsea1,time,201,nwa,nwb)
;; EAJS.
  jetsme1 = fme->$var1$ ; jet speed
  jetsme = ISO_Filter(jetsme1,time,201,nwa,nwb)

;-------------------------------------------------------
;; temp.
  mlatn = 31
  mlats = 23
  mlonw = 60
  mlone = 100

  elatn = 35
  elats = 25
  elonw = 160
  elone = 200

  ft = addfile(tdir,"r")

  tm  = ft->t(:,{mlats:mlatn},{mlonw:mlone})
  te  = ft->t(:,{elats:elatn},{elonw:elone})

  rad = atan(1.0)*4/180.
  clatm = cos(tm&latitude*rad)
  tm := wgt_areaave_Wrap(tm, clatm, 1.0, 1)
  clate = cos(te&latitude*rad)
  te := wgt_areaave_Wrap(te, clate, 1.0, 1)
  delete([/clatm,clate/])
;-------------------------------------------------------
;; zonal wind.
  mlatn = 33
  mlats = 28
  mlonw = 70
  mlone = 100

  elatn = 40
  elats = 30
  elonw = 150
  elone = 190
  
  fu = addfile(udir,"r")

  um  = fu->u(:,{mlats:mlatn},{mlonw:mlone})
  ue  = fu->u(:,{elats:elatn},{elonw:elone})

  rad = atan(1.0)*4/180.
  clatm = cos(um&latitude*rad)
  um := wgt_areaave_Wrap(um, clatm, 1.0, 1)
  clate = cos(ue&latitude*rad)
  ue := wgt_areaave_Wrap(ue, clate, 1.0, 1)

;first conver from 12 to 3.
  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.31)
  reamp   := reamp(n1:n2)
  rephase := rephase(n1:n2)

  jetsna := jetsna(n1:n2)
  jetsea := jetsea(n1:n2)
  jetsme := jetsme(n1:n2)

  tm := tm(n1:n2)
  te := te(n1:n2)

  um := um(n1:n2)
  ue := ue(n1:n2)

  mon  := mon(n1:n2)
  time := time(n1:n2)

  nd   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  reamp  := reamp(nd)
  in1  = ispan(0,dimsizes(time)-1,1)
  in2  = in1(nd)

  nt  = ind(reamp.gt.1.)
  nt1 = in2(nt)
;; calc esccr.
  sme = tm 
  sea = te
  mxlag = 20
  me_lead_ea = esccr(sme(nd), sea(nd), mxlag)
  ea_lead_me = esccr(sea(nd), sme(nd), mxlag)

  ccr = new( 2*mxlag+1, float)
  ccr(0:mxlag-1) = me_lead_ea(1:mxlag:-1)
  ccr(mxlag:) = ea_lead_me(0:mxlag)

  print(ccr)

;;plot 
  dayt = ispan(-mxlag, mxlag, 1)

  wks = gsn_open_wks("x11",pdir)
  rts           = True
  rts@gsnDraw   = False       
  rts@gsnFrame  = False       
  rts@gsnScale  = True  

  rts@trYMinF = -1.0
  rts@trYMaxF =  1.0             

  rts@xyLineThicknessF = 2
  rts@gsnYRefLine  = 0.                 
  rts@gsnXRefLine  = 0.                 
  rts@gsnCenterString  = "r(MEJS,EAJS)"

  plot = gsn_csm_xy (wks,dayt,ccr ,rts)

  draw(plot)

  i = 0
  do while (i.eq.0)
    i = 0
  end do