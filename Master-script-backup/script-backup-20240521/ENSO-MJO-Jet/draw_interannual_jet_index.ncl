;-----------------------------------------------------------
;-- MEJS: (45-75E,25-35N)
;-- EAJS: (170-210E, 25-40N)
;-- NAJS: (-60 - -25W, 30-55N) (300-335)
;-----------------------------------------------------------
  yrStrt = 19791201
  yrLast = 20220201

  nlats = 30
  nlatn = 55
  nlonw = 300
  nlone = 335

  mlats = 25
  mlatn = 35
  mlonw = 45
  mlone = 75

  elats = 25
  elatn = 40
  elonw = 170
  elone = 210

  pi = atan(1.0)*4.

  ddir = "/home/sunming/data5/cuixy/ENSO-MJO-Jet/data"

  f  = addfile (ddir+"/ERA5_u200_mon_mean_1979-2022_r1.nc", "r")

  YYYY   = cd_calendar(f->time,-2)   
  iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast) 
  ;u  = f->uwnd(iYYYY,{200},{lats:latn},{lonw:lone})

  ;clat = cos(u&lat*pi/180.)
  ;ut = wgt_areaave_Wrap(u, clat, 1.0, 1)

  u  = f->uwnd(iYYYY,{200},:,:)

  ut = u(:,{mlats:mlatn},{mlonw:mlone})
  clat = cos(ut&lat*pi/180.)
  um = wgt_areaave(ut, clat, 1.0, 1)
  delete([/ut,clat/])

  ut = u(:,{elats:elatn},{elonw:elone})
  clat = cos(ut&lat*pi/180.)
  ue = wgt_areaave(ut, clat, 1.0, 1)
  delete([/ut,clat/])

  ut = u(:,{nlats:nlatn},{nlonw:nlone})
  clat = cos(ut&lat*pi/180.)
  un = wgt_areaave(ut, clat, 1.0, 1)
  delete([/ut,clat/])

  ymd  = cd_calendar(u&time, 0)
  year = toint(ymd(:,0))
  mon  = toint(ymd(:,1))

  nt = ind(mon.eq.12)
  year(nt) = year(nt)+1
  ;ua = new((/43/),float)
  ;do i = 1980,2022
  ;  ni = ind(year.eq.i.and.(mon.eq.2.or.mon.eq.1.or.mon.eq.12))
  ;  ua(i-1980) = dim_avg(ut(ni))
  ;  delete(ni)
  ;end do

  ua = new((/3,43/),float)
  do i = 1980,2022
    ni = ind(year.eq.i.and.(mon.eq.2.or.mon.eq.1.or.mon.eq.12))
    ua(0,i-1980) = dim_avg(um(ni))
    ua(1,i-1980) = dim_avg(ue(ni))
    ua(2,i-1980) = dim_avg(un(ni))
    delete(ni)
  end do

  ua := dim_standardize_n(ua, 0,1)
  ;ua := runave_n(ua, 5, 0,1)

  ninodir = "/home/sunming/data5/cuixy/data/clim_indices/Nino_3.4_index.csv"
  pdir = "/home/sunming/data5/cuixy/ENSO-MJO-Jet/figure/bar_nino_index"

  delim  = ","
  ninof  = asciiread(ninodir,-1,"string")
  iyear  = str_get_field(ninof,1,delim)
  iJan   = str_get_field(ninof,2,delim)
  iFeb   = str_get_field(ninof,3,delim)
  iDec   = str_get_field(ninof,4,delim)

  iyear:= toint(iyear(1::))
  iJan := tofloat(iJan(1::))
  iFeb := tofloat(iFeb(1::))
  iDec := tofloat(iDec(1::))

  nino = new((/43/),float)
  nino = 0.
  do i = 0,42
    nino(i) = (iJan(i)+iFeb(i)+iDec(i))/3.
  end do

  COR = esccr(nino, ua(0,:), 0)
  print(COR)

  ;cf = ezfftf_n(ua, 0)
  ;cf(:,0:2) = 0.
  ;printVarSummary(cf)
  ;uc = ezfftb_n(cf, cf@xbar, 1)

;--plot               
  res   = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False

  res@vpWidthF  = 0.7
  res@vpHeightF = res@vpWidthF*0.618

  res@tiYAxisString = ""     
  res@tiXAxisString = ""
  res@tiMainString  = "" 

  res@trYMinF = -1.0
  res@trYMaxF =  1.0

  res@xyLineColors = (/"red","blue","green"/)
  res@xyDashPatterns = (/0,0,0/)
  res@xyLineThicknesses = (/3.,3.,3./)
  res@gsnYRefLine   = 0.   

  res1 = res
  res1@xyDashPattern = 1 
              
  wks   = gsn_open_wks ("x11", "jet_index" )  

  ;plot = gsn_csm_xy2(wks, ispan(1980, 2022, 1), ua, uc, res,res1)

  plot = gsn_csm_xy(wks,ispan(1980, 2022, 1),ua,res)

  resp = True
  resp@gsLineColor = "black"
  resp@gsLineThicknessF = 2.5

  dum1 = gsn_add_polyline(wks, plot, (/1900,2200/), (/0.7,0.7/), resp)
  dum2 = gsn_add_polyline(wks, plot, (/1900,2200/), (/-0.7,-0.7/), resp)

  draw(plot)

  i = 0
  do while (i.eq.0)
    i = 0
  end do
