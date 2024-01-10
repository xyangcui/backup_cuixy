load "/home/sunming/data5/cuixy/script/ncl_func/calc_spectral.ncl"
;--------------------------------------------------------
;; spectral analysis.
;--------------------------------------------------------
pi = atan(1.0)*4.0
undef("deseason")
function deseason(x,time)
local time
begin
  TIME1 = cd_calendar(time, 0)
  year1 = toint(TIME1(:,0))
  mon1  = toint(TIME1(:,1))
  day1  = toint(TIME1(:,2))
  ddd   = day_of_year(year1, mon1, day1)
  yyyyddd = year1*1000 + ddd

  xtmp  = new((/dimsizes(x),2,2/),float)
  xtmp  = 0.0
  xtmp(:,0,0) = x
  xClim = clmDayTLL(xtmp, yyyyddd)
  xAnom = calcDayAnomTLL(xtmp, yyyyddd, xClim)
  xr    = xAnom(:,0,0)

  return(xr)
  delete([/TIME1,year1,mon1,day1,xtmp,xAnom,xr/])

end

;-----------------------------------------------------------------
;; mean_spectral.
;-----------------------------------------------------------------
undef("spec")
function spec(y)
begin
  
  nday = 90
  nsea = 43
  d    = 0
  sm   = 1
  pct  = 0.10

  Time = y&time
  TIME = cd_calendar(Time,0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  nt = ind(mon.eq.1.or.mon.eq.2)
  year(nt) = year(nt) - 1
  delete(nt)

  spcavg = new(nday/2,typeof(y))
  spcavg = 0.0
  xtsVar = 0.0
  r1zsum = 0.0

  do n = 0, 42
    iyear = n + 1979
    nt    = ind(year.eq.iyear)
    yz    = dim_standardize(y(nt),0)
    dof   = specx_anal(yz,d,sm,pct)
    spcavg = spcavg + dof@spcx
    r1     = dof@xlag1
    xtsVar = xtsVar + dof@xvari
    r1zsum = r1zsum + 0.5*log((1.+r1)/(1.-r1))
    delete([/nt,yz/])
  end do

  r1z = r1zsum/nsea
  r1  = (exp(2.*r1z)-1.)/(exp(2.*r1z)+1.)
  xvari   = xtsVar/nsea
  spcavg  = spcavg/nsea                   
  
  df      = nsea*2.0
  df@spcx = spcavg     
  df@frq  = dof@frq
  df@xlag1= r1   
  df@bw   = dof@bw     
  df@xvari= xvari
  splt = specx_ci(df, 0.90, 0.95)  
  splt@frq = dof@frq

  return(splt) 

end 

;--------------------------------------------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/global_jets/data/jet_indices_"+reg+".nc"
  udir = "/home/sunming/data5/cuixy/global_jets/data/ERA5_daily_u250_1979-2022_all.nc"
  pdir = "/home/sunming/data5/cuixy/global_jets/plot/"+reg+"spec"

;; uwind
  fu = addfile(udir,"r")
  u  = fu->u(:,0,{lats:latn},{lonw:lone})
  clat = cos(u&latitude*pi/180.)
  uc   = wgt_areaave_Wrap(u, clat, 1.0, 1)
  uc   := deseason(uc,u&time)
  delete(u)
;; Jet indices.
  f  = addfile(ddir,"r")
  jets1 = f->jets
  jetl1 = f->jetl
  jets1 := deseason(jets1,f->time)
  jetl1 := deseason(jetl1,f->time)

  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.10.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.4.and.day.eq.1)

  year  := year(n1:n2)
  mon   := mon(n1:n2)
  jets1 := jets1(n1:n2)
  jetl1 := jetl1(n1:n2)
  uc    := uc(n1:n2)

  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  jets  = jets1(nt)
  jetl  = jetl1(nt)
  u     = uc(nt)
;---------------------------------------------------------------------------------------
  xspec1  = spec(u)
  xspec   = xspec1*conform(xspec1,xspec1@frq,1)
  x1spec1 = spec(jets)
  x1spec  = x1spec1*conform(x1spec1,x1spec1@frq,1)
  x2spec1 = spec(jetl)
  x2spec  = x2spec1*conform(x2spec1,x2spec1@frq,1)

  wks = gsn_open_wks("x11",pdir)
  res = True
  res@gsnFrame  = False
  res@gsnDraw   = False
  res@gsnMaximize  = False
  res@tiYAxisString       = "Power" 
  
  res@trXReverse = True
  res@trXLog  = True 
  res@trXMinF = 2      
  
  res@xyLineThicknesses   =  (/3.0,2.5,2.5,2.5/) 
  res@xyDashPatterns      = (/0,0,1,1/)
  res@xyLineColors        = (/"black","red4","blue4","green4"/) 

  res@tmXBMode   = "Explicit"
  res@tmXBValues = (/ 140, 100, 60, 40,30, 20, 10, 5 /)
  res@tmXBLabels = (/"140","100","60","40","30","20","10","5" /)       

  fr  = xspec1@frq
  p   = 1/fr                                    
  p!0 = "f"
  p&f = fr
  p@long_name = "period"
  p@units     = "day"
  res@tiXAxisString = "Period (days)"                   
  ip   = ind(p.le.150)                       

  res@gsnCenterString = "Jets"
  xc = x1spec
  plot1 = gsn_csm_xy(wks,p(ip),xc(:,ip),res)
  res@gsnCenterString = "Jetl"
  xc = x2spec
  plot2 = gsn_csm_xy(wks,p(ip),xc(:,ip),res)
  res@gsnCenterString = "u"
  xc = xspec
  plot3 = gsn_csm_xy(wks,p(ip),xc(:,ip),res)

;; draw addition.
  respoly = True
  respoly@gsLineDashPattern = 2
  respoly@gsLineColor = "gray"
  respoly@gsLineThicknessF = 2.0
  dum1 = gsn_add_polyline(wks, plot1, (/100,100/), (/-5,10000/), respoly)
  dum2 = gsn_add_polyline(wks, plot1, (/60,60/), (/-5,10000/), respoly)
  dum4 = gsn_add_polyline(wks, plot1, (/30,30/), (/-5,10000/), respoly)
  dum5 = gsn_add_polyline(wks, plot1, (/20,20/), (/-5,10000/), respoly)
  dum6 = gsn_add_polyline(wks, plot1, (/10,10/), (/-5,10000/), respoly)

  dum70 = gsn_add_polyline(wks, plot2, (/100,100/), (/-5,10000/), respoly)
  dum80 = gsn_add_polyline(wks, plot2, (/60,60/), (/-5,10000/), respoly)
  dum40 = gsn_add_polyline(wks, plot2, (/30,30/), (/-5,10000/), respoly)
  dum50 = gsn_add_polyline(wks, plot2, (/20,20/), (/-5,10000/), respoly)
  dum60 = gsn_add_polyline(wks, plot2, (/10,10/), (/-5,10000/), respoly)

  dum100 = gsn_add_polyline(wks, plot3, (/100,100/), (/-5,10000/), respoly)
  dum200 = gsn_add_polyline(wks, plot3, (/60,60/), (/-5,10000/), respoly)
  dum400 = gsn_add_polyline(wks, plot3, (/30,30/), (/-5,10000/), respoly)
  dum500 = gsn_add_polyline(wks, plot3, (/20,20/), (/-5,10000/), respoly)
  dum600 = gsn_add_polyline(wks, plot3, (/10,10/), (/-5,10000/), respoly)

  resP = True
  resP@gsnMaximize      = False
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = False
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02
  resP@gsnPanelRowSpec  = True
  gsn_panel(wks, (/plot3,plot1,plot2/), 3, resP)

  i = 1
  do while (i.eq.1)
    i = 1
  end do
