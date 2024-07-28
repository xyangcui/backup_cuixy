load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
;-----------------------------------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/global_jets/data"
  udir = ddir+"/ERA5_daily_u250_1979-2022_r2.5.nc"
  stddir  = ddir+"/ERA5_u250_10-120std.nc"
  plotdir = "/home/sunming/data5/cuixy/global_jets/plot/Fig1_spec"

  mlats = 25
  mlatn = 35
  mlonw = 10;10
  mlone = 70;60

  elats = 30
  elatn = 40
  elonw = 150;150;110
  elone = 210;210;180

  nlats =  30
  nlatn =  55
  nlonw =  270;-90
  nlone =  330;0
;-----------------------------------------------------------------------------
; Main code.
;-----------------------------------------------------------------------------
  fs   = addfile(stddir,"r")
  ustd = fs->u(0,:,:)

  f    = addfile(udir,"r")
  ut = f->u(:,0,:,:)
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)

  mon := mon(n1:n2)
  ut  := ut(n1:n2,:,:)
  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  u = dim_avg_n_Wrap(ut(nt,:,:),0)
  delete([/n1,n2,nt/])
; u250 resources.
;--plot
  res1 = True
  res1@gsnFrame     = False
  res1@gsnDraw      = False
  res1@gsnMaximize  = False
  res1@gsnAddCyclic = True
  res1@gsnMajorLatSpacing = 15
  res1@gsnMajorLonSpacing = 30

  res1@gsnLeftStringOrthogonalPosF  = 0.012
  res1@gsnRightStringOrthogonalPosF = 0.012
  res1@gsnLeftStringFontHeightF  = 0.012
  res1@gsnRightStringFontHeightF = 0.012
  res1@tmXBLabelFontHeightF  = 0.012
  res1@tmYLLabelFontHeightF  = 0.012

  res1@vpXF     = 0.15
  res1@vpYF     = 0.9
  res1@vpWidthF = 0.7
  res1@vpHeightF= 0.4

  res1@mpFillOn = True
  res1@mpLandFillColor = "white"
  res1@mpCenterLonF      = 180.
  res1@mpMinLatF         = 0

  cmap1  = read_colormap_file("MPL_jet")
  cmap2 = cmap1(49::,:)
  res1@cnFillOn = True
  res1@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
  res1@cnFillPalette        = "WhiteBlueGreenYellowRed";cmap2
  res1@cnLinesOn            = False
  res1@cnLineLabelsOn       = False
  res1@cnLevelSelectionMode = "Manuallevels"
  ;res1@cnMinLevelValF       =  5.
  ;res1@cnMaxLevelValF       =  70.
  ;res1@cnLevelSpacingF      =  1.

  res1@lbBoxLinesOn          = True
  res1@lbBoxSeparatorLinesOn = False
  ;res1@lbLabelStride         = 10
  res1@lbTitleString         = "~F21~(unit: m s~S~-1~N~)"
  res1@lbTitlePosition       = "Bottom"
  res1@lbTitleFontHeightF    = 0.01
  res1@lbOrientation         = "Vertical"
  res1@pmLabelBarHeightF     = 0.2
  res1@pmLabelBarWidthF      = 0.1
  ;res1@pmLabelBarParallelPosF   = -0.003
  res1@pmLabelBarOrthogonalPosF = -0.01

  res1@tmXBLabelFont = 21
  res1@tmYRLabelFont = 21

  res1@gsnLeftString    = "~F21~(a) ~F21~250hPa Filtered_U_std"
  res1@gsnRightString   = "~F21~DJF"

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

  wks = gsn_open_wks("x11",plotdir)
  plot  = new(1,graphic)
  plotu = plot
  plot  = gsn_csm_contour_map(wks, ustd, res1)
  plotu = gsn_csm_contour(wks,u,resc)
  overlay(plot,plotu)

  respoly = True
  respoly@gsLineColor      = "black"
  respoly@gsLineThicknessF = 3.0
;; EAJS
  lonmin = elonw
  lonmax = elone
  latmin = elats
  latmax = elatn
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum1 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)
;; MEJS
  lonmin = mlonw
  lonmax = mlone
  latmin = mlats
  latmax = mlatn
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum2 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)
;; NAJS
  lonmin = nlonw
  lonmax = nlone
  latmin = nlats
  latmax = nlatn
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum3 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)

  draw(plot)
  delete([/udir,u,time,TIME,mon,cmap1,cmap2/])

;;----------------------------------------------------------------------
;; draw spec.
;-----------------------------------------------------------------------
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

;--------------------------------------------------------------------------------------
  udir = "/home/sunming/data5/cuixy/global_jets/data/ERA5_daily_u250_1979-2022_all.nc"
;; uwind
  fu = addfile(udir,"r")
  u  = fu->u(:,0,{mlats:mlatn},{mlonw:mlone})
  clat = cos(u&latitude*pi/180.)
  ume  = wgt_areaave_Wrap(u, clat, 1.0, 1)
  ume := deseason(ume,u&time)
  delete([/u,clat/])

  u  = fu->u(:,0,{elats:elatn},{elonw:elone})
  clat = cos(u&latitude*pi/180.)
  uea  = wgt_areaave_Wrap(u, clat, 1.0, 1)
  uea := deseason(uea,u&time)
  delete([/u,clat/])

  u  = fu->u(:,0,{nlats:nlatn},{nlonw:nlone})
  clat = cos(u&latitude*pi/180.)
  una  = wgt_areaave_Wrap(u, clat, 1.0, 1)
  una := deseason(una,u&time)
  delete([/u,clat/])

;; Jet indices.
  time = fu->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)

  mon   := mon(n1:n2)
  ume  := ume(n1:n2)
  uea  := uea(n1:n2)
  una  := una(n1:n2)

  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  ume  := ume(nt)
  uea  := uea(nt)
  una  := una(nt)

  ume  := dim_standardize(ume,0) 
  uea  := dim_standardize(uea,0) 
  una  := dim_standardize(una,0) 

;---------------------------------------------------------------------------------------
  d   = 0
  sm  = 41
  tp  = 0.10

  x     = specx_anal(ume, d, sm, tp)
  xspec1= specx_ci(x,0.5,0.9)
  xspec = xspec1*conform(xspec1,x@frq,1)
  x1      = specx_anal(uea, d, sm, tp)
  x1spec1  = specx_ci(x1,0.5,0.9)
  x1spec = x1spec1*conform(x1spec1,x1@frq,1)
  x2      = specx_anal(una, d, sm, tp)
  x2spec1  = specx_ci(x2,0.5,0.9)
  x2spec = x2spec1*conform(x2spec1,x2@frq,1)

  res = True
  res@gsnFrame  = False
  res@gsnDraw   = False
  res@gsnMaximize  = False
  res@gsnLeftStringFontHeightF  = 0.012
  res@gsnRightStringFontHeightF = 0.012
  res@tiYAxisString = "Power" 
  res@tiXAxisString = "Period (days)" 

  res@trXReverse = True
  res@trXLog     = True 
  res@trXMinF    = 5     

  res@xyLineThicknesses   =  (/3.0,2.5,2.5,2.5/) 
  res@xyDashPatterns      = (/0,0,1,1/)
  res@xyLineColors        = (/"black","Transparent","Transparent","green4"/) 

  res@tmXBMode   = "Explicit"
  res@tmXBValues = (/ 120, 100, 60, 40,30, 20, 10, 5 /)
  res@tmXBLabels = (/"120","100","60","40","30","20","10","5" /)       

  i = 20
  fr  = x@frq
  p   = 1/fr                                    
  p!0 = "f"
  p&f = fr
  p@long_name = "period"
  p@units     = "day"
  p := p(::i)            
  ip = ind(p.le.150)                       

  res@gsnLeftString  = "~F21~(b) MEJS"
  res@gsnRightString = "~F21~"+mlats+"-"+mlatn+"N,"+mlonw+"-"+mlone+"E"
  xc = xspec(:,::i)

  res@vpWidthF  =  0.25
  res@vpHeightF =  0.25
  res@vpXF = 0.1
  res@vpYF = 0.53
  res@tiYAxisString = "~F21~Power"
  plot1 = gsn_csm_xy(wks,p(ip),xc(:,ip),res)

  res@tiYAxisString = ""
  res@gsnLeftString = "~F21~(c) EAJS"
  res@gsnRightString = "~F21~"+elats+"-"+elatn+"N,"+elonw+"-"+elone+"E"
  xc = x1spec(:,::i)
  res@vpXF = 0.4
  plot2 = gsn_csm_xy(wks,p(ip),xc(:,ip),res)

  res@gsnLeftString = "~F21~(d) NAJS"
  res@gsnRightString = "~F21~"+nlats+"-"+nlatn+"N,"+nlonw+"-"+nlone+"E"
  xc = x2spec(:,::i)
  res@vpXF = 0.7

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

  draw(plot1)
  draw(plot3)
  draw(plot2)

  i = 0
  do while(i.eq.0)
    i = 0
  end do
