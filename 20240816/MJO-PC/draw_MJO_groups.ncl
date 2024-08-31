
;--------------------------------------------------------------------------
;function: 返回位相3对应的日期 格式：YYYYMMDD
;--------------------------------------------------------------------------
function return_p3(timt:integer,xdir,pht:integer,amt:float)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,pt,at,ip
begin
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nf = asciiread(xdir,-1,"string")
  is = str_get_field(nf,1,delim)
  ie = str_get_field(nf,2,delim)
  is := toint(is(1::))
  ie := toint(ie(1::))
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ;event all days.
  nda := ndtooned(nda)
  nda := nda(ind(.not.ismissing(nda)))  ;index for z.
  pt = pht(nda)
  at = amt(nda)
  pt  := pt(ind(at.gt.1.))
  nda := nda(ind(at.gt.1.))
  ;select amp >1. and phase.
  m = 3
  ip = ind(pt.eq.3)
  ip@_FillValue = -999
  if (any(ip.eq.ip@_FillValue).ne."Missing")
    nda := nda(ip)
  end if
  return nda

end

  lats = -15
  latn =  15
  lonw =  10
  lone =  180

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  odir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
;;jet events
  ndir = ddir+"/"+"IOnP_date_omi.csv"
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"
  idir = ddir+"/"+"IOiP_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
;-----------------------------------------------------------------
; main code.
;-----------------------------------------------------------------
  f = addfile(odir,"r")
  olr = f->olr(:,{lats:latn},{lonw:lone})
  time = olr&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  timi = cd_calendar(time, -2) 
;;OMI
  delim = ","
  rmmf  = asciiread(omidir,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
  phaset = str_get_field(rmmf,4,delim)
  ampt = str_get_field(rmmf,7,delim)

  iyear := toint(iyear(1::))
  imon  := toint(imon(1::))
  iday  := toint(iday(1::))
  phaset := tointeger(phaset(1::))
  ampt := tofloat(ampt(1::))

  phase = new((/dimsizes(time)/),integer)
  amp   = new((/dimsizes(time)/),float)
  phase = 0
  amp   = 0.

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i)\
             .and.iday.eq.day(i))
    phase(i) = phaset(nt)
    amp(i) = ampt(nt)
    delete(nt)
  end do
  delete([/phaset,ampt,iyear,imon,iday/])
  print("complete OMI selling "+ systemfunc("date"))

;;MJO events.
;;jet events (year,mon,day) (events,days,lat,lon)
  print("jet events begin. "+ systemfunc("date"))
  IOnP = return_p3(timi,ndir,phase,amp) ;IOnP
  IOiP = return_p3(timi,idir,phase,amp) ;IOiP
  IOfP = return_p3(timi,fdir,phase,amp) ;IOfP
  IOsP = return_p3(timi,sdir,phase,amp) ;IOsP
  print("jet events done. "+ systemfunc("date"))

  olrt = dim_avg_n_Wrap(olr, 1)
  hov  = new((/4,46,dimsizes(olr&lon)/),float)

  do i = 0, 45
    nt1 = IOnP-15+i
    hov(0,i,:) = dim_avg_n(olrt(nt1,:),0)
    nt2 = IOiP-15+i
    hov(1,i,:) = dim_avg_n(olrt(nt2,:),0)
    nt3 = IOfP-15+i
    hov(2,i,:) = dim_avg_n(olrt(nt3,:),0)
    nt4 = IOsP-15+i
    hov(3,i,:) = dim_avg_n(olrt(nt4,:),0)
  end do
  copy_VarCoords(olrt(0,:), hov(0,0,:))
  hov!1 = "time"
  hov&time = ispan(-15,30,1)
;--------------------------------------------------------------------
  wks  = gsn_open_wks ("x11", "hov" )  

  cmap = read_colormap_file("MPL_RdBu")         

  res             = True  
  res@gsnDraw     = False
  res@gsnMaximize = False
  res@gsnFrame    = False
  res@gsnRightStringFontHeightF = 0.015
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringOrthogonalPosF = 0.0
  res@gsnLeftStringOrthogonalPosF  = 0.0

  res@tmYLLabelFontHeightF = 0.012
  res@tmXBLabelFontHeightF = 0.012

  res@tmYROn = False
  res@tmXTOn = False

  ;res@tiXAxisString = "~F21~RMM phase"
  res@tiYAxisString = "~F21~Day"
  res@tiYAxisFontHeightF = 0.015
  res@tiXAxisFontHeightF = 0.015

  ;res@tmXBMode = "Explicit" 
  ;res@tmXBValues = z&lat
  ;res@tmXBLabels = (/"8","1","2","3","4","5","6","7","8" /)

  res@cnFillOn       = True              
  res@cnFillPalette  = "BlueWhiteOrangeRed";cmap(::-1,:)  
  res@cnLinesOn = False    
  res@cnLevelSelectionMode = "ManualLevels"     
  res@cnMinLevelValF       = -25         
  res@cnMaxLevelValF       =  25             
  res@cnLevelSpacingF      =  5

  res@lbOrientation  = "Vertical" 
  res@lbLabelStride  = 2
     
  res@vpWidthF  = 0.25
  res@vpHeightF = 0.35

  label = (/"(a) IOnP","(b) IOiP","(c) IOfP","(d) IOsP"/)

  res@vpXF = 0.18
  res@vpYF = 0.95
  res@gsnLeftString = "~F21~"+label(0)
  plot1 = gsn_csm_hov(wks, hov(0,:,:), res) 

  res@vpXF = 0.6
  res@vpYF = 0.95
  res@gsnLeftString = "~F21~"+label(1)
  plot2 = gsn_csm_hov(wks, hov(1,:,:), res) 

  res@vpXF = 0.18
  res@vpYF = 0.47
  res@gsnLeftString = "~F21~"+label(2)
  plot3 = gsn_csm_hov(wks, hov(2,:,:), res) 

  res@vpXF = 0.6
  res@vpYF = 0.47
  res@gsnLeftString = "~F21~"+label(3)
  plot4 = gsn_csm_hov(wks, hov(3,:,:), res) 
  
  draw(plot1)
  draw(plot2)
  draw(plot3)
  draw(plot4)

  i = 0
  do while(i.eq.0)
    i = 0
  end do
