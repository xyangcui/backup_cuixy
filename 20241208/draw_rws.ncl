load "/home/sunming/data5/cuixy/script/ncl_func/calc_rws.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
;----------------------------------------------------------------------
  lats = 10
  latn = 60
  lonw = 80
  lone = 240

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  zdir = ddir+"/"+"rws.nc"

  umdir = ddir+"/ERA5_u200_daily_mean_1979-2022_r2.5.nc"
  vmdir = ddir+"/ERA5_v200_daily_mean_1979-2022_r2.5.nc"

  rmdir = ddir+"/ERA5_rwsv_daily_mean_1979-2022_r2.5.nc"
  ridir = ddir+"/ERA5_rwsv_daily_anom_1979-2022_r2.5.nc"
;;jet events
  ndir = ddir+"/"+"IOnP_date_omi.csv"
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"
  idir = ddir+"/"+"IOiP_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
;-----------------------------------------------------------
; main code.
;-----------------------------------------------------------
;calc u and v mean.
  fm = addfile(rmdir,"r")
  fi = addfile(ridir,"r")
  vom = fm->vor({lats:latn},{lonw:lone})
  um  = fm->u({lats:latn},{lonw:lone})
  vm  = fm->v({lats:latn},{lonw:lone})
;read in ui and calc rws.
  voi = fi->vor(:,{lats:latn},{lonw:lone})
  ui = fi->u(:,{lats:latn},{lonw:lone})
  vi = fi->v(:,{lats:latn},{lonw:lone})
;;times
  f = addfile(zdir,"r")
  z = f->rws(:,{lats:latn},{lonw:lone},0)
  printVarSummary(z)
  time = z&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.10.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.5.and.day.eq.31)

  mon := mon(n1:n2)
  voi := voi(n1:n2,:,:)
  ui := ui(n1:n2,:,:)
  vi := vi(n1:n2,:,:)
  z  := z(n1:n2,:,:)
  nt = ind(mon.eq.10.or.mon.eq.11.or.mon.eq.12.or.mon.eq.1\
           .or.mon.eq.2.or.mon.eq.3.or.mon.eq.4.or.mon.eq.5)
  voi := voi(nt,:,:)
  ui  := ui(nt,:,:)
  vi  := vi(nt,:,:)
  z   := z(nt,:,:)
  delete([/TIME,time,year,mon,day,nt/])
;; ONDJFMAM
  time = ui&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  timi = cd_calendar(time,-2)
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
;----------------------------------------------------------------------
;;jet events (year,mon,day) (events,days,lat,lon)
  nf  = asciiread(idir,-1,"string")
  is  = str_get_field(nf,1,delim)
  ie  = str_get_field(nf,2,delim)
  ist = toint(is(1::))
  iet = toint(ie(1::))

  nf := asciiread(ndir,-1,"string")
  is := str_get_field(nf,1,delim)
  ie := str_get_field(nf,2,delim)

  ist := array_append_record(ist, toint(is(1::)), 0)
  iet := array_append_record(iet, toint(ie(1::)), 0)

  nf := asciiread(fdir,-1,"string")
  is := str_get_field(nf,1,delim)
  ie := str_get_field(nf,2,delim)

  isf = toint(is(1::))
  ief = toint(ie(1::))

  ist := array_append_record(ist, toint(is(1::)), 0)
  iet := array_append_record(iet, toint(ie(1::)), 0)

  nf := asciiread(sdir,-1,"string")
  is := str_get_field(nf,1,delim)
  ie := str_get_field(nf,2,delim) 

  iss = toint(is(1::))
  ies = toint(ie(1::))

  ist := array_append_record(ist, toint(is(1::)), 0)
  iet := array_append_record(iet, toint(ie(1::)), 0)

  delete([/ie,is/])

  is = ist
  ie = iet

  delim = ","
  nint = ispan(0,dimsizes(timi),1) ; time index.
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timi.eq.is(i))
    ite = ind(timi.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ;event all days.
  nda  := ndtooned(nda)
  igood = ind(.not.ismissing(nda))
  nda  := nda(igood)  ;index for z.
  ;nda = nint
  pt = phase(nda)
  at = amp(nda)
  pt  := pt(ind(at.gt.1.))
  nda := nda(ind(at.gt.1.))
  delete(igood)
  ;select amp >1. and phase.
  print("phase "+ np +" begin.")
  ip = ind(pt.eq.6.or.pt.eq.7)
  ip@_FillValue = -999
  if (any(ip.eq.ip@_FillValue).eq."Missing")
    continue
  else
    zt  := dim_avg_n(z(nda(ip),:,:),0)
    print("phase "+ np +" done.")
  end if
  rws := zt*10^(10)
  copy_VarCoords(ui(0,:,:), rws)
  delete([/zt,ip,nda,pt,at,nda,nint,is,ie/])
;ss
  is = iss
  ie = ies

  delim = ","
  nint = ispan(0,dimsizes(timi),1) ; time index.
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timi.eq.is(i))
    ite = ind(timi.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ;event all days.
  nda  := ndtooned(nda)
  igood = ind(.not.ismissing(nda))
  nda  := nda(igood)  ;index for z.
  ;nda = nint
  pt = phase(nda)
  at = amp(nda)
  pt  := pt(ind(at.gt.1.))
  nda := nda(ind(at.gt.1.))
  delete(igood)
  ;select amp >1. and phase.
  print("phase "+ np +" begin.")
  ip = ind(pt.eq.6.or.pt.eq.7)
  ip@_FillValue = -999
  if (any(ip.eq.ip@_FillValue).eq."Missing")
    continue
  else
    zt  := dim_avg_n(z(nda(ip),:,:),0)
    print("phase "+ np +" done.")
  end if
  rwss := zt*10^(10)
  copy_VarCoords(ui(0,:,:), rwss)
  delete([/zt,ip,nda,pt,at,nda,nint,is,ie/])

; ff 
  is = isf
  ie = ief

  delim = ","
  nint = ispan(0,dimsizes(timi),1) ; time index.
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timi.eq.is(i))
    ite = ind(timi.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ;event all days.
  nda  := ndtooned(nda)
  igood = ind(.not.ismissing(nda))
  nda  := nda(igood)  ;index for z.
  ;nda = nint
  pt = phase(nda)
  at = amp(nda)
  pt  := pt(ind(at.gt.1.))
  nda := nda(ind(at.gt.1.))
  delete(igood)
  ;select amp >1. and phase.
  print("phase "+ np +" begin.")
  ip = ind(pt.eq.6.or.pt.eq.7)
  ip@_FillValue = -999
  if (any(ip.eq.ip@_FillValue).eq."Missing")
    continue
  else
    zt  := dim_avg_n(z(nda(ip),:,:),0)
    print("phase "+ np +" done.")
  end if
  rwsf := zt*10^(10)
  copy_VarCoords(ui(0,:,:), rwsf)
  delete([/zt,ip,nda,pt,at,nda,nint,is,ie/])

  rwsd = rwss-rwsf
  copy_VarCoords(ui(0,:,:), rwsd)

  wks = gsn_open_wks("x11", "map")
  cmap = read_colormap_file("MPL_RdBu")

  resu = True
  resu@gsnFrame           = False
  resu@gsnDraw            = False
  resu@gsnMaximize        = False
  resu@gsnAddCyclic       = False
  resu@gsnMajorLatSpacing = 10
  resu@gsnMajorLonSpacing = 15
  resu@gsnRightStringFontHeightF = 0.015 
  resu@gsnLeftStringFontHeightF  = 0.015

  resu@gsnRightStringOrthogonalPosF = 0.01
  resu@gsnLeftStringOrthogonalPosF  = 0.01

  resu@tmXTOn = False
  resu@tmYROn = False
  resu@tmXBLabelFont = 21
  resu@tmYLLabelFont = 21
  resu@tmXBLabelFontHeightF = 0.01
  resu@tmYLLabelFontHeightF = 0.01

  resu@mpMinLatF = 15
  resu@mpMaxLatF = 55
  resu@mpMinLonF = 80
  resu@mpMaxLonF = 180
  resu@mpFillOn        = False
  resu@mpLandFillColor = "white"
  resu@mpGridAndLimbDrawOrder = "PostDraw"
  resu@mpGridLineColor        = "Gray60"

  resu@cnFillOn               = True
  resu@cnFillPalette          = cmap(::-1,:)
  resu@cnLinesOn              = False
  resu@cnLineLabelsOn         = False
  resu@cnLevelSelectionMode   = "Manuallevels"
  resu@cnMinLevelValF         = -2.0
  resu@cnMaxLevelValF         =  2.0
  resu@cnLevelSpacingF        =  0.2
  resu@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  resu@lbLabelBarOn           = True

  resw = True
  resw@gsnDraw  = False
  resw@gsnFrame = False
  resw@vcRefMagnitudeF = 0.3
  resw@vcRefLengthF    = 0.01

  plot  = gsn_csm_contour_map(wks, rwsd, resu)
  ;plot2 = gsn_csm_vector(wks, rwsd(::2,::2), rwsd(::2,::2), resw)
  ;overlay(plot, plot2)

  respoly = True
  respoly@gsLineColor      = "black"
  respoly@gsLineThicknessF = 3.0

  txres               = True 
  txres@txFontHeightF = 0.012
  txres@txFont = 22
;; EAJS
  lonmin = 135
  lonmax = 150
  latmin = 25
  latmax = 35
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum1 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)
  tex1 = gsn_add_text(wks, plot, "WNP", 140, 37, txres)
;; MEJS
  lonmin = 105
  lonmax = 120
  latmin = 30
  latmax = 40
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum2 = gsn_add_polyline(wks, plot, polyx, polyy, respoly)
  txt2 = gsn_add_text(wks, plot, "EIO", 110,43, txres)

  draw(plot)

  i = 0 
  do while(i.eq.0)
    i = 0
  end do