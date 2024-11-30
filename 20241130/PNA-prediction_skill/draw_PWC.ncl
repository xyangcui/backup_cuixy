;----------------------------------------------------------------------
  lats = 10
  latn = 60
  lonw = 80
  lone = 240

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  zdir = ddir+"/vs/"+"ERA5_uwlevel_zonalavg_r2.5.nc"
;;jet events
  ndir = ddir+"/"+"IOnP_date_omi.csv"
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"
  idir = ddir+"/"+"IOiP_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
;-----------------------------------------------------------
; main code.
;-----------------------------------------------------------
  f = addfile(zdir,"r")
  u = f->u(:,{1000:100},{120:220})
  w = f->w(:,{1000:100},{120:220})
  e = f->e(:,{1000:100},{120:220})
  h = tofloat(f->h(:,{1000:100},{120:220}))*10^(2)
  copy_VarCoords(u,h)

  time = u&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.10.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.5.and.day.eq.31)

  mon := mon(n1:n2)
  u := u(n1:n2,:,:)
  w := w(n1:n2,:,:)
  e := e(n1:n2,:,:)
  h := h(n1:n2,:,:)
  nt = ind(mon.eq.10.or.mon.eq.11.or.mon.eq.12.or.mon.eq.1\
           .or.mon.eq.2.or.mon.eq.3.or.mon.eq.4.or.mon.eq.5)
  u := u(nt,:,:)
  w := w(nt,:,:)
  e := e(nt,:,:)
  h := h(nt,:,:)
  delete([/TIME,time,year,mon,day,nt/])
;; ONDJFMAM
  time = u&time
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

  is = isf
  ie = ief

  dims = dimsizes(w)
  dim2 = dims(1)
  dim3 = dims(2)
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

  pt = phase(nda)
  at = amp(nda)
  pt  := pt(ind(at.gt.1.))
  nda := nda(ind(at.gt.1.))
  delete(igood)
  ;select amp >1. and phase.
  i = 0
  np = 6
  print("phase "+ np +" begin.")
  ip = ind(pt.eq.np)
  ip@_FillValue = -999
  if (any(ip.eq.ip@_FillValue).eq."Missing")
    continue
  else
    out  = new((/4,dimsizes(ip),dim2,dim3/),float)
    out(0,:,:,:) = w(nda(ip)+i,:,:)
    out(1,:,:,:) = u(nda(ip)+i,:,:)
    out(2,:,:,:) = e(nda(ip)+i,:,:)
    out(3,:,:,:) = h(nda(ip)+i,:,:)

    out := dim_avg_n(out, 1) 
    print("phase "+ np +" done.")
  end if
  copy_VarCoords(u(0,:,:), out(0,:,:))

  delete([/w,u,h,e/])
  u = out(1,:,:)
  w = out(0,:,:)*10*(-1)
  h = smth9(out(3,:,:),0.5,0.25,False)
  e = smth9(out(2,:,:),0.5,0.25,False)
  copy_VarCoords(out(0,:,:),w)
  copy_VarCoords(out(0,:,:),u)
  copy_VarCoords(out(0,:,:),e)
  copy_VarCoords(out(0,:,:),h)
  print(u({850},{160:220}))

  wks = gsn_open_wks("x11", "map")
  cmap = read_colormap_file("MPL_RdBu")

  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@gsnMaximize = False

  res@gsnYAxisIrregular2Linear = True

  res@tiMainString         = "Pressure/Height Vector Example"
  res@cnFillPalette        = cmap(::-1,:)
  res@cnFillOn             = True                
  res@cnLevelSelectionMode = "ManualLevels"      
  res@cnMinLevelValF       = -1.6               
  res@cnMaxLevelValF       =  1.6               
  res@cnLevelSpacingF      =  0.2                 
  ;res@lbLabelStride        = 2     

  res@trXMaxF = 210
  res@gsnMajorLonSpacing = 10
  res@trYMinF = 200
        
  res@vcRefAnnoOrthogonalPosF  = -0.165          
  res@pmLabelBarOrthogonalPosF = -0.05           
  res@vcRefMagnitudeF          = 0.5        
  res@vcRefLengthF             = 0.01        
  res@vcMapDirection           = False

  res@vpHeightF = 0.25
  res@vpWidthF  = 0.5

  wt = w;*(-1)
  plot  = gsn_csm_pres_hgt_vector(wks, h, u, wt, res)

  draw(plot)

  i = 0 
  do while(i.eq.0)
    i = 0
  end do