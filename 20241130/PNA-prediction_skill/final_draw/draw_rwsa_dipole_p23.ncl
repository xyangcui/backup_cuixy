load "/home/sunming/data5/cuixy/script/ncl_func/calc_rws.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"

undef("ErrorBar")
procedure ErrorBar(wks1,plotebar,xx,yy,yup,ydown,corr)
local wks1,xx,yy,yup,ydown,up1,down1,polyres,polyres2,\
dum1,plotebar
begin
  mstring =       "-"
  fontnum =       19
  xoffset =       0.0
  yoffset =       0.0
  aspect  =       1E-10
  size    =       1.0
  angle   =       0.0

  new_index = NhlNewMarker(wks1, mstring, fontnum, xoffset, yoffset,\
   aspect, size, angle)

  COLOR    = corr
  LThick   = 3
  MarkSiz  = 0.015
  Dpattern = 0
  polyres                       = True
  polyres@gsLineColor           = COLOR
  polyres@gsLineDashPattern     = Dpattern
  polyres@gsLineThicknessF      = LThick
  polyres@tfPolyDrawOrder       = "PostDraw"

  res_lines = True
  res_lines@gsLineDashPattern = Dpattern                   
  res_lines@gsLineThicknessF  = LThick                    
  res_lines@gsLineColor       = COLOR     

  dum1   = gsn_add_polyline(wks1,plotebar,(/ydown,yup/),(/xx,xx/),polyres)
  up1    = gsn_add_polyline(wks1,plotebar,(/yup,yup/),(/xx-0.065,xx+0.065/),res_lines)
  down1  = gsn_add_polyline(wks1,plotebar,(/ydown,ydown/),(/xx-0.065,xx+0.065/),res_lines)

  str1 = unique_string("dum")
  str2 = unique_string("up")
  str3 = unique_string("down")
  plotebar@$str1$ = dum1
  plotebar@$str2$ = up1
  plotebar@$str3$ = down1

end


function stats_frequency(timt,pha,ampt,is,ie,x,y)
local is,ie,x,y
begin
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ;event all days.
  nda  := ndtooned(nda)
  igood = ind(.not.ismissing(nda))
  nda  := nda(igood)  ;index for z.
  pt = pha(nda)
  at = ampt(nda)
  ;pt  := pt(ind(at.gt.1.))
  ;nda := nda(ind(at.gt.1.))
  delete(igood)
  ;select amp >1. and phase.
  ip = ind(pt.eq.2)
  ip@_FillValue = -999
  if (any(ip.eq.ip@_FillValue).eq."Missing")
    continue
  else
  out  = new(dimsizes(ip),float)
  do k = 0,dimsizes(ip)-1
    r1 = x(nda(ip(k)))
    r2 = y(nda(ip(k)))
    if (r1.lt.0.and.r2.gt.0) then
      out(k) = 1
    else
      out(k) = 0
    end if
  end do
  out := dim_sum(out)
  rws_f = 100.*dim_sum(out)/dimsizes(ip)
  end if 
  return rws_f
end
;----------------------------------------------------------------------
  lats = 10
  latn = 60
  lonw = 80
  lone = 165

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  zdir = ddir+"/"+"rws.nc"
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
  z = f->rws(:,{lats:latn},{lonw:lone},0)
  ;z := dim_sum_n_Wrap(z, 3)
  printVarSummary(z)
  time = z&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.10.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.5.and.day.eq.31)

  mon := mon(n1:n2)
  z  := z(n1:n2,:,:)
  nt = ind(mon.eq.10.or.mon.eq.11.or.mon.eq.12.or.mon.eq.1\
           .or.mon.eq.2.or.mon.eq.3.or.mon.eq.4.or.mon.eq.5)
  z   := z(nt,:,:)
  delete([/TIME,time,year,mon,day,nt/])
;; ONDJFMAM
  time = z&time
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
  delete([/phaset,ampt,iyear,imon,iday,TIME,year,mon,day/])
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
  ;pt  := pt(ind(at.gt.1.))
  ;nda := nda(ind(at.gt.1.))
  delete(igood)
  ;select amp >1. and phase.
  ip = ind(pt.eq.3.or.pt.eq.2)
  ip@_FillValue = -999
  if (any(ip.eq.ip@_FillValue).eq."Missing")
    continue
  else
    zt := dim_avg_n(z(nda(ip),:,:),0)
    zs = z(nda(ip),:,:)*10^(10)
  end if
  rwss := zt*10^(10)
  copy_VarCoords(z(0,:,:), rwss)
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
  ;pt  := pt(ind(at.gt.1.))
  ;nda := nda(ind(at.gt.1.))
  delete(igood)
  ;select amp >1. and phase.
  ip = ind(pt.eq.3.or.pt.eq.2)
  ip@_FillValue = -999
  if (any(ip.eq.ip@_FillValue).eq."Missing")
    continue
  else
    zt  := dim_avg_n(z(nda(ip),:,:),0)
    zf = z(nda(ip),:,:)*10^(10)
  end if
  rwsf := zt*10^(10)
  copy_VarCoords(z(0,:,:), rwsf)
  delete([/zt,ip,nda,pt,at,nda,nint,is,ie/])

  rwsd = rwsf-rwss
  copy_VarCoords(z(0,:,:), rwsd)
; rws averaged.
  latne = 37;30
  latse = 27;37
  lonwe = 120;105
  lonee = 135;120

  rad  = 4.0*atan(1.0)/180.0
  z1 = z(:,{latse:latne},{lonwe:lonee})
  lat1 = z1&lat
  clat1 = cos(lat1*rad)
  rws1  = wgt_areaave_Wrap(z1, clat1, 1.0, 1)

  latnw = 25;25
  latsw = 35;35
  lonww = 140;130
  lonew = 155;145

  z2 = z(:,{latsw:latnw},{lonww:lonew})
  lat2 = z2&lat
  clat2 = cos(lat2*rad)
  rws2  = wgt_areaave_Wrap(z2, clat2, 1.0, 1)  

;rws1>0 rws2<0.
  mt = 1
  rff = stats_frequency(timi,phase,amp,isf,ief,rws1,rws2)
  rfs = stats_frequency(timi,phase,amp,iss,ies,rws1,rws2)

  print(rff)
  print(rfs)

  odir = ddir+"/"+"boostrap_rws_omi_p23_region1.nc"
  fout = addfile(odir,"r")           
  rs   = fout->rf
  rs := rs*100
  nrand = dimsizes(rs(0,:))
  rs_me = rs
  do i = 0,nrand - 1
    rs_me(0,i) = rs(0,i) - rff
    rs_me(1,i) = rs(1,i) - rfs
  end do
  rs_sorted = dim_pqsort_n(rs_me, 2, 1)

  ind_up = toint(nrand*0.95-1)
  ind_dn = toint(nrand*0.05-1)

  sdf_up = rs_me(0,ind_up)
  sdf_dn = rs_me(0,ind_dn)

  sds_up = rs_me(1,ind_up)
  sds_dn = rs_me(1,ind_dn)

  print(sdf_up)
  print(sdf_dn)
  xf = (/1.2,1.7/)
  wks = gsn_open_wks("x11", "bar")
  fcolor = "dodgerblue1";chocolate1"
  scolor = "firebrick1";"deepskyblue1"

  res             = True  
  res@gsnDraw     = False
  res@gsnMaximize = False
  res@gsnFrame    = False
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnLeftStringOrthogonalPosF = -0.15


  res@gsnRightString = "~F21~Frequency [%]  "
  res@gsnRightStringFontHeightF  = 0.015
  res@gsnRightStringOrthogonalPosF = -0.95
  ;res@tiXAxisFontHeightF = 0.012
  ;res@tiYAxisFontHeightF = 0.015
  ;res@tiXAxisOffsetYF = 0.008

  res@tmYROn = False
  res@tmXTOn = False
  res@tmXBLabelFontHeightF = 0.012
  res@tmYLLabelFontHeightF = 0.015
  res@tmYLMode = "Explicit" 
  res@tmYLValues = xf
  res@tmYLLabels = (/"~F21~IOfP","~F21~IOsP"/)

  res@vpXF = 0.58
  res@vpYF = 0.86+0.065
  res@vpWidthF  = 0.35
  res@vpHeightF = 0.18

  res@trYMaxF = 2.35
  res@trYMinF = 0.85
  res@trXMinF = 10.
  res@trXMaxF = 50.
  res@gsnXRefLine  = 10.

  res@gsnXYBarChart = True
  res@gsnXYBarChartBarWidth = 0.28
  res@gsnXYBarChartColors = (/fcolor,scolor/)

  res@gsnLeftString = "~F21~ (b) Freq. (EA~B~RWS~N~<0&WP~B~RWS~N~>0)"  
  plot1 = gsn_csm_xy(wks,(/rff,rfs/),xf,res)  

  ErrorBar(wks,plot1,xf(0),rff,rff-sdf_dn,rff-sdf_up,"dodgerblue4")
  ErrorBar(wks,plot1,xf(1),rfs,rfs-sds_dn,rfs-sds_up,"firebrick4")

;;底图 RWSa
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
  resu@gsnLeftStringOrthogonalPosF  = -0.15

  resu@tmXTOn = False
  resu@tmYROn = False
  resu@tmXBLabelFont = 21
  resu@tmYLLabelFont = 21
  resu@tmXBLabelFontHeightF = 0.012
  resu@tmYLLabelFontHeightF = 0.012

  resu@mpMinLatF = 15
  resu@mpMaxLatF = 60
  resu@mpMinLonF = 100
  resu@mpMaxLonF = 245
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
  resu@lbTitleOn = True
  resu@lbTitleString = "[10~S~-10~N~ s~S~-2~N~]"
  resu@lbTitlePosition  = "Right"
  resu@lbTitleDirection = "Across"
  resu@lbTitleFontHeightF = 0.015
  resu@lbBoxLinesOn       = True
  resu@lbLabelStride = 2
  resu@lbLabelFontHeightF = 0.012
  resu@pmLabelBarOrthogonalPosF = 0.15

  resu@vpXF = 0.15
  resu@vpYF = 0.95
  ;res@vpWidthF  = 0.3
  ;res@vpHeightF = 0.185
  resu@gsnLeftString = "~F21~ (a) RWS (IOfP - IOsP)" 
  plot2 = gsn_csm_contour_map(wks, rwsd, resu)

  respoly = True
  respoly@gsLineColor      = "black"
  respoly@gsLineThicknessF = 3.0

  txres               = True 
  txres@txFontHeightF = 0.015
  txres@txFont = 22
;; EAJS
  lonmin = lonww
  lonmax = lonew
  latmin = latsw
  latmax = latnw
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum1 = gsn_add_polyline(wks, plot2, polyx, polyy, respoly)
  tex1 = gsn_add_text(wks, plot2, "WP~B~RWS~N~", 150, 22, txres)
;; MEJS
  lonmin = lonwe
  lonmax = lonee
  latmin = latse
  latmax = latne
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum2 = gsn_add_polyline(wks, plot2, polyx, polyy, respoly)
  txt2 = gsn_add_text(wks, plot2, "EA~B~RWS~N~", 128,24, txres)

  draw(plot2)
  draw(plot1)
;--------------------------------------------------------------------------
;function: 计算每种事件并返回矩阵的程序。
;--------------------------------------------------------------------------
function calc_pro(timt:integer,xdir,zt:float,pht:integer,ampt:float)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,igood,pt,at,ip,out,pct1
begin
  dims = dimsizes(zt)
  dim2 = dims(1)
  dim3 = dims(2)
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nf = asciiread(xdir,-1,"string")
  is = str_get_field(nf,1,delim)
  ie = str_get_field(nf,2,delim)
  is := toint(is(1::))
  ie := toint(ie(1::))
  ind_time = ispan(0,dimsizes(timt)-1,1)
  ;event by event.
  pct1 = new((/dimsizes(is),dim2,dim3/),float) 
  pct1 = 0.
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    pt = pht(its:ite)
    at = ampt(its:ite)
    nda = ind_time(its:ite)
    ;pt  := pt(ind(at.gt.1.))
    ;nda := nda(ind(at.gt.1.))
    ip = ind(pt.eq.3.or.pt.eq.2)
    ip@_FillValue = -999
    if (any(ip.eq.ip@_FillValue).eq."Missing")
      delete([/ip,its,ite,at,pt,nda/])
      continue
    else
      if (dimsizes(ip).eq.1)
        pct1(i,:,:) = zt(nda(ip),:,:)
      else
        pct1(i,:,:) = dim_avg_n(zt(nda(ip),:,:),0)
      end if
    end if
    delete([/ip,its,ite,at,pt,nda/])
  end do
  return pct1

end
;------------------------------------------------------------------

  lats = -25
  latn =  25
  lonw =  50
  lone =  230

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  odir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
;;jet events
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
;-----------------------------------------------------------------
; main code.
;-----------------------------------------------------------------
  f = addfile(odir,"r")
  olr = f->olr(:,{lats:latn},{lonw:lone})
  undef("time")
  time = olr&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  undef("timi")
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

  undef("phase")
  undef("amp")
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
  IOsP = calc_pro(timi,sdir,olr,phase,amp) ;IOsP
  IOfP = calc_pro(timi,fdir,olr,phase,amp) ;IOfP
  print("jet events done. "+ systemfunc("date"))

  IOsP!0 = "time"
  IOsP!1 = "lat"
  IOsP!2 = "lon"

  IOfP!0 = "time"
  IOfP!1 = "lat"
  IOfP!2 = "lon"
;;value for plot.
  If = dim_avg_n(IOfP,0)
  Is = dim_avg_n(IOsP,0)
  Id = Is - If

  copy_VarCoords(olr(0,:,:),If)
  copy_VarCoords(olr(0,:,:),Is)
  copy_VarCoords(olr(0,:,:),Id)
;;for ttest.
  undef("nf")
  undef("ns")
  nf = dimsizes(IOsP(:,0,0))
  ns = dimsizes(IOfP(:,0,0))

;(t-test)
  doff = equiv_sample_size(IOfP(lat|:,lon|:,time|:), 0.1, 0)
  dff  =  toint(wgt_areaave(doff, 1., 1., 0))

  dofs = equiv_sample_size(IOsP(lat|:,lon|:,time|:), 0.1, 0)
  dfs  =  toint(wgt_areaave(dofs, 1., 1., 0))
;;sst
  f_var=dim_variance_n(IOfP, 0)
  f_std=sqrt(f_var/nf)
  f_std=where(f_std.eq.0,f_std@_FillValue,f_std);

  t_f=If/f_std
  confi_f=If
  confi_f=student_t(t_f, nf-1)

  s_var=dim_variance_n(IOsP, 0)
  s_std=sqrt(s_var/ns)
  s_std=where(s_std.eq.0,s_std@_FillValue,s_std)

  t_s=Is/s_std
  confi_s=Is
  confi_s=student_t(t_s, ns-1)

  copy_VarCoords(If, confi_f)
  copy_VarCoords(If, confi_s)
;;diff
  confi_d = ttest(If,f_var,nf-1, Is,s_var,ns-1, False, False)
  alpha = 100.*(1. - confi_d) 

  Is = smth9_Wrap(Is, 0.5, 0.25, False)
  If = smth9_Wrap(If, 0.5, 0.25, False)
  Id = smth9_Wrap(Id, 0.5, 0.25, False)
;-------------------------------------------------------------------- 
  undef("cmap")
  cmap = read_colormap_file("BlueWhiteOrangeRed")         
;shading u plot
  undef("res")
  res = True
  res@gsnFrame           = False
  res@gsnDraw            = False
  res@gsnMaximize        = False
  res@gsnAddCyclic       = False
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnLeftStringOrthogonalPosF  = -0.15

  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBLabelFont = 21
  res@tmYLLabelFont = 21
  res@tmXBLabelFontHeightF = 0.012
  res@tmYLLabelFontHeightF = 0.012

  res@mpMinLatF = lats
  res@mpMaxLatF = latn
  res@mpMinLonF = lonw
  res@mpMaxLonF = lone
  res@mpCenterLonF = 140
  res@mpFillOn        = False
  res@mpLandFillColor = "white"
  res@mpGridAndLimbDrawOrder = "PostDraw"
  res@mpGridLineColor        = "Gray60"

  res@cnFillOn               = True
  res@cnFillPalette          = cmap
  res@cnLinesOn              = False
  res@cnLineLabelsOn         = False
  res@cnLevelSelectionMode   = "Manuallevels"
  res@cnMinLevelValF         = -24
  res@cnMaxLevelValF         =  24
  res@cnLevelSpacingF        =  4.
  res@cnLabelBarEndStyle     = "ExcludeOuterBoxes"

  res@lbTitleOn = True
  res@lbTitleString = "[W m~S~-2~N~]"
  res@lbTitlePosition  = "Right"
  res@lbTitleDirection = "Across"
  res@lbTitleFontHeightF = 0.015
  res@lbBoxLinesOn       = True
  res@lbLabelFontHeightF = 0.012
  res@pmLabelBarOrthogonalPosF = 0.15

  undef("opt")
  opt                      =True
  opt@gsnShadeFillType     ="pattern"
  opt@gsnShadeLow          =17;打点
  opt@gsnAddCyclic         =True
  opt@cnFillDotSizeF       =0.001
  opt@cnFillColor          ="gray64"

  undef("sres")
  sres                     =True
  sres@gsnDraw             =False
  sres@gsnFrame            =False
  sres@cnLinesOn           =False
  sres@gsnLeftString       =""
  sres@gsnRightString      =""
  sres@cnLevelSelectionMode="ExplicitLevels"
  sres@cnLevels            =(/0.15,0.1,0.05,0.01/)
  sres@cnFillPalette       ="GMT_gray"
  sres@cnFillColor         ="gray64"
  sres@cnLineLabelsOn      =False
  sres@cnInfoLabelOn       =False
  sres@lbLabelBarOn        =False
  sres@cnFillScaleF        = 2.0
  sres@cnFillDotSizeF      = 0.005

  n = 0
  res@vpXF = 0.15
  res@vpYF = 0.59
  res@lbLabelBarOn  = False
  res@gsnLeftString = "~F21~ (c) OLR (IOsP P2&3)"
  plots = gsn_csm_contour_map(wks, Is, res) 
  bases=gsn_csm_contour(wks,confi_s, sres)
  bases=gsn_contour_shade(bases,0.05,1,opt)
  overlay(plots,bases)

  n = 1
  res@vpXF = 0.15
  res@vpYF = 0.325
  res@lbLabelBarOn  = True
  res@gsnLeftString = "~F21~ (d) OLR (IOfP P2&3)"
  plotf = gsn_csm_contour_map(wks, If, res) 
  basef=gsn_csm_contour(wks,confi_f, sres)
  basef=gsn_contour_shade(basef,0.05,1,opt)
  overlay(plotf,basef)

  resc = True
  resc@gsnFrame     = False
  resc@gsnDraw      = False
  resc@gsnMaximize  = False
  resc@gsnAddCyclic = False
  resc@gsnContourNegLineDashPattern = 0
  resc@gsnLeftString      = ""
  resc@gsnRightString     = ""
  resc@cnMonoLineColor    = True
  resc@cnLinesOn          = True
  resc@cnLineLabelsOn     = False
  resc@cnLineColor        = "blue";magenta
  resc@cnLineThicknessF   = 2.5
  resc@cnLineDashPattern  = 0
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels = (/-2.,-4.,-6.,-8./)

  ;plotc = gsn_csm_contour(wks, Id, resc)
  ;overlay(plot,plotc)
  ;draw(plot)

  draw(plots)
  draw(plotf)

  i = 0 
  do while(i.eq.0)
    i = 0
  end do