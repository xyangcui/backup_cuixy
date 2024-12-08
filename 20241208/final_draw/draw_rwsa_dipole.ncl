
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

function stats_frequency67(timt,pha,ampt,is,ie,x,y)
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
  ip = ind(pt.eq.6)
  ip@_FillValue = -999
  if (any(ip.eq.ip@_FillValue).eq."Missing")
    continue
  else
  out  = new(dimsizes(ip),float)
  do k = 0,dimsizes(ip)-1
    r1 = x(nda(ip(k)))
    r2 = y(nda(ip(k)))
    if (r1.gt.0.and.r2.lt.0) then
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
  printVarSummary(z)
  ;z := dim_sum_n_Wrap(z, 3)
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
  ip = ind(pt.eq.6)
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
  pt = phase(nda)
  at = amp(nda)
  ;pt  := pt(ind(at.gt.1.))
  ;nda := nda(ind(at.gt.1.))
  delete(igood)
  ;select amp >1. and phase.
  ip = ind(pt.eq.6)
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

  rwsd = rwss-rwsf
  copy_VarCoords(z(0,:,:), rwsd)

  zf!0 = "time"
  zf!1 = "lat"
  zf!2 = "lon"

  zs!0 = "time"
  zs!1 = "lat"
  zs!2 = "lon"

; rws averaged.
  latne = 30;27
  latse = 40;37
  lonwe = 110;105
  lonee = 125;120

  rad  = 4.0*atan(1.0)/180.0
  z1 = z(:,{latse:latne},{lonwe:lonee})
  lat1 = z1&lat
  clat1 = cos(lat1*rad)
  rws1  = wgt_areaave_Wrap(z1*10^(10), clat1, 1.0, 1)

  latnw = 32;30;25
  latsw = 42;40;35
  lonww = 145;130
  lonew = 160;145

  z2 = z(:,{latsw:latnw},{lonww:lonew})
  lat2 = z2&lat
  clat2 = cos(lat2*rad)
  rws2  = wgt_areaave_Wrap(z2*10^(10), clat2, 1.0, 1)  

;rws1>0 rws2<0.
  rff = stats_frequency67(timi,phase,amp,isf,ief,rws1,rws2)
  rfs = stats_frequency67(timi,phase,amp,iss,ies,rws1,rws2)

  print(rff)
  print(rfs)

  odir = ddir+"/"+"boostrap_rws_omi_p67_region1.nc"
  fout = addfile(odir,"r")           
  rs   = fout->rf
  rs := rs*100

  nrand = dimsizes(rs(0,:))
  rs_me = rs
  do i = 0,nrand - 1
    rs_me(0,i) = rs_me(0,i) - rff
    rs_me(1,i) = rs_me(1,i) - rfs
  end do
  rs_sorted = dim_pqsort_n(rs_me, 2, 1)

  ind_up = toint(nrand*0.95-1)
  ind_dn = toint(nrand*0.05-1)

  sdf_up = rs_me(0,ind_up)
  sdf_dn = rs_me(0,ind_dn)

  sds_up = rs_me(1,ind_up)
  sds_dn = rs_me(1,ind_dn)

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

  res@tmYROn = False
  res@tmXTOn = False
  res@tmXBLabelFontHeightF = 0.012
  res@tmYLLabelFontHeightF = 0.015
  res@tmYLMode = "Explicit" 
  res@tmYLValues = xf
  res@tmYLLabels = (/"~F21~IOfP","~F21~IOsP"/)

  res@vpXF = 0.58-0.025
  res@vpYF = 0.86+0.065
  res@vpWidthF  = 0.35
  res@vpHeightF = 0.18

  res@trYMaxF = 2.35
  res@trYMinF = 0.85
  res@trXMinF = 20.
  res@trXMaxF = 65.
  res@gsnXRefLine  = 10.

  res@gsnXYBarChart = True
  res@gsnXYBarChartBarWidth = 0.28
  res@gsnXYBarChartColors = (/fcolor,scolor/)
  res@xyLineColors = (/fcolor,scolor/)
  res@gsnLeftString = "~F21~ (b) Freq. (EA~B~RWS~N~>0&WP~B~RWS~N~<0)"  
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

  resu@vpXF = 0.15-0.025
  resu@vpYF = 0.95
  ;res@vpWidthF  = 0.3
  ;res@vpHeightF = 0.185
  resu@gsnLeftString = "~F21~ (a) RWS (IOsP)" 
  plot2 = gsn_csm_contour_map(wks, rwss, resu)

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
  tex1 = gsn_add_text(wks, plot2, "EA~B~RWS~N~", lonwe+8, latne-3, txres)
;; MEJS
  lonmin = lonwe
  lonmax = lonee
  latmin = latse
  latmax = latne
  polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
  polyy = (/latmin,latmin,latmax,latmax,latmin/)
  dum2 = gsn_add_polyline(wks, plot2, polyx, polyy, respoly)
  txt2 = gsn_add_text(wks, plot2, "WP~B~RWS~N~", lonww+8,latnw-3, txres)

  draw(plot2)
  draw(plot1)

  delete([/res,cmap,resu,amp,phase,time,timi/])
;;jet
;---------------------------------------------------------
;function: 返回位相3对应的日期 格式：YYYYMMDD
;---------------------------------------------------------
function return_days(timt:integer,xdir)
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
  return nda

end
;----------------------------------------------------------
  udir = ddir+"/"+"ERA5_u200_daily_mean_1979-2022_r2.5.nc"
;;jet events
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"

  lats = 0
  latn = 80
  lonw = 60
  lone = 300
;----------------------------------------------------------
; main code.
;----------------------------------------------------------
  f = addfile(udir,"r")
  u = f->u(:,{200},{lats:latn},{lonw:lone})
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  timi = cd_calendar(time,-2)
  u := runave_n_Wrap(u,5,0,0)
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

  IOfP = return_days(timi,fdir) ;IOfP
  IOsP = return_days(timi,sdir) ;IOsP

  uf = dim_avg_n_Wrap(u(IOfP,:,:),0)
  us = dim_avg_n_Wrap(u(IOsP,:,:),0)

  nf = 31-1
  ns = 14-1

  ud = us - uf
  copy_VarCoords(uf, ud)

  uuf = u(IOfP,:,:)
  doff = equiv_sample_size(uuf(lat|:,lon|:,time|:), 0.1, 0)
  dff  =  toint(wgt_areaave(doff, 1., 1., 0))

  uus = u(IOsP,:,:)
  dofs = equiv_sample_size(uus(lat|:,lon|:,time|:), 0.1, 0)
  dfs  =  toint(wgt_areaave(dofs, 1., 1., 0))

  f_var=dim_variance_n(uuf, 0)
  s_var=dim_variance_n(uus, 0)

;;diff
  undef("confi_d")
  confi_d = ttest(uf,f_var,nf, us,s_var,ns, True, False)
  copy_VarCoords(ud,confi_d)
; u250 resources.
;--plot
  cmap = read_colormap_file("MPL_RdBu") 

  res = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False
  res@gsnAddCyclic = False
  res@gsnMajorLatSpacing  = 20
  res@gsnMajorLonSpacing  = 30
  res@gsnLeftString = ""
  res@gsnRightString = ""

  res@gsnLeftStringOrthogonalPosF  = 0.013
  res@gsnRightStringOrthogonalPosF = 0.013
  res@gsnLeftStringOrthogonalPosF = -0.2
  res@gsnLeftStringFontHeightF  = 0.013
  res@gsnRightStringFontHeightF = 0.013
  res@tmXBLabelFontHeightF  = 0.01
  res@tmYLLabelFontHeightF  = 0.01

  res@mpFillOn = True
  res@mpLandFillColor = "white"
  res@mpMinLatF       = lats
  res@mpMaxLatF       = latn
  res@mpMinLonF       = lonw
  res@mpMaxLonF       = lone

  res@cnFillOn       = True              
  res@cnFillPalette  = cmap(::-1,:)  
  res@cnLinesOn = True    
  res@cnLevelSelectionMode = "ManualLevels"     
  res@cnMinLevelValF       = -5         
  res@cnMaxLevelValF       =  5            
  res@cnLevelSpacingF      =  1
  res@cnInfoLabelPerimOn = False
  res@cnInfoLabelOn   = False
  res@cnLineLabelsOn   = False

  res@lbLabelBarOn = False
  res@pmLabelBarOrthogonalPosF = 0.23
  res@lbTitleOn = True
  res@lbTitleString = "[m s~S~-1~N~]"
  res@lbTitlePosition  = "Right"
  res@lbTitleDirection = "Across"
  res@lbTitleFontHeightF = 0.015
 
  res@tmXBLabelFont = 21
  res@tmYRLabelFont = 21

;dort
  opt                      =True
  opt@gsnShadeFillType     ="pattern"
  opt@gsnShadeLow          =17;打点
  opt@gsnAddCyclic         =True
  opt@cnFillDotSizeF       =0.001

  sres                     =True
  sres@gsnDraw             =False
  sres@gsnFrame            =False
  sres@cnLinesOn           =False
  sres@gsnLeftString       =""
  sres@gsnRightString      =""
  sres@cnLevelSelectionMode="ExplicitLevels"
  sres@cnLevels            =(/0.05,0.01/)
  sres@cnFillPalette       ="GMT_gray"
  sres@cnFillColors        =(/5,7,-1/)
  sres@cnLineLabelsOn      =False
  sres@cnInfoLabelOn       =False
  sres@lbLabelBarOn        =False
  sres@cnFillScaleF        = 1.3
  sres@cnFillDotSizeF      = 0.005

  label = (/"~F21~(c) U200(IOsP)","~F21~(d) U200(IOfP)","~F21~(e) U200(IOsP - IOfP)"/)
  res@vpWidthF = 0.4
  res@vpHeightF= 0.2

  res@cnMinLevelValF       = -80        
  res@cnMaxLevelValF       =  80          
  res@cnLevelSpacingF      =  10
  res@vpXF     = 0.1
  res@vpYF     = 0.62
  res@gsnLeftString  = label(0)
  plotb1 = gsn_csm_contour_map(wks, us, res)

  res@vpXF     = 0.1
  res@vpYF     = 0.45
  res@gsnLeftString  = label(1)
  res@lbLabelBarOn = True
  plotb2 = gsn_csm_contour_map(wks, uf, res)

  res@cnMinLevelValF       = -20        
  res@cnMaxLevelValF       =  20            
  res@cnLevelSpacingF      =  2
  res@cnLinesOn = False 
  res@lbLabelBarOn = True
  res@vpXF     = 0.1
  res@vpYF     = 0.24
  res@gsnLeftString  = label(2)
  plotb3 = gsn_csm_contour_map(wks, ud, res)
  ;base =gsn_csm_contour(wks,confi_d, sres)
  ;base =gsn_contour_shade(base,0.1,1,opt)
  ;overlay(plotb3,base)

  draw(plotb1)
  draw(plotb2)
  draw(plotb3)

  delete([/TIME,year,mon,day,phase,amp,IOsP,IOfP/])
;draw olr
;--------------------------------------------------------------------------
;function: 计算每种事件并返回矩阵的程序。
;--------------------------------------------------------------------------
function calc_pro(timt:integer,is,ie,zt:float,pht:integer,ampt:float)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,igood,pt,at,ip,out,pct1
begin
  dims = dimsizes(zt)
  dim2 = dims(1)
  dim3 = dims(2)
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
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
    ;ip = ind(pt.eq.3.or.pt.eq.2)
    ip = ind(pt.eq.6)
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
;--------------------------------------------------------------------------
;function: 计算每种事件并返回矩阵的程序。
;--------------------------------------------------------------------------
function calc_olr(timt:integer,xdir,ztt:float,pht:integer,ampt:float)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,igood,pt,at,ip,out,pct1,lats,latn,lonw,lone
begin
  lats = -15
  latn =  15
  lonw =  135
  lone =  165
  zt = wgt_areaave_Wrap(ztt(:,{lats:latn},{lonw:lone}), 1., 1., 0)
  ;zt = dim_standardize_n(ztd, 0, 0)
  dims = dimsizes(zt)

  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nf = asciiread(xdir,-1,"string")
  is = str_get_field(nf,1,delim)
  ie = str_get_field(nf,2,delim)
  is := toint(is(1::))
  ie := toint(ie(1::))
  ind_time = ispan(0,dimsizes(timt)-1,1)
  ;event by event.
  pct1 = new((/dimsizes(is),26/),float) 
  pct1 = 0.
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    pt  = pht(its:ite)
    nda = ind_time(its:ite)
    ip  = ind(pt.eq.5)
    ip@_FillValue = -999
    if (any(ip.eq.ip@_FillValue).eq."Missing")
      delete([/ip,its,ite,pt,nda/])
      continue
    else
      do j =0,25
        if (dimsizes(ip).eq.1)
          pct1(i,j) = zt(nda(ip)+j)
        else
          pct1(i,j) = dim_avg_n(zt(nda(ip)+j),0)
        end if
      end do

    end if
    delete([/ip,its,ite,pt,nda/])
  end do

  return pct1

end

;------------------------------------------------------------------
  lats = -25
  latn =  30
  lonw =  50
  lone =  210

  odir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
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
  olr_iosp = calc_olr(timi,sdir,olr,phase,amp)
  olr_iofp = calc_olr(timi,fdir,olr,phase,amp)

  olr_iosp := dim_avg_n(olr_iosp,0)
  olr_iofp := dim_avg_n(olr_iofp,0)

  print(olr_iosp)
  print(olr_iofp)

  IOsP = calc_pro(timi,iss,ies,olr,phase,amp) 

  IOsP!0 = "time"
  IOsP!1 = "lat"
  IOsP!2 = "lon"

;;value for plot.
  Is = dim_avg_n(IOsP,0)
  copy_VarCoords(olr(0,:,:),Is)
;;for ttest.
  undef("ns")
  ns = dimsizes(IOsP(:,0,0))
;(t-test)
  undef("s_var")
  s_var=dim_variance_n(IOsP, 0)
  s_std=sqrt(s_var/ns)
  s_std=where(s_std.eq.0,s_std@_FillValue,s_std)

  undef("t_s")
  t_s=Is/s_std
  confi_s=Is
  confi_s=student_t(t_s, ns-1)

  copy_VarCoords(Is, confi_s)
  Is = smth9_Wrap(Is, 0.5, 0.25, False)
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
  res@gsnLeftStringFontHeightF  = 0.013
  res@gsnLeftStringOrthogonalPosF  = -0.2

  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBLabelFont = 21
  res@tmYLLabelFont = 21
  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01

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
  res@cnMinLevelValF         = -20
  res@cnMaxLevelValF         =  20
  res@cnLevelSpacingF        =  4.
  res@cnLabelBarEndStyle     = "ExcludeOuterBoxes"

  res@lbTitleOn = True
  res@lbTitleString = "[W m~S~-2~N~]"
  res@lbTitlePosition  = "Right"
  res@lbTitleDirection = "Across"
  res@lbTitleFontHeightF = 0.015
  res@lbBoxLinesOn       = True
  res@lbLabelFontHeightF = 0.012
  res@pmLabelBarOrthogonalPosF = 0.2
  res@lbLabelStride = 1
  res@pmLabelBarWidthF = 0.4

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
  sres@cnFillScaleF        = 1.3
  sres@cnFillDotSizeF      = 0.003

  res@vpXF = 0.57
  res@vpYF = 0.6
  res@vpHeightF = 0.25
  res@vpWidthF  = 0.4
  res@lbLabelBarOn  = True
  res@gsnLeftString = "~F21~ (f) OLR (IOsP P6)"
  plots = gsn_csm_contour_map(wks, Is, res) 
  bases=gsn_csm_contour(wks,confi_s, sres)
  bases=gsn_contour_shade(bases,0.1,1,opt)
  overlay(plots,bases)

  resp = True
  resp@gsLineColor = "black"
  resp@gsLineThicknessF = 3.
  resp@gsLineDashPattern = 0

  polyx = (/130,170,170,130,130/)
  polyy = (/-15,-15,15,15,-15/)

  dum1 = gsn_add_polyline(wks, plots, polyx, polyy, resp)

  draw(plots)

;;curves.
  undef("res")
  res = True                      
  res@gsnDraw  = False
  res@gsnFrame = False     
  res@gsnLeftStringFontHeightF    = 0.013
  res@gsnLeftStringOrthogonalPosF = -0.2

  res@tmYROn = False
  res@tmXTOn = False
  res@tmXBLabelFontHeightF = 0.012
  res@tmYLLabelFontHeightF = 0.012

  res@tiXAxisString = "Lag day"
  res@tiYAxisString = "[W m~S~-2~N~]"

  res@trYMaxF =  15.
  res@trYMinF = -15.
  ;res@trXMinF = 20.
  ;res@trXMaxF = 65.
  res@gsnYRefLine = 0.

  res@xyLineThicknesses = (/ 3.0, 3.0/)          
  res@xyLineColors      = (/"gray","gray"/)          

  res@xyDashPattern  = 0                         
  res@xyMarkLineMode = "MarkLines"               
  res@xyMarkers      = (/16,16/)               
  res@xyMarkerColors = (/fcolor,scolor/)

  res@vpXF = 0.6
  res@vpYF = 0.32
  res@vpWidthF  = 0.36
  res@vpHeightF = 0.2

  res@gsnLeftString = "~F21~(g) OLR after P5 "
  plotolr = gsn_csm_xy(wks,ispan(0,dimsizes(olr_iofp)-1,1),\
                       (/olr_iofp,olr_iosp/),res) 

  draw(plotolr)

  i = 0 
  do while(i.eq.0)
    i = 0
  end do