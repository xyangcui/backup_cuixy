
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
  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  udir = ddir+"/"+"ERA5_u200_daily_mean_1979-2022_r2.5.nc"
;;jet events
  ndir = ddir+"/"+"IOnP_date_omi.csv"
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"
  idir = ddir+"/"+"IOiP_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"

  lats = 0
  latn = 80
  lonw = 60
  lone = 300

  pdir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/plot/"
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

  nf = 30
  ns = 13

  ud = us - uf
  copy_VarCoords(uf, ud)

  uuf = u(IOfP,:,:)
  doff = equiv_sample_size(uuf(lat|:,lon|:,time|:), 0.05, 0)
  dff  =  toint(wgt_areaave(doff, 1., 1., 0))

  uus = u(IOsP,:,:)
  dofs = equiv_sample_size(uus(lat|:,lon|:,time|:), 0.05, 0)
  dfs  =  toint(wgt_areaave(dofs, 1., 1., 0))

  print(dff)
  print(dfs)

  f_var=dim_variance_n(uuf, 0)
  s_var=dim_variance_n(uus, 0)
;;diff
  confi_d = ttest(uf,f_var,nf, us,s_var,ns, False, False)
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

  res@gsnLeftStringOrthogonalPosF  = 0.012
  res@gsnRightStringOrthogonalPosF = 0.012
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringFontHeightF = 0.015
  res@tmXBLabelFontHeightF  = 0.012
  res@tmYLLabelFontHeightF  = 0.012

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
  res@cnInfoLabelOns   = False
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

  wks = gsn_open_wks("eps",pdir+"fig.u200_clim_two_episodes")

  label = (/"~F21~(a) U200(IOsP)","~F21~(b) U200(IOfP)","~F21~(c) U200(IOsP - IOfP)"/)
  res@vpWidthF = 0.6
  res@vpHeightF= 0.2

  res@cnMinLevelValF       = -80        
  res@cnMaxLevelValF       =  80          
  res@cnLevelSpacingF      =  10
  res@vpXF     = 0.2
  res@vpYF     = 0.95
  res@gsnLeftString  = label(0)
  plot1 = gsn_csm_contour_map(wks, us, res)

  res@vpXF     = 0.2
  res@vpYF     = 0.68
  res@gsnLeftString  = label(1)
  res@lbLabelBarOn = True
  plot2 = gsn_csm_contour_map(wks, uf, res)

  res@cnMinLevelValF       = -20        
  res@cnMaxLevelValF       =  20            
  res@cnLevelSpacingF      =  2
  res@cnLinesOn = False 
  res@lbLabelBarOn = True
  res@vpXF     = 0.2
  res@vpYF     = 0.35
  res@gsnLeftString  = label(2)
  plot3 = gsn_csm_contour_map(wks, ud, res)
  base =gsn_csm_contour(wks,confi_d, sres)
  base =gsn_contour_shade(base,0.1,1,opt)
  overlay(plot3,base)

  draw(plot1)
  draw(plot2)
  draw(plot3)

  ;i = 0
  ;do while (i.eq.0)
  ;  i = 0
  ;end do
