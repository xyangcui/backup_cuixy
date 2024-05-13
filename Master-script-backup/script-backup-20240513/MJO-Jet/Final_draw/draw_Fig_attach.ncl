;--------------------------------------------------------
; draw uwind anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour. 
; waf vector
; This is for test sig. 95% by monte carlo method.
; time: 2023-12-18
; 2.5x2.5.
;--------------------------------------------------------

;;----------------------------------------------------------
function siglevel(ovalue,mvalue)
begin

  rv = ovalue
  lati = rv&lat
  loni = rv&lon

  do i = 0, dimsizes(lati)-1
    do j = 0, dimsizes(loni)-1
      if(ovalue(i,j).lt.mvalue(0,i,j).or.ovalue(i,j).gt.mvalue(1,i,j)) then
        rv(i,j) = 2.
      else
        rv(i,j) = -2.
      end if
    end do
  end do

  return rv

end
;-------------------------------------------------------------
;function calc_lagvalue
;-------------------------------------------------------------
function calc_lagvalue(x,phaset,ampt,ind_DJFt,p1t,p2t,nday)
local llx,x,phaset,ampt,nt,nind,ntt,nindd
begin
  llx = new((/2,dimsizes(x&lat),dimsizes(x&lon)/),float)
  copy_VarCoords(x(0,:,:), llx(0,:,:))
  llx = 0.

  nt = ind(phaset.eq.p1t.and.ampt.gt.1.)
  nind = ind_DJFt(nt)
  ntt   = nday-2
  nindd = nind + ntt
  do ntt = nday-1,nday+2
    nit = nind + ntt
    nindd := array_append_record(nindd, nit, 0)
    delete(nit)
  end do
  llx(0,:,:)   = dim_avg_n(x(nindd,:,:), 0)
  delete([/nt,nind,ntt,nindd/])

  nt = ind(phaset.eq.p2t.and.ampt.gt.1.)
  nind = ind_DJFt(nt)
  ntt   = nday-2
  nindd = nind + ntt
  do ntt = nday-1,nday+2
    nit = nind + ntt
    nindd := array_append_record(nindd, nit, 0)
    delete(nit)
  end do
  llx(1,:,:) = dim_avg_n(x(nindd,:,:), 0)

  return(llx)
end
;-----------------------------------------------------------

  ddir   = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  udir   = ddir+"/"+"ERA5_daily_anom_u250_1979-2022_2.5x2.5.nc"
  vdir   = ddir+"/"+"ERA5_daily_anom_v250_1979-2022_2.5x2.5.nc"

  olrdir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
  ucdir  = ddir+"/"+"ERA5_daily_u250_1979-2022.nc"

  sigdir = "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/Fig3_uv_sf_sig.nc"

  mlats = -15
  mlatn =  60
  mlonw =   0
  mlone =  105

  elats = -15
  elatn =  60
  elonw =  105
  elone =  190

  nlats = -15
  nlatn =  60
  nlonw =  190
  nlone =  340

;--------------------------------------------------------
; u250
;--------------------------------------------------------
  f    = addfile(ucdir,"r")

  time = f->time
  TIME = cd_calendar(time, 0)
  mon  = toint(TIME(:,1))
  nt   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)

  uc    = dim_avg_n_Wrap(f->u(nt,0,:,:),0)

;contour plot
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
  resc@cnLineColor        = "black";magenta
  resc@cnLineThicknessF   = 2.5
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels      = (/30.,40.,46.,60.,70./)

  delete([/time,f,TIME,mon,nt/])

  print("draw clim.u complete. "+ systemfunc("date"))

;; common part
  fu   = addfile(udir,"r")
  fv   = addfile(vdir,"r")
  folr = addfile(olrdir,"r")
;-------------------------------------------------------
;; selling MJO.
  time = fu->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  delim  = ","
  rmmf   = asciiread(rmmdir,-1,"string")
  iyear  = str_get_field(rmmf,1,delim)
  imon   = str_get_field(rmmf,2,delim)
  iday   = str_get_field(rmmf,3,delim)
  phaset = str_get_field(rmmf,6,delim)
  ampt   = str_get_field(rmmf,7,delim)

  iyear  := toint(iyear(1::))
  imon   := toint(imon(1::))
  iday   := toint(iday(1::))
  phaset := toint(phaset(1::))
  ampt   := tofloat(ampt(1::))

  amp   = new((/dimsizes(time)/),float)
  phase = new((/dimsizes(time)/),integer)
  amp   = 0.
  phase = 0

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
    amp(i)   = ampt(nt)
    phase(i) = phaset(nt)
    delete(nt)
  end do

  nn = ispan(0, dimsizes(time)-1, 1)
  amp!0   = "timeind"
  phase!0 = "timeind"
  amp&timeind   = nn
  phase&timeind = nn
  delete([/nn,ampt,phaset,iyear,imon,iday/])

  print("complete RMM selling "+ systemfunc("date"))
;-------------------------------------------------------------

  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.2.and.day.eq.28)

  amp   := amp(n1:n2)
  phase := phase(n1:n2)
  mon   := mon(n1:n2)
  time  := time(n1:n2)

  nt     = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  amp   := amp(nt)
  phase := phase(nt)

  delete(nt)
  ind_DJF = amp&timeind   ; indice of DJF between 1979-2022 DJF.
;;res

;----------------------------------------------------------------------
;; resources
  cmap = read_colormap_file("MPL_RdBu")  ;MPL_RdBu
;shading u plot
  resu = True
  resu@gsnFrame           = False
  resu@gsnDraw            = False
  resu@gsnMaximize        = False
  resu@gsnAddCyclic       = False
  resu@gsnMajorLatSpacing = 15
  resu@gsnMajorLonSpacing = 30
  resu@gsnRightStringFontHeightF = 0.03
  resu@gsnLeftStringFontHeightF  = 0.03

  resu@gsnRightStringOrthogonalPosF = 0.01
  resu@gsnLeftStringOrthogonalPosF  = 0.01

  resu@tmXTOn = False
  resu@tmYROn = False
  resu@tmXBLabelFont = 21
  resu@tmYLLabelFont = 21
  resu@tmXBLabelFontHeightF = 0.02
  resu@tmYLLabelFontHeightF = 0.02
  resu@mpFillOn        = False
  resu@mpLandFillColor = "white"
  resu@mpGridAndLimbDrawOrder = "PostDraw"
  resu@mpGridLineColor        = "Gray60"

  resu@cnFillOn               = True
  resu@cnFillPalette          = cmap(::-1,:);(66:146,:)
  resu@cnLinesOn              = False
  resu@cnLineLabelsOn         = False
  resu@cnLevelSelectionMode   = "Manuallevels"
  resu@cnMinLevelValF         = -6.
  resu@cnMaxLevelValF         =  6.
  resu@cnLevelSpacingF        =  1.
  resu@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  resu@lbLabelBarOn           = False
; contour olr plot.
  resolr = True
  resolr@gsnDraw      = False
  resolr@gsnMaximize  = False
  resolr@gsnFrame     = False
  resolr@gsnAddCyclic = False
  resolr@gsnContourNegLineDashPattern = 1
  resolr@cnMonoLineColor    = True
  resolr@cnLinesOn          = True
  resolr@cnLineLabelsOn     = False
  resolr@cnLineColor        = "red3"
  resolr@cnLineThicknessF   = 3.5
  resolr@cnInfoLabelPerimOn = False
  resolr@cnInfoLabelOn      = False
  resolr@cnLevelSelectionMode = "Explicitlevels"
  resolr@cnLevels      = (/-25,-20,-15,-10,-5,5,10,15,20,25/)

  sres   = True 
  sres@cnLineDrawOrder    = "PostDraw"                    
  sres@gsnDraw            = False                          
  sres@gsnFrame           = False                        
  sres@cnLineLabelsOn     = False                              
  sres@cnLinesOn          = False                           
  sres@cnInfoLabelOn      = False
  sres@cnFillOn           = False                       
  sres@cnInfoLabelOn      = False
  sres@cnFillScaleF       = 0.5
  sres@cnFillDotSizeF     = 0.003
  sres@gsnLeftString      =" "
  sres@gsnRightString     =" "
  sres@cnMonoFillColor    = True
  sres@cnFillColor        = "gray"

  opt = True
  opt@gsnShadeFillType   = "pattern"
  opt@gsnShadeHigh       = 17
  opt@gsnShadeFillScaleF = 1.5

  respoly = True
  respoly@gsLineColor      = "purple"
  respoly@gsLineThicknessF = 4.5

  print("complete res. "+ systemfunc("date"))
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------
  olr  = folr->olr
  u    = fu->u
  print("complete data selling "+ systemfunc("date"))
;---------------------------------------------------------------------
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
;; MEJS
  lats = 15
  latn = 60
  lonw = mlonw
  lone = mlone
  p1 = 7
  p2 = 3

  llum = calc_lagvalue(u(:,{lats:latn},{lonw:lone}),phase,amp,ind_DJF,p1,p2,5)
  lats = -20
  latn =  20
  llom = calc_lagvalue(olr(:,{lats:latn},{lonw:lone}),phase,amp,ind_DJF,p1,p2,5)
  llom := smth9_Wrap(llom, 0.5, 0.25, False)
;; EAJS
  lats = 15
  latn = 60
  lonw = elonw
  lone = elone
  p1 = 7
  p2 = 3

  llue = calc_lagvalue(u(:,{lats:latn},{lonw:lone}),phase,amp,ind_DJF,p1,p2,20)
  lats = -20
  latn =  20
  lloe = calc_lagvalue(olr(:,{lats:latn},{lonw:lone}),phase,amp,ind_DJF,p1,p2,20)
  lloe := smth9_Wrap(lloe, 0.5, 0.25, False)
;; NAJS
  lats = 15
  latn = 60
  lonw = nlonw
  lone = nlone
  p1 = 4
  p2 = 8

  llun = calc_lagvalue(u(:,{lats:latn},{lonw:lone}),phase,amp,ind_DJF,p1,p2,5)
  lats = -20
  latn =  20
  llon = calc_lagvalue(olr(:,{lats:latn},{lonw:lone}),phase,amp,ind_DJF,p1,p2,5)
  llon := smth9_Wrap(llon, 0.5, 0.25, False)
  print("complete calc leadlag data. "+ systemfunc("date"))
;----------------------------------------------------------------------
;; plot.
  wks = gsn_open_wks("x11", pdir+"/"+"Fig3_U_SF_T_MEJS")

  plotu1 = new(2,graphic)
  ploto1 = plotu1

  plotu2  = plotu1
  ploto2  = ploto1

  plotu3  = plotu1
  ploto3  = ploto1

  Label  = (/"~F21~(b) 250hPa_U&OLR","~F21~(a) 250hPa_U&OLR"/)
  R1 = (/"~F21~p7 day5","~F21~p3 day5"/)
  R2 = (/"~F21~p7 day20","~F21~p3 day20"/)
  R3 = (/"~F21~p4 day5","~F21~p8 day5"/)

;;MEJS
  lats = mlats
  latn = mlatn
  lonw = mlonw
  lone = mlone

  do n = 0,1

    y  = llum(n,:,:)
    z  = llom(n,:,:)

    resu@mpMinLatF  = lats
    resu@mpMaxLatF  = latn
    resu@mpMinLonF  = lonw
    resu@mpMaxLonF  = lone
    resu@gsnLeftString = Label(n)
    resu@gsnRightString = R1(n)

    resu@cnMinLevelValF  = -5.
    resu@cnMaxLevelValF  =  5.
    resu@cnLevelSpacingF =  1.

    plotu1(n) = gsn_csm_contour_map(wks, y, resu)

    ploto1(n) = gsn_csm_contour(wks, z, resolr)

    plotut = gsn_csm_contour(wks, uc, resc)

    overlay(plotu1(n),ploto1(n))
    overlay(plotu1(n),plotut)
    delete([/y,z/])
  end do

;;EAJS
  lats = elats
  latn = elatn
  lonw = elonw
  lone = elone

  resu@tmYLOn = False

  do n = 0,1

    y  = llue(n,:,:)
    z  = lloe(n,:,:)

    resu@mpMinLatF  = lats
    resu@mpMaxLatF  = latn
    resu@mpMinLonF  = lonw
    resu@mpMaxLonF  = lone

    resu@cnMinLevelValF  = -5.
    resu@cnMaxLevelValF  =  5.
    resu@cnLevelSpacingF =  1.

    resu@gsnLeftString = ""
    resu@gsnRightString = R2(n)

    plotu2(n) = gsn_csm_contour_map(wks, y, resu)

    ploto2(n) = gsn_csm_contour(wks, z, resolr)

    plotut = gsn_csm_contour(wks, uc, resc)

    overlay(plotu2(n),ploto2(n))
    overlay(plotu2(n),plotut)
    delete([/y,z/])
  end do

;;NAJS
  lats = nlats
  latn = nlatn
  lonw = nlonw
  lone = nlone

  do n = 0,1

    y  = llun(n,:,:)
    z  = llon(n,:,:)

    resu@mpMinLatF  = lats
    resu@mpMaxLatF  = latn
    resu@mpMinLonF  = lonw
    resu@mpMaxLonF  = lone

    resu@cnMinLevelValF  = -5.
    resu@cnMaxLevelValF  =  5.
    resu@cnLevelSpacingF =  1.

    resu@gsnLeftString = ""
    resu@gsnRightString = R3(n)

    plotu3(n) = gsn_csm_contour_map(wks, y, resu)

    ploto3(n) = gsn_csm_contour(wks, z, resolr)

    plotut = gsn_csm_contour(wks, uc, resc)

    overlay(plotu3(n),ploto3(n))
    overlay(plotu3(n),plotut)
    delete([/y,z/])
  end do

  res1 = True

  res1@gsnAttachPlotsXAxis = False 
  ;res1@gsnAttachBorderOn   = False 

  res2 = True
  ;res2@gsnAttachBorderOn   = False

  id1 = gsn_attach_plots(plotu1(1),(/plotu2(1),plotu3(1)/),res1,res2)
  id2 = gsn_attach_plots(plotu1(0),(/plotu2(0),plotu3(0)/),res1,res2)

  resP = True
  resP@gsnMaximize      = True
  resP@gsnFrame         = False
  resP@gsnDraw          = True

  resP@pmLabelBarOrthogonalPosF = 0.0
  resP@pmLabelBarParallelPosF   = 0.0
  resP@pmLabelBarWidthF         = 0.5

  resP@lbTitleOn = True
  resP@lbTitleString = "~F21~[m/s]"
  resP@lbTitlePosition = "Right"
  resP@lbTitleDirection = "Across"
  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012

  resP@gsnPanelYWhiteSpacePercent = 5

  resP@gsnPanelLabelBar = True
  gsn_panel(wks,(/plotu1(1),plotu1(0)/),(/1,1/),resP)

  i = 0
  do while(i.eq.0)
    i= 0
  end do