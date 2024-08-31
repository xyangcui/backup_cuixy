;-----------------------------------------------------
; draw coPhase for All, EI, LA and NE.
; MJO phases color table
; colors = (/"indianred3","darkorange2","darkgoldenrod2","lightgoldenrod2",\
; "lightgoldenrodyellow","darkseagreen3","darkturquoise","deepskyblue3"/)
;
; time: 20240528
;-----------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/global_jets/data"
  ucdir = ddir+"/ERA5_daily_u250_1979-2022_r2.5.nc"
  sdir  = ddir+"/monte_carlo"+"/MLR_sig.nc"

;;----------------------------------------------------------
function siglevel(ovalue,mvalue)
begin
  rv = ovalue
  lati = rv&lat
  loni = rv&lon
  do i = 0, dimsizes(lati)-1
    do j = 0, dimsizes(loni)-1

      if(ovalue(i,j).gt.mvalue(i,j)) then
        rv(i,j) =  2.
      else
        rv(i,j) = -2.
      end if

    end do
  end do
  return rv
end
;-----------------------------------------------------
; Main code.
;-----------------------------------------------------
  fuc = addfile(ucdir,"r")
  uc  = fuc->u(:,0,:,:)

  time = uc&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)
  year := year(n1:n2)
  mon  := mon(n1:n2)
  uc   := uc(n1:n2,:,:)

  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  year := year(nt)
  mon  := mon(nt)
  uc   := uc(nt,:,:)

;; sig value.
  fs = addfile(sdir,"r")
  R = fs->R2
  R95 = R(474,:,:)
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
  resc@cnLineColor        = "magenta";magenta
  resc@cnLineThicknessF   = 2.5
  resc@cnInfoLabelPerimOn = False
  resc@cnInfoLabelOn      = False
  resc@cnLevelSelectionMode = "Explicitlevels"
  resc@cnLevels = (/30.,45.,60.,75./)

  f1 = addfile(ddir+"/MLRegression.nc","r")
  coAmp1 = f1->coAmp({0:60},:)
  R2 = f1->R2({0:60},:)
;; draw plot  
;;res.
  res = True
  res@gsnFrame     = False
  res@gsnDraw      = False
  res@gsnMaximize  = False
  res@gsnAddCyclic = True
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30

  res@gsnLeftStringOrthogonalPosF  = 0.012
  res@gsnRightStringOrthogonalPosF = 0.012
  res@gsnLeftStringFontHeightF  = 0.016
  res@gsnRightStringFontHeightF = 0.016
  res@tmXBLabelFontHeightF  = 0.012
  res@tmYLLabelFontHeightF  = 0.012

  res@mpFillOn = True
  res@mpLandFillColor = "white"
  res@mpCenterLonF    = 180.
  res@mpMinLatF       = 0
  res@mpMaxLatF       = 60
  res@mpMinLonF       = 0
  res@mpMaxLonF       = 360

  res@cnFillOn = True
  res@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
  res@cnFillPalette        = "MPL_Blues";
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "Manuallevels"
  res@cnMinLevelValF       =  0;0.
  res@cnMaxLevelValF       =  3;3.5
  res@cnLevelSpacingF      =  0.25;0.25
  res@lbLabelBarOn = True
  res@pmLabelBarOrthogonalPosF = 0.25

  res@tmXBLabelFont = 21
  res@tmYRLabelFont = 21
  res@gsnRightString = ""

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


  wks  = gsn_open_wks("x11","maponly")
;; NE
  u = dim_avg_n_Wrap(uc,0)

  res@gsnLeftString  = "~F21~(b) coAmp"
  res@vpYF = 0.95
  plot1  = gsn_csm_contour_map(wks, coAmp1, res) 
  plotc1 = gsn_csm_contour(wks, u, resc)

  xt = siglevel(R2,R95)
  dum1 = gsn_csm_contour(wks,xt,sres)
  dum1 = gsn_contour_shade(dum1,1.9,-999.,opt)

  overlay(plot1,dum1)
  overlay(plot1,plotc1)

  draw(plot1)  

  i = 0
  do while (i.eq.0)
    i = 0
  end do