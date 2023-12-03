;--------------------------------------------------------
; draw Hovmueller. 
; shading: u anomalies at 35N.
; vontour: phi anomalies at 20N.
; time: 2023-12-01
;--------------------------------------------------------
;-----------------------------------------------------------

ddir   = "/home/sunming/data5/cuixy/global_jets/data"
pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

phidir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022.nc"
udir   = ddir+"/"+"ERA5_daily_anom_u250_1979-2022.nc"

;--------------------------------------------------------

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

fphi = addfile(phidir,"r")
fu   = addfile(udir,"r")

phi  = fphi->phi(:,{15:25},{0:250})
u    = fu->u(:,{25:35},{0:250})

phi := dim_avg_n_Wrap(phi, 1)
u   := dim_avg_n_Wrap(u, 1)
print("complete data selling "+ systemfunc("date"))

time = phi&time
TIME = cd_calendar(time, 0)
year = toint(TIME(:,0))
mon  = toint(TIME(:,1))
day  = toint(TIME(:,2))

;-------------------------------------------------------
;; selling MJO.

delim = ","
rmmf  = asciiread(rmmdir,-1,"string")
iyear = str_get_field(rmmf,1,delim)
imon  = str_get_field(rmmf,2,delim)
iday  = str_get_field(rmmf,3,delim)
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

;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
llphi = new((/8,41,dimsizes(phi&longitude)/),float)
copy_VarCoords(phi(0,:), llphi(0,0,:))
llphi = 0.
llu   = llphi

do i = 0, 7
    ip = i + 1
    nt = ind(phase.eq.ip.and.amp.gt.1.)
    nind = ind_DJF(nt)
    do j = 0, 40
        ntt   = -15 + j
        nindd = nind + ntt
        llphi(i,j,:) = dim_avg_n(phi(nindd,:), 0)
        llu(i,j,:)   = dim_avg_n(u(nindd,:), 0)
        delete([/ntt,nindd/])
    end do
    delete([/ip,nt,nind/])
end do

llphi = llphi*10^(-6)

print("complete calc leadlag data. "+ systemfunc("date"))

;; resources

;shading T plot
cmap = read_colormap_file("MPL_RdBu")

rest = True
rest@gsnFrame           = False
rest@gsnDraw            = False
rest@gsnMaximize        = False
;rest@gsnAddCyclic       = True
;rest@gsnMajorLatSpacing = 15
rest@gsnMajorLonSpacing  = 30
;rest@gsnCenterStringFontHeightF = 0.025 

rest@vpWidthF  = 0.45
rest@vpHeightF = 0.2

rest@gsnRightString = "u(25-35"+"~S~o~N~N"+")"+" & phi(15-25"+"~S~o~N~N"+")"
rest@gsnRightStringFontHeightF  = 0.012
rest@gsnLeftStringFontHeightF   = 0.015

rest@tmXTOn = False
rest@tmYROn = False

rest@tiYAxisString  = "Lag days"

rest@tmYLMode   = "Explicit"
rest@tmYLValues = (/ 0. , 5., 10., 15., 20., 25., 30., 35., 40. /)
rest@tmYLLabels = (/"-15","-10","-5","0","5","10","15","20","25" /)

;rest@mpFillOn        = False
;rest@mpFillDrawOrder = "Postdraw"
;resolr@mpProjection= "Robinson"
;rest@mpMinLatF   = -30
;rest@mpMaxLatF   = 40
;rest@mpMinLonF   = 0
;rest@mpCenterLonF= 180
;rest@mpMaxLonF   = 360
;rest@mpLandFillColor = "white"
;rest@mpGridAndLimbDrawOrder = "PostDraw"
;rest@mpGridLineColor        = "Gray60"
;rest@mpGridLineThicknessF   = 1.0

rest@cnFillOn               = True
;rest@cnFillMode             = "CellFill"
rest@cnFillPalette          = cmap(::-1,:)
;rest@cnFillOpacityF         = 0.6
rest@cnLinesOn              = False
rest@cnLineLabelsOn         = False
rest@cnInfoLabelOn          = False
rest@cnLevelSelectionMode   = "Manuallevels"
rest@cnMinLevelValF         = -5.
rest@cnMaxLevelValF         =  5.
rest@cnLevelSpacingF        =  1.
rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
rest@lbLabelBarOn           = False
rest@lbTitleString          = "u250 anom (m/s)"
rest@lbTitlePosition        = "Right"
rest@lbTitleDirection       = "Across"
rest@lbTitleFontHeightF     = 0.015
rest@lbBoxLinesOn           = True
rest@lbBoxSeparatorLinesOn  = False
rest@lbLabelFontHeightF     = 0.012

; contour phi plot.
resphi = True
resphi@gsnDraw      = False
resphi@gsnMaximize  = False
resphi@gsnFrame     = False
resphi@gsnAddCyclic = False
resphi@gsnContourNegLineDashPattern = 1
resphi@cnMonoLineColor    = True
resphi@cnLinesOn          = True
resphi@cnLineLabelsOn     = False
resphi@cnLineColor        = "black"
resphi@cnLineThicknessF   = 2.5
resphi@cnInfoLabelPerimOn = False
resphi@cnInfoLabelOn      = False
resphi@cnLevelSelectionMode = "Explicitlevels"
resphi@cnLevels      = (/-10,-8,-6,-4,-2,2,4,6,8,10/)

resP = True
resP@gsnMaximize      = False
resP@gsnFrame         = False
resP@gsnDraw          = True
resP@gsnPanelLabelBar = True
resP@gsnPanelTop      = 0.98
resP@gsnPanelBottom   = 0.02

resP@pmLabelBarOrthogonalPosF = 0.0
resP@pmLabelBarWidthF         = 0.85

resP@lbTitleString     = "u250 anom (m/s)"
resP@lbTitlePosition   = "Right"
resP@lbTitleDirection  = "Across"
resP@lbTitleFontHeightF= 0.015
resP@lbBoxLinesOn      = True
resP@lbBoxSeparatorLinesOn  = False
;resP@lbLabelStride      = 2
;resP@lbFillOpacityF     = 0.6
resP@gsnPanelRowSpec    = True
resP@lbLabelFontHeightF = 0.012

print("complete res. "+ systemfunc("date"))

;; plot.

wks = gsn_open_wks("eps", pdir+"/"+"Hovmueller_u_phi")

;plott   = new(1,graphic)
;plotphi = plott

;plott   = gsn_csm_hov(wks, llu(i,:,:), rest)
;plotphi = gsn_csm_hov(wks, llphi(i,:,:), resphi)
;overlay(plott,plotphi)
;ployx = (/10,10/)
;ployy = (/0.,40./)

;str = unique_string("fill_lft")
;resp = True
;ployx = (/10,10/)
;resp@gsLineColor = (/0.,79.,0.,1./)
;resp@gsLineThicknessF = 4.0
;plott@$str$ = gsn_add_polyline(wks, plott, ployx, ployy, resp)

;str = unique_string("fill_lft")
;ployx = (/70,70/)
;plott@$str$ = gsn_add_polyline(wks, plott, ployx, ployy, resp)

;str = unique_string("fill_lft")
;ployx = (/110,110/)
;resp@gsLineDashPattern = 1
;plott@$str$ = gsn_add_polyline(wks, plott, ployx, ployy, resp)

;str = unique_string("fill_lft")
;ployx = (/180,180/)
;plott@$str$ = gsn_add_polyline(wks, plott, ployx, ployy, resp)

;draw(plott)


label = (/"(a) Phase 1","(b) Phase 2","(c) Phase 3","(d) Phase 4",\
         "(e) Phase 5","(f) Phase 6","(g) Phase 7","(h) Phase 8"/)

print("Prepare labels done. "+ systemfunc("date"))

plott   = new(8,graphic)
plotphi = plott

do i = 0, 7
  rest@gsnLeftString = label(i)
  plott(i)  = gsn_csm_hov(wks, llu(i,:,:), rest)
  plotphi(i) = gsn_csm_hov(wks, llphi(i,:,:), resphi)
  overlay(plott(i),plotphi(i))
  ployx = (/10,10/)
  ployy = (/0.,40./)

  resp = True
  ployx = (/10,10/)
  resp@gsLineColor = (/0.,79.,0.,1./)
  resp@gsLineDashPattern = 0
  resp@gsLineThicknessF = 4.0
  dum1 = gsn_add_polyline(wks, plott(i), ployx, ployy, resp)
  overlay(plott(i),dum1)

  ployx = (/70,70/)
  dum2  = gsn_add_polyline(wks, plott(i), ployx, ployy, resp)
  overlay(plott(i),dum2)
  
  ployx = (/110,110/)
  resp@gsLineDashPattern = 1
  dum3  = gsn_add_polyline(wks, plott(i), ployx, ployy, resp)
  overlay(plott(i),dum3)
  
  ployx = (/180,180/)
  dum4  = gsn_add_polyline(wks, plott(i), ployx, ployy, resp)
  overlay(plott(i),dum4)
end do
;;
gsn_panel(wks, (/plott(0),plott(4),plott(1),plott(5),\
                plott(2),plott(6),plott(3),plott(7)/), (/2,2,2,2/), resP)
print("Draw done. "+ systemfunc("date"))

delete(wks)
