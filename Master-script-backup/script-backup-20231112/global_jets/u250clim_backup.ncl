
;----------------------------------------------------------------------------------
; plot DJF zonal wind at 250 hPa.
; time: 2023-11-05
;----------------------------------------------------------------------------------
datadir = "/home/sunming/data5/cuixy/global_jets/data/ERA5_daily_u250_1979-2022.nc"
plotdir = "/home/sunming/data5/cuixy/global_jets/plot/"

;----------------------------------------------------------------------------------
; Main code.
;----------------------------------------------------------------------------------
begin
f    = addfile(datadir,"r")

time = f->time
TIME = cd_calendar(time, 0)
mon  = toint(TIME(:,1))
nt   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)

u    = dim_avg_n_Wrap(f->u(nt,0,:,:),0)

;--plot
res = True
res@gsnFrame     = False
res@gsnDraw      = False
res@gsnMaximize  = False
res@gsnAddCyclic = True
res@gsnMajorLatSpacing = 15
res@gsnMajorLonSpacing = 30

res@vpXF     = 0.2
res@vpYF     = 0.85
res@vpWidthF = 0.7
res@vpHeightF= 0.8

res@mpFillOn = True
res@mpLandFillColor = "white"
res@mpGridLineColor       = "gray"
res@mpGridLineDashPattern = 1
res@mpCenterLonF  = 150

cmap  = read_colormap_file("MPL_jet")
cmap2 = cmap(49::,:)
res@cnFillOn = True
res@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
res@cnFillPalette        = cmap2
res@cnLinesOn            = False
res@cnLineLabelsOn       = False
res@cnLevelSelectionMode = "Manuallevels"
res@cnMinLevelValF       =  5.
res@cnMaxLevelValF       =  70.
res@cnLevelSpacingF      =  1.

res@lbBoxLinesOn          = True
res@lbBoxSeparatorLinesOn = False
res@lbLabelStride         = 5
res@lbTitleString         = "~F21~(unit: m/s)"
res@lbTitlePosition       = "Bottom"
res@lbTitleFontHeightF    = 0.015
res@pmLabelBarOrthogonalPosF = 0.1

res@gsnRightString  = ""
res@gsnLeftString   = "~F21~Zonal wind at 250 hPa"
;contour
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
resc@cnLineColor        = "black"
resc@cnLineThicknessF   = 2.0
resc@cnInfoLabelPerimOn = False
resc@cnInfoLabelOn      = False
resc@cnLevelSelectionMode = "Explicitlevels"
resc@cnLevels      = (/30.,35.,40.,45./)

wks = gsn_open_wks("eps","/home/sunming/data5/cuixy/global_jets/plot/u250clim")

plot  = gsn_csm_contour_map(wks, u, res)
plotc = gsn_csm_contour(wks, u, resc)
overlay(plot, plotc)
draw(plot)


end