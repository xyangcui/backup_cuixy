
;----------------------------------------------------------------------------------
; plot DJF zonal wind at 250 hPa.
; time: 2023-11-05
;----------------------------------------------------------------------------------
datadir = "/home/sunming/data5/cuixy/global_jets/data/ERA5_daily_u250_1979-2022.nc"
plotdir = getenv("plotu")

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

;X    = u
;dims = dimsizes(X)  
;x1d  = ndtooned(X)      ; convert 2D array to 1D for use in maxind
;inds = ind_resolve(maxind (x1d), dims)    ; convert 1D array back to 2D  
;ilat = inds(0,0)        ; select the latitude index where the X array is at its' maximum   
;ilon = inds(0,1)        ; select the longitude index where the X array is at its' maximum  
;lat_max = X&lat(ilat)   ; insert the latitude index into the lat coordinate variable
;lon_max = X&lon(ilon)   ; insert the longitude index into the lon coordinate variable  
;print("Maximum value located at "+lat_max+", "+lon_max)

;--plot
res1 = True
res1@gsnFrame     = False
res1@gsnDraw      = False
res1@gsnMaximize  = False
res1@gsnAddCyclic = True
res1@gsnMajorLatSpacing = 15
res1@gsnMajorLonSpacing = 30

;res1@vpXF     = 0.2
;res1@vpYF     = 0.98
;res1@vpWidthF = 0.7
;res1@vpHeightF= 0.2

res1@mpFillOn = True
res1@mpLandFillColor = "white"
res1@mpGridLineColor       = "gray"
res1@mpGridLineDashPattern = 1
res1@mpMinLatF     = -20
res1@mpMaxLatF     = 85
res1@mpMinLonF     = -95
res1@mpMaxLonF     = 70

cmap1  = read_colormap_file("MPL_jet")
cmap2 = cmap1(49::,:)
res1@cnFillOn = True
res1@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
res1@cnFillPalette        = cmap2
res1@cnLinesOn            = False
res1@cnLineLabelsOn       = False
res1@cnLevelSelectionMode = "Manuallevels"
res1@cnMinLevelValF       =  5.
res1@cnMaxLevelValF       =  70.
res1@cnLevelSpacingF      =  1.

res1@lbBoxLinesOn          = True
res1@lbBoxSeparatorLinesOn = False
res1@lbLabelStride         = 2
res1@lbTitleString         = "~F21~(unit: m/s)"
res1@lbTitlePosition       = "Bottom"
res1@lbTitleFontHeightF    = 0.012
res1@lbOrientation         = "Vertical"
;res1@pmLabelBarHeightF     = 0.2

res1@gsnRightString  = ""
res1@gsnLeftString   = "~F22~a ~F21~Zonal wind at 250 hPa"


wks = gsn_open_wks("eps",plotdir)

plot  = gsn_csm_contour_map(wks, u, res1)

respoly = True
respoly@gsLineColor      = "navyblue"
respoly@gsLineThicknessF = 2.0

lonmin = -90
lonmax = 0
latmin = 15
latmax = 75
polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
polyy = (/latmin,latmin,latmax,latmax,latmin/)
dum1  = gsn_add_polyline(wks, plot, polyx, polyy, respoly)

respoly@gsLineColor      = "mediumvioletred"

lonmin = -75
lonmax = 50
latmin = 10
latmax = 80
polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
polyy = (/latmin,latmin,latmax,latmax,latmin/)
dum2  = gsn_add_polyline(wks, plot, polyx, polyy, respoly)

draw(plot)

end