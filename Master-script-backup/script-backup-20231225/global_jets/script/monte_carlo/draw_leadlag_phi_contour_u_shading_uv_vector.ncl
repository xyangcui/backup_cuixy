;--------------------------------------------------------
; draw uwind anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour.
; This is for test sig. 95% by monte carlo method.
; time: 2023-12-17
; 2.5x2.5.
;--------------------------------------------------------

ddir   = "/home/sunming/data5/cuixy/global_jets/data"
pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

phidir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022_2.5x2.5.nc"
udir   = ddir+"/"+"ERA5_daily_anom_u250_1979-2022_2.5x2.5.nc"
vdir   = ddir+"/"+"ERA5_daily_anom_v250_1979-2022.nc"
olrdir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
ucdir  = ddir+"/"+"ERA5_daily_u250_1979-2022.nc"
sigdir = ddir+"/"+"monte_carlo/phiu_sig_phase15.nc"
usigdir = ddir+"/"+"monte_carlo/u_sig_phase15.nc"
phisigdir = ddir+"/"+"monte_carlo/phi_sig_phase15.nc"

latmin = -30
latmax = 90
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
resc@cnLineColor        = "magenta"
resc@cnLineThicknessF   = 3.0
resc@cnInfoLabelPerimOn = False
resc@cnInfoLabelOn      = False
resc@cnLevelSelectionMode = "Explicitlevels"
resc@cnLevels      = (/30.,40./)

delete([/time,f,TIME,mon,nt/])

print("draw clim.u complete. "+ systemfunc("date"))
;--------------------------------------------------------

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

fphi = addfile(phidir,"r")
fu   = addfile(udir,"r")
fv   = addfile(vdir,"r")
folr = addfile(olrdir,"r")

olr  = folr->olr(time|:,lat|:,lon|:)
phi  = fphi->phi(:,{latmin:latmax},:)
u    = fu->u(:,{latmin:latmax},:)
v    = fv->v(:,{latmin:latmax},:)

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

;; sig
sigl = 0.05
sigh = 0.95

;fsig = addfile(sigdir,"r")
;m_value_phi = fsig->m_valuephi(:,:,:,:,:)
;m_value_u   = fsig->m_valueu(:,:,:,:,:)

fsigphi = addfile(phisigdir,"r")
fsigu   = addfile(usigdir,"r")
m_value_phi = fsigphi->m_valuephi(:,:,:,:,:)
m_value_u   = fsigu->m_valueu(:,:,:,:,:)

m_valuephi  = new((/2,2,8,dimsizes(phi&lat),dimsizes(phi&lon)/),float) 
m_valuephi  = 0.
m_valueu  = new((/2,2,8,dimsizes(u&lat),dimsizes(u&lon)/),float) 
m_valueu  = 0.

nl  = toint(dimsizes(m_value_phi(:,0,0,0,0))*sigl) - 1
nh  = toint(dimsizes(m_value_phi(:,0,0,0,0))*sigh) - 1
m_valuephi(0,:,:,:,:) = m_value_phi(nl,:,:,:,:)
m_valuephi(1,:,:,:,:) = m_value_phi(nh,:,:,:,:)

m_valuephi = m_valuephi*10^(-6)

nl = toint(dimsizes(m_value_u(:,0,0,0,0))*sigl) - 1
nh = toint(dimsizes(m_value_u(:,0,0,0,0))*sigh) - 1
m_valueu(0,:,:,:,:) = m_value_u(nl,:,:,:,:)
m_valueu(1,:,:,:,:) = m_value_u(nh,:,:,:,:)

delete([/m_value_phi,m_value_u/])
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
llphi = new((/8,8,dimsizes(phi&lat),dimsizes(phi&lon)/),float)
copy_VarCoords(phi(0,:,:), llphi(0,0,:,:))
llphi = 0.

llu = new((/8,8,dimsizes(u&lat),dimsizes(u&lon)/),float)
copy_VarCoords(u(0,:,:), llu(0,0,:,:))
llu = 0.

llolr = new((/8,8,dimsizes(olr&lat),dimsizes(olr&lon)/),float)
copy_VarCoords(olr(0,:,:), llolr(0,0,:,:))
llolr = 0.

do i = 0, 7
    ip = i + 1
    nt = ind(phase.eq.ip.and.amp.gt.1.)
    nind = ind_DJF(nt)
    do j = 0, 7
        ntt   = -15 + j*5
        nindd = nind + ntt
        llphi(i,j,:,:) = dim_avg_n(phi(nindd,:,:), 0)
        llu(i,j,:,:)   = dim_avg_n(u(nindd,:,:), 0)
        llolr(i,j,:,:) = dim_avg_n(olr(nindd,:,:), 0)
        delete([/ntt,nindd/])
    end do
    delete([/ip,nt,nind/])
end do

llphi = llphi*10^(-6)

print("complete calc leadlag data. "+ systemfunc("date"))

;; resources
cmap = read_colormap_file("MPL_RdBu")
;shading u plot
rest = True
rest@gsnFrame           = False
rest@gsnDraw            = False
rest@gsnMaximize        = False
rest@gsnAddCyclic       = True
rest@gsnMajorLatSpacing = 15
rest@gsnMajorLonSpacing = 30
;rest@gsnCenterStringFontHeightF = 0.025 

rest@tmXTOn = False
rest@tmYROn = False

rest@mpFillOn        = False
rest@mpFillDrawOrder = "Postdraw"
;resolr@mpProjection= "Robinson"
rest@mpMinLatF   = -30
;rest@mpMaxLatF   = 40
rest@mpMinLonF   = 0
rest@mpCenterLonF= 180
rest@mpMaxLonF   = 360
rest@mpLandFillColor = "white"
rest@mpGridAndLimbDrawOrder = "PostDraw"
rest@mpGridLineColor        = "Gray60"
;rest@mpGridLineThicknessF   = 1.0

rest@cnFillOn               = True
;rest@cnFillMode             = "CellFill"
rest@cnFillPalette          = cmap(::-1,:)
;rest@cnFillOpacityF         = 0.6
rest@cnLinesOn              = False
rest@cnLineLabelsOn         = False
rest@cnLevelSelectionMode   = "Manuallevels"
rest@cnMinLevelValF         = -8.
rest@cnMaxLevelValF         =  8.
rest@cnLevelSpacingF        =  1.
rest@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
rest@lbLabelBarOn           = False

; contour phi plot.
resphi = True
resphi@gsnDraw      = False
resphi@gsnMaximize  = False
resphi@gsnFrame     = False
resphi@gsnAddCyclic = True
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

; contour olr plot.
resolr = True
resolr@gsnDraw      = False
resolr@gsnMaximize  = False
resolr@gsnFrame     = False
resolr@gsnAddCyclic = True
resolr@gsnContourNegLineDashPattern = 1
resolr@cnMonoLineColor    = True
resolr@cnLinesOn          = True
resolr@cnLineLabelsOn     = False
resolr@cnLineColor        = "firebrick4"
resolr@cnLineThicknessF   = 2.0
resolr@cnInfoLabelPerimOn = False
resolr@cnInfoLabelOn      = False
resolr@cnLevelSelectionMode = "Explicitlevels"
resolr@cnLevels      = (/-25,-20,-10,-5,5,10,20,25/)

;; vector res.
resw = True
resw@gsnFrame     = False
resw@gsnDraw      = False
resw@gsnMaximize  = False
resw@gsnAddCyclic = True
resw@vcGlyphStyle            = "LineArrow"
resw@vcLineArrowThicknessF   = 1.
;resw@vcMinDistanceF          = 0.01
resw@vcRefMagnitudeF         = 5.
resw@vcRefLengthF            = 0.05
resw@vcRefAnnoOrthogonalPosF = -1.0
resw@vcRefAnnoParallelPosF   = 1.18
resw@vcLineArrowColor        = "black"
resw@vcRefAnnoString2On      = False
resw@vcRefAnnoString1        = "$VMG$" +" m/s"
resw@vcRefAnnoPerimSpaceF    = 0.5
resw@vcRefAnnoPerimOn        = False
resw@vcRefAnnoFontHeightF    = 0.015
resw@vcRefAnnoSide           = "Left"
resw@vcRefAnnoJust           = "TopRight"
resw@vcRefAnnoBackgroundColor = -1

ressf = True
ressf@gsnFrame     = False
ressf@gsnDraw      = False
ressf@gsnMaximize  = False
ressf@gsnAddCyclic = True
ressf@stArrowLengthF     = 0.004                    ; size of the arrows.
ressf@stMinArrowSpacingF = 0.004                    ; arrow spacing.
ressf@stMinDistanceF     = 0.01
ressf@stArrowStride      = 3                        ; arrows start every third

resP = True
resP@gsnMaximize      = True
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
;resP@lbLabelStride     = 2
;resP@lbFillOpacityF    = 0.6
resP@gsnPanelRowSpec    = True
resP@lbLabelFontHeightF = 0.012
print("complete res. "+ systemfunc("date"))

;; plot.
i  = 1
ip = 5
wks = gsn_open_wks("x11", pdir+"/"+"t-shading_phi_contour_olr_contour_phase"+ip)

plott   = new(8,graphic)
plotphi = plott
plotolr = plott
plotvc  = plott

Label1 = (/"(a) phase 1 lag -15","(e) phase 1 lag 5",\
                              "(b) phase 1 lag -10","(f) phase 1 lag 10",\
                              "(c) phase 1 lag -5","(g) phase 1 lag 15",\
                              "(d) phase 1 lag 0","(h) phase 1 lag 20"/)
do n = 0,7
  x = llu(4,n,:,:)
  y = llphi(4,n,:,:)
  x = mask(x,x.le.m_valueu(0,i,n,:,:).or.x.ge.m_valueu(1,i,n,:,:),True)
  ;y = mask(y,y.le.m_valuephi(0,i,n,:,:).or.y.ge.m_valuephi(1,i,n,:,:),True)
  plott(n)   = gsn_csm_contour_map(wks, x, rest)
  plotphi(n) = gsn_csm_contour(wks, y, resphi)
  plotolr(n) = gsn_csm_contour(wks, llolr(4,n,:,:), resolr)
  plotu = gsn_csm_contour(wks, uc, resc)
  overlay(plott(n),plotu)
  overlay(plott(n),plotphi(n))
  overlay(plott(n),plotolr(n))
end do
gsn_panel(wks, (/plott(0),plott(4),plott(1),plott(5),\
                plott(2),plott(6),plott(3),plott(7)/), (/2,2,2,2/), resP)

i = 1
do while (i.eq.1)
    i = 1
end do
delete([/plott,plotphi,plotolr,wks/])

Label1 = (/"(a) phase 1 lag -15","(e) phase 1 lag 5",\
                              "(b) phase 1 lag -10","(f) phase 1 lag 10",\
                              "(c) phase 1 lag -5","(g) phase 1 lag 15",\
                              "(d) phase 1 lag 0","(h) phase 1 lag 20"/)
Label2 = (/"(a) phase 2 lag -15","(e) phase 2 lag 5",\
                              "(b) phase 2 lag -10","(f) phase 2 lag 10",\
                              "(c) phase 2 lag -5","(g) phase 2 lag 15",\
                              "(d) phase 2 lag 0","(h) phase 2 lag 20"/)
Label3 = (/"(a) phase 3 lag -15","(e) phase 3 lag 5",\
                              "(b) phase 3 lag -10","(f) phase 3 lag 10",\
                              "(c) phase 3 lag -5","(g) phase 3 lag 15",\
                              "(d) phase 3 lag 0","(h) phase 3 lag 20"/)
Label4 = (/"(a) phase 4 lag -15","(e) phase 4 lag 5",\
                              "(b) phase 4 lag -10","(f) phase 4 lag 10",\
                              "(c) phase 4 lag -5","(g) phase 4 lag 15",\
                              "(d) phase 4 lag 0","(h) phase 4 lag 20"/)
Label5 = (/"(a) phase 5 lag -15","(e) phase 5 lag 5",\
                              "(b) phase 5 lag -10","(f) phase 5 lag 10",\
                              "(c) phase 5 lag -5","(g) phase 5 lag 15",\
                              "(d) phase 5 lag 0","(h) phase 5 lag 20"/)
Label6 = (/"(a) phase 6 lag -15","(e) phase 6 lag 5",\
                              "(b) phase 6 lag -10","(f) phase 6 lag 10",\
                              "(c) phase 6 lag -5","(g) phase 6 lag 15",\
                              "(d) phase 6 lag 0","(h) phase 6 lag 20"/)
Label7 = (/"(a) phase 7 lag -15","(e) phase 7 lag 5",\
                              "(b) phase 7 lag -10","(f) phase 7 lag 10",\
                              "(c) phase 7 lag -5","(g) phase 7 lag 15",\
                              "(d) phase 7 lag 0","(h) phase 7 lag 20"/)
Label8 = (/"(a) phase 8 lag -15","(e) phase 8 lag 5",\
                              "(b) phase 8 lag -10","(f) phase 8 lag 10",\
                              "(c) phase 8 lag -5","(g) phase 8 lag 15",\
                              "(d) phase 8 lag 0","(h) phase 8 lag 20"/)

Label = new((/8,8/),string)
Label(0,:)  = Label1
Label(1,:)  = Label5
;Label(2,:) = Label3
;Label(3,:) = Label4
;Label(4,:) = Label5
;Label(5,:) = Label6
;Label(6,:) = Label7
;Label(7,:) = Label8
ph = (/1,5/)
print("Prepare labels done. "+ systemfunc("date"))

do i = 0, 1
    ip = ph(i)
    wks = gsn_open_wks("eps", pdir+"/"+"u-shading_phi_contour_olr_contour_uv_vector_phase"+ip)

    plott   = new(8,graphic)
    plotphi = plott
    plotolr = plott
    plotvc  = plott

    do n = 0,7
        x = llu(i,n,:,:)
        y = llphi(i,n,:,:)
        x = mask(x,x.le.m_valueu(0,i,n,:,:).or.x.ge.m_valueu(1,i,n,:,:),True)
        ;y = mask(y,y.le.m_valueu(0,i,n,:,:).or.y.ge.m_valueu(1,i,n,:,:),True)
        plott(n)   = gsn_csm_contour_map(wks, x, rest)
        plotphi(n) = gsn_csm_contour(wks, y, resphi)
        plotolr(n) = gsn_csm_contour(wks, llolr(i,n,:,:), resolr)
        plotu = gsn_csm_contour(wks, uc, resc)
        overlay(plott(n),plotu)
        overlay(plott(n),plotphi(n))
        overlay(plott(n),plotolr(n))
    end do
;;
    resP@gsnPanelFigureStringsPerimOn = True
    resP@gsnPanelFigureStringsJust    = "Bottom Right"
    resP@gsnPanelFigureStringsFontHeightF = 0.012
    resP@gsnPanelFigureStrings = Label(i,:)
    gsn_panel(wks, (/plott(0),plott(4),plott(1),plott(5),\
                plott(2),plott(6),plott(3),plott(7)/), (/2,2,2,2/), resP)
    delete([/plott,plotphi,plotolr,plotvc/])
    print("Draw Phase"+ip+" done. "+ systemfunc("date"))
end do
