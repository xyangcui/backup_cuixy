;--------------------------------------------------------
; draw temperature anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour.
; time: 2023-11-22
;--------------------------------------------------------

;--------------------------------------------------------
; merge MJO dates.
;--------------------------------------------------------
undef("setMJO")
procedure setMJO(Rmmdir,Time,Year,Mon,Day)
begin

  

end 
;-----------------------------------------------------------

ddir   = "/home/sunming/data5/cuixy/global_jets/data"
pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

phidir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022.nc"
tdir   = ddir+"/"+"ERA5_daily_anom_tvert_1979-2022_2.5x2.5.nc"
olrdir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
udir   = ddir+"/"+"ERA5_daily_u250_1979-2022.nc"


;--------------------------------------------------------
; u250
;--------------------------------------------------------
f    = addfile(udir,"r")

time = f->time
TIME = cd_calendar(time, 0)
mon  = toint(TIME(:,1))
nt   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)

u    = dim_avg_n_Wrap(f->u(nt,0,:,:),0)

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
resc@cnLineThicknessF   = 2.0
resc@cnInfoLabelPerimOn = False
resc@cnInfoLabelOn      = False
resc@cnLevelSelectionMode = "Explicitlevels"
resc@cnLevels      = (/30.,40.,50./)

delete([/time,f,TIME,mon,nt/])

print("draw clim.u complete. "+ systemfunc("date"))
;--------------------------------------------------------

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------

fphi = addfile(phidir,"r")
ft   = addfile(tdir,"r")
folr = addfile(olrdir,"r")

olr  = folr->olr(time|:,lat|:,lon|:)
phi  = fphi->phi
t    = ft->t

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
llphi = new((/8,8,dimsizes(phi&latitude),dimsizes(phi&longitude)/),float)
copy_VarCoords(phi(0,:,:), llphi(0,0,:,:))
llphi = 0.

llt = new((/8,8,dimsizes(t&lat),dimsizes(t&lon)/),float)
copy_VarCoords(t(0,:,:), llt(0,0,:,:))
llt = 0.

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
        llt(i,j,:,:)   = dim_avg_n(t(nindd,:,:), 0)
        llolr(i,j,:,:) = dim_avg_n(olr(nindd,:,:), 0)
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
rest@gsnAddCyclic       = True
rest@gsnMajorLatSpacing = 15
rest@gsnMajorLonSpacing = 30
;rest@gsnCenterStringFontHeightF = 0.025 

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
rest@cnMinLevelValF         = -1.5
rest@cnMaxLevelValF         =  1.5
rest@cnLevelSpacingF        =  0.3
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
resphi@cnLineThicknessF  = 3.0
resphi@cnInfoLabelPerimOn = False
resphi@cnInfoLabelOn      = True
resphi@cnLevelSelectionMode = "Explicitlevels"
resphi@cnLevels      = (/-13,-10,-7,-4,-1,1,4,7,10,13/)

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
resolr@cnLineColor        = "green4"
resolr@cnLineThicknessF   = 4.0
resolr@cnInfoLabelPerimOn = False
resolr@cnInfoLabelOn      = False
resolr@cnLevelSelectionMode = "Explicitlevels"
resolr@cnLevels      = (/-20,-10,10,20/)

resP = True
resP@gsnMaximize      = True
resP@gsnFrame         = False
resP@gsnDraw          = True
resP@gsnPanelLabelBar = True
resP@gsnPanelTop      = 0.98
resP@gsnPanelBottom   = 0.02

resP@pmLabelBarOrthogonalPosF = 0.0
resP@pmLabelBarWidthF         = 0.85

resP@lbTitleString     = "Temp anom (K)"
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
;i  = 0
;ip = i + 1
;wks = gsn_open_wks("x11", pdir+"/"+"t-shading_phi_contour_olr_contour_phase"+ip)
;
;plott   = new(8,graphic)
;plotphi = plott
;plotolr = plott

;do n = 0,7
 ;   plott(n)   = gsn_csm_contour_map(wks, llt(i,n,:,:), rest)
 ;   plotphi(n) = gsn_csm_contour(wks, llphi(i,n,:,:), resphi)
 ;   plotolr(n) = gsn_csm_contour(wks, llolr(i,n,:,:), resolr)
 ;   plotu = gsn_csm_contour(wks, u, resc)
 ;   overlay(plott(n),plotu)
;    overlay(plott(n),plotphi(n))
;    overlay(plott(n),plotolr(n))
;end do

;resP@gsnPanelFigureStringsPerimOn = True
;resP@gsnPanelFigureStringsJust    = "Bottom Right"
;resP@gsnPanelFigureStringsFontHeightF = 0.012
;resP@gsnPanelFigureStrings = (/"(a) phase 1 lag -15","(e) phase 1 lag 5",\
;                              "(b) phase 1 lag -10","(f) phase 1 lag 10",\
;                              "(c) phase 1 lag -5","(g) phase 1 lag 15",\
;                              "(d) phase 1 lag 0","(h) phase 1 lag 20"/)
;gsn_panel(wks, (/plott(0),plott(4),plott(1),plott(5),\
;                plott(2),plott(6),plott(3),plott(7)/), (/2,2,2,2/), resP)
;delete([/plott,plotphi,plotolr/])

;i = 1
;do while (i.eq.1)
    i = 1
;end do

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
Label(0,:) = Label1
Label(1,:) = Label2
Label(2,:) = Label3
Label(3,:) = Label4
Label(4,:) = Label5
Label(5,:) = Label6
Label(6,:) = Label7
Label(7,:) = Label8

print("Prepare labels done. "+ systemfunc("date"))

do i = 0, 7
    ip = i + 1
    wks = gsn_open_wks("eps", pdir+"/"+"t-shading_phi_contour_olr_contour_phase"+ip)

    plott   = new(8,graphic)
    plotphi = plott
    plotolr = plott

    do n = 0,7
        plott(n)   = gsn_csm_contour_map(wks, llt(i,n,:,:), rest)
        plotphi(n) = gsn_csm_contour(wks, llphi(i,n,:,:), resphi)
        plotolr(n) = gsn_csm_contour(wks, llolr(i,n,:,:), resolr)
        plotu = gsn_csm_contour(wks, u, resc)
        overlay(plott(n),plotu)
        ;overlay(plott(n),plotphi(n))
        overlay(plott(n),plotolr(n))
    end do
;;
    resP@gsnPanelFigureStringsPerimOn = True
    resP@gsnPanelFigureStringsJust    = "Bottom Right"
    resP@gsnPanelFigureStringsFontHeightF = 0.012
    resP@gsnPanelFigureStrings = Label(i,:)
    gsn_panel(wks, (/plott(0),plott(4),plott(1),plott(5),\
                plott(2),plott(6),plott(3),plott(7)/), (/2,2,2,2/), resP)
    delete([/plott,plotphi,plotolr/])
    print("Draw Phase"+ip+" done. "+ systemfunc("date"))
end do
