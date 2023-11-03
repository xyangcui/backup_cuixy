;--------------------------------------------------
; plot lag days for every phase. u250(DJF,:,:)
; lag 12 days.
;--------------------------------------------------
INDIR   = getenv("outdir")
rmmdir  = getenv("rmmdir")
plotdir = getenv("plotdir")
var     = getenv("var")
l       = toint(getenv("i"))
latmin = tofloat(getenv("latmin"))
latmax = tofloat(getenv("latmax"))
lonmin = tofloat(getenv("lonmin"))
lonmax = tofloat(getenv("lonmax"))
;--------------------------------------------------
; Main code.
;--------------------------------------------------
;-- RMM phases.
f = addfile(INDIR,"r")

time = f->time
TIME = cd_calendar(time, 0)
year = toint( TIME(:,0) )
mon  = toint( TIME(:,1) )
day  = toint( TIME(:,2) )
unit = time@units

RMMf  = rmmdir
delim = ","
rmmf  = asciiread(RMMf,-1,"string")
iyear = str_get_field(rmmf,1,delim)
imon  = str_get_field(rmmf,2,delim)
iday  = str_get_field(rmmf,3,delim)
phase = str_get_field(rmmf,6,delim)
amp   = str_get_field(rmmf,7,delim)


iyear := toint(iyear(1::))
imon  := toint(imon(1::))
iday  := toint(iday(1::))
phase := toint(phase(1::))
amp   := tofloat(amp(1::))
reamp    = new((/dimsizes(time)/),float)
rephase1 = new((/dimsizes(time)/),integer)
reamp    = 0.
rephase1 = 0

do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
    reamp(i)   = amp(nt)
    rephase1(i) = phase(nt)
    delete(nt)
end do

x    = f->$var$
lat  = x&lat
lon  = x&lon
nlat = dimsizes(lat)
nlon = dimsizes(lon)
nd   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
x       := x(nd,:,:)
reamp   := reamp(nd)
ndd  = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
rephase  = rephase1(ndd)
;--lag values.
plag = new((/8,3,nlat,nlon/),float)
copy_VarCoords(x(0,:,:), plag(0,0,:,:))
do l = 1,8
  j = 0
    nt  = ind(reamp.gt.1.5.and.rephase.eq.l)
    nc  = nt + j
    print(dimsizes(nc))
    plag(l-1,j,:,:) = dim_avg_n(x(nc,:,:),0)
    delete([/nt,nc/])
  j = 1
    nt  = ind(rephase.eq.l.and.reamp.gt.1.5)
    nc  = nt + j
    plag(l-1,j,:,:) = dim_avg_n(x(nc,:,:),0)
    delete([/nt,nc/])
  j = 2
    nt  = ind(rephase.eq.l.and.reamp.gt.1.5)
    nc  = nt + j
    plag(l-1,j,:,:) = dim_avg_n(x(nc,:,:),0)
    delete([/nt,nc/])
end do

;--plot resources.
labels = "Lag "+ispan(0, 11, 1)+"day"
res = True
res@gsnDraw              = False        
res@gsnFrame             = False       
res@gsnMaximize          = False
res@gsnAddCyclic         = False
res@gsnMajorLatSpacing   = 10
res@gsnMajorLonSpacing   = 20

res@vpWidthF             = 0.8
res@vpHeightF            = 0.45 

res@mpFillOn             = True
res@mpMinLatF            = latmin         
res@mpMaxLatF            = latmax
res@mpMinLonF            = lonmin
res@mpMaxLonF            = lonmax    
res@mpLandFillColor      = "white"

res@cnFillOn             = True         
res@cnFillPalette        = "GMT_polar" 
res@cnLinesOn            = False        
res@cnLineLabelsOn       = False   
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF       =  -5.0        
res@cnMaxLevelValF       =  5.0    
res@cnLevelSpacingF      =  1.0    
res@lbLabelBarOn         = False        

res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025                                                  


;-- plot.
wks = gsn_open_wks("eps", plotdir)
plot = new(24, graphic)
do l = 0,7
  do k = 0, 2
    res@gsnLeftString = labels(k)
    plot(k+l*3) = gsn_csm_contour_map(wks, plag(l,k,:,:), res)
  end do
end do

resP = True     
resP@gsnDraw                  = True        
resP@gsnMaximize              = True 
resP@gsnFrame                 = False        
resP@gsnPanelLabelBar         = True 
resP@gsnPanelTop              = 0.95     
resP@gsnPanelBottom           = 0.03
resP@gsnPanelMainFontHeightF  = 0.015
resP@gsnPanelRowSpec          = True

resP@pmLabelBarOrthogonalPosF = -0.02         
resP@pmLabelBarWidthF         = 0.3 

gsn_panel(wks,plot,(/3,3,3,3,3,3,3,3/),resP) 
delete(wks)

