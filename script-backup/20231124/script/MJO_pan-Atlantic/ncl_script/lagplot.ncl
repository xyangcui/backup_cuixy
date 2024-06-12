;--------------------------------------------------
; plot lag days for every phase. u250(DJF,:,:)
; lag 12 days.
;--------------------------------------------------
INDIR   = getenv("outdir")
rmmdir  = getenv("rmmdir")
plotdir = getenv("plotdir")
var     = getenv("var")
am       = tofloat(getenv("i"))
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
reamp   := reamp(nd)
ndd  = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
rephase  = rephase1(ndd)

;indices
in1  = ispan(0,dimsizes(time)-1,1)
in2  = in1(nd)

;--lag values.
plag = new((/8,3,nlat,nlon/),float)
std1 = plag
n1   = plag(:,:,0,0)
copy_VarCoords(x(0,:,:), plag(0,0,:,:))

do l = 1,8
  j = 0
    nt  = ind(rephase.eq.l.and.reamp.gt.am)
    nt1 = in2(nt)
    nc  = nt1 + j
    plag(l-1,0,:,:) = dim_avg_n(x(nc,:,:),0)
    std1(l-1,0,:,:) = dim_stddev_n(x(nc,:,:),0)
    n1(l-1,0) = dimsizes(nt)
    delete([/nt,nc,nt1/])
  j = 5
    nt  = ind(rephase.eq.l.and.reamp.gt.am)
    nt1 = in2(nt)
    nc  = nt1 + j
    plag(l-1,1,:,:) = dim_avg_n(x(nc,:,:),0)
    std1(l-1,1,:,:) = dim_stddev_n(x(nc,:,:),0)
    n1(l-1,1) = dimsizes(nt)
    delete([/nt,nc,nt1/])
  j = 10
    nt  = ind(rephase.eq.l.and.reamp.gt.am)
    nt1 = in2(nt)
    nc  = nt1 + j
    plag(l-1,2,:,:) = dim_avg_n(x(nc,:,:),0)
    std1(l-1,2,:,:) = dim_stddev_n(x(nc,:,:),0)
    n1(l-1,2) = dimsizes(nt)
    delete([/nt,nc,nt1/])
end do

;--plot resources.

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


;--test
sig  = 0.05
ndi  = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2.or.mon.eq.3)
ave2 = dim_avg_n(x(ndi,:,:), 0)

res1 = True
res1@gsnDraw              = False       
res1@gsnAddCyclic         = False
res1@gsnFrame             = False
res1@cnLinesOn            = False       
res1@cnLineLabelsOn       = False
res1@cnFillScaleF         = 1.0        
res1@cnLevelSelectionMode = "ManualLevels" 
res1@cnMinLevelValF       = 0.0      
res1@cnMaxLevelValF       = 1.05     
res1@cnLevelSpacingF      = 0.05     
res1@cnInfoLabelPerimOn   = False
res1@cnInfoLabelOn        = False
res1@trGridType = "TriangularMesh"
res1@cnLevelSelectionMode = "ExplicitLevels"                                                                                             
res1@cnFillDrawOrder="postDraw"                                                                                                     
res1@cnLevels    = (/0,0.05/)                                                                                                
res1@cnMonoFillPattern    = False            ; want multiple patterns                                                               
res1@cnFillPatterns       = (/17,-1,17/)     ; the patterns                                                                        
res1@cnMonoFillScale      = False            ; want different densities                                                            
res1@cnFillScales         = (/1,1,1/)        ; change densities                                                                        
res1@cnMonoFillColor      = True                                                                                                     
res1@cnFillDotSizeF       = 0.002 



labels = "lag "+(/0,5,10/)+"day"
;-- plot.
wks    = gsn_open_wks("eps", plotdir)
plot   = new(24, graphic)
pshade = plot
do l = 0,7
  do k = 0, 2
    res@gsnLeftString = "Phase"+(l+1)+" "+labels(k)
    plot(k+l*3) = gsn_csm_contour_map(wks, plag(l,k,:,:), res)
    df   = n1(l,k) - 1
    t_ts = (plag(l,k,:,:)-ave2)/(std1(l,k,:,:)/sqrt(n1(l,k)))
    prob = student_t(t_ts, df)
    print(num(prob.lt.sig))
    print(min(prob))
    copy_VarCoords(x(0,:,:),prob)
    pshade(k+l*3) = gsn_csm_contour(wks,prob,res1)
    overlay(plot(k+l*3),pshade(k+l*3))
    delete([/prob,t_ts/])
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
resP@lbLabelFontHeightF       = 0.015
resP@pmLabelBarOrthogonalPosF = -0.02         
resP@pmLabelBarWidthF         = 0.3 

gsn_panel(wks,plot,(/3,3,3,3,3,3,3,3/),resP) 
delete(wks)
