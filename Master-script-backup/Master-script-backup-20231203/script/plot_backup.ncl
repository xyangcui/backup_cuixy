load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
;-----------------------------------------------------------------------------
;This script is used to plot phase-lag composite of 3 jet indices.
;     |
;     |
;phase|
;     |
;     |
;     ----------------------------------
;                   lag days
;reference: Jenney et al. (2019) JGR:A
;Information
;author: cui xiangyang
;time: 2023-10-21
;release: v1
;-----------------------------------------------------------------------------
datadir = getenv("outdir")
var1    = getenv("var1")
var2    = getenv("var2")

udir    = getenv("udir")

rmmdir  = getenv("rmmdir")
plotdir = getenv("plotdir")

latmin  = tofloat(getenv("latmin"))
latmax  = tofloat(getenv("latmax"))
lonmin  = tofloat(getenv("lonmin"))
lonmax  = tofloat(getenv("lonmax"))

nw    = tofloat(getenv("n1"))   ; nwa
nb    = tofloat(getenv("n2"))   ; nwb

clat  = tofloat(getenv("clat"))
clon  = tofloat(getenv("clon"))
;----------------------------------------------------------------------------
; draw u250.
;----------------------------------------------------------------------------
f    = addfile(udir,"r")

time = f->time
TIME = cd_calendar(time, 0)
mon  = toint(TIME(:,1))
nt   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)

u    = dim_avg_n_Wrap(f->u(nt,0,:,:),0)
; u250 resources.
;--plot
res1 = True
res1@gsnFrame     = False
res1@gsnDraw      = False
res1@gsnMaximize  = False
res1@gsnAddCyclic = True
res1@gsnMajorLatSpacing = 15
res1@gsnMajorLonSpacing = 30

res1@gsnRightString = "("+abs(clat)+"~S~o~N~N,"+clon+"~S~o~N~E)"

res1@gsnLeftStringFontHeightF = 0.012
res1@gsnRightStringFontHeightF = 0.012
res1@tmXBLabelFontHeightF  = 0.012
res1@tmYLLabelFontHeightF  = 0.012

res1@vpXF     = 0.15
res1@vpYF     = 0.95
res1@vpWidthF = 0.7
res1@vpHeightF= 0.2

res1@mpFillOn = True
res1@mpLandFillColor = "white"
res1@mpGridLineColor       = "gray"
res1@mpGridLineDashPattern = 1
res1@mpMinLatF     = latmin - 15
res1@mpMaxLatF     = latmax + 15
res1@mpMinLonF     = lonmin - 30
res1@mpMaxLonF     = lonmax + 30

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
res1@lbLabelStride         = 10
res1@lbTitleString         = "~F21~(unit: m/s)"
res1@lbTitlePosition       = "Bottom"
res1@lbTitleFontHeightF    = 0.01
res1@lbOrientation         = "Vertical"
res1@pmLabelBarHeightF     = 0.2

res1@gsnLeftString   = "~F22~a ~F21~Zonal wind at 250 hPa"

wks = gsn_open_wks("eps",plotdir)

plot  = gsn_csm_contour_map(wks, u, res1)

respoly = True
respoly@gsLineColor      = "mediumvioletred"
respoly@gsLineThicknessF = 2.0
polyx = (/lonmin,lonmax,lonmax,lonmin,lonmin/)
polyy = (/latmin,latmin,latmax,latmax,latmin/)
dum  = gsn_add_polyline(wks, plot, polyx, polyy, respoly)

draw(plot)
delete([/u,time,TIME,mon,nt,cmap1,cmap2/])

;-----------------------------------------------------------------------------
; Function to filter.
;-----------------------------------------------------------------------------
undef("ISO_Filter")
function ISO_Filter(x,time,npt,nwa,nwb)
local time
begin
  TIME1 = cd_calendar(time, 0)
  year1 = toint(TIME1(:,0))
  mon1  = toint(TIME1(:,1))
  day1  = toint(TIME1(:,2))
  ddd   = day_of_year(year1, mon1, day1)
  yyyyddd = year1*1000 + ddd
  wts     = filwgts_lanczos (npt, 2, 1./nwa,1./nwb , 1.0)
  xtmp  = new((/dimsizes(x),2,2/),float)
  xtmp  = 0.0
  xtmp(:,0,0) = x
  xClim = clmDayTLL(xtmp, yyyyddd)
  xAnom = calcDayAnomTLL(xtmp, yyyyddd, xClim)
  xr    = xAnom(:,0,0)
  xr := wgt_runave_Wrap(xr, wts, 0)

  return(xr)
  delete([/TIME1,year1,mon1,day1,xtmp,xAnom,xr/])

end

;-----------------------------------------------------------------------------
; Function to deseanalized.
;-----------------------------------------------------------------------------
undef("deseason")
function deseason(x,time)
local time
begin
  TIME1 = cd_calendar(time, 0)
  year1 = toint(TIME1(:,0))
  mon1  = toint(TIME1(:,1))
  day1  = toint(TIME1(:,2))
  ddd   = day_of_year(year1, mon1, day1)
  yyyyddd = year1*1000 + ddd

  xtmp  = new((/dimsizes(x),2,2/),float)
  xtmp  = 0.0
  xtmp(:,0,0) = x
  xClim = clmDayTLL(xtmp, yyyyddd)
  xAnom = calcDayAnomTLL(xtmp, yyyyddd, xClim)
  xr    = xAnom(:,0,0)

  return(xr)
  delete([/TIME1,year1,mon1,day1,xtmp,xAnom,xr/])

end

;-----------------------------------------------------------------------------
; Function to draw polygons in a viewpoint space.
;-----------------------------------------------------------------------------
undef("add_polygon")
function add_polygon(wks:graphic,plot:graphic,values[*][*],ranges,color)
local plot,values,xmin,xmax,ymin,ymax,nx,ny,fil_color
begin
    dims    = dimsizes(values)
    ny      = dims(1)
    nx      = dims(0)
;loop throw and draw polygons.
    getvalues plot
      "trXMinF" : xmin
      "trXMaxF" : xmax
      "trYMinF" : ymin
      "trYMaxF" : ymax
    end getvalues

    xsize  = 1.
    ysize  = 1.

    gnres                 = True
    gnres@tfPolyDrawOrder = "PostDraw"  
; First row.
    j = 0
    ypos = ymin
    do i=1,nx-2
        fil_color = get_color_rgba(color, ranges, values(i,j))
        gnres@gsFillColor = fil_color;fil_color(i,j,:)
        xpos = xmin + (xsize/2. * (2*i-1))
        xbox = (/xpos,xpos+xsize,xpos+xsize,xpos,xpos/)
        ybox = (/ypos,ypos,ypos+ysize/2.,ypos+ysize/2.,ypos/)    
        str = unique_string("fill_lft")     
        plot@$str$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)   
    end do
; Left
    i = 0
    fil_color = get_color_rgba(color, ranges, values(i,j))
    gnres@gsFillColor = fil_color;fil_color(i,j,:)
    ;gnres@gsFillColor = fil_color(i,j,:)
    xpos = xmin 
    xbox = (/xpos,xpos+xsize/2.,xpos+xsize/2.,xpos,xpos/)
    ybox = (/ypos,ypos,ypos+ysize/2.,ypos+ysize/2.,ypos/)
    str = unique_string("fill_lft")     
    plot@$str$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)
; The rightmost.
    i = nx-1
    fil_color = get_color_rgba(color, ranges, values(i,j))
    gnres@gsFillColor = fil_color;fil_color(i,j,:)
    ;gnres@gsFillColor = fil_color(i,j,:)
    xpos = xmin + (xsize/2. * (2*i-1))
    xbox = (/xpos,xmax,xmax,xpos,xpos/)
    ybox = (/ypos,ypos,ypos+ysize/2.,ypos+ysize/2.,ypos/)    
    str = unique_string("fill_lft")     
    plot@$str$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)

; Top
    j = ny - 1
    ypos = ymax - ysize/2.
    do i=1,nx-2
        fil_color = get_color_rgba(color, ranges, values(i,j))
        gnres@gsFillColor = fil_color;fil_color(i,j,:)
        ;gnres@gsFillColor = fil_color(i,j,:)
        xpos = xmin + (xsize/2. * (2*i-1))
        xbox = (/xpos,xpos+xsize,xpos+xsize,xpos,xpos/)
        ybox = (/ypos,ypos,ymax,ymax,ypos/)    
        str = unique_string("fill_lft")     
        plot@$str$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)    
    end do
; Left
    i = 0
    fil_color = get_color_rgba(color, ranges, values(i,j))
    gnres@gsFillColor = fil_color;fil_color(i,j,:)
    ;gnres@gsFillColor = fil_color(i,j,:)
    xpos = xmin 
    xbox = (/xpos,xpos+xsize/2.,xpos+xsize/2.,xpos,xpos/)
    ybox = (/ypos,ypos,ymax,ymax,ypos/)    
    str = unique_string("fill_lft")     
    plot@$str$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)
; The rightmost.
    i = nx-1
    fil_color = get_color_rgba(color, ranges, values(i,j))
    gnres@gsFillColor = fil_color;fil_color(i,j,:)
    ;gnres@gsFillColor = fil_color(i,j,:)
    xpos = xmin + (xsize/2. * (2*i-1))
    xbox = (/xpos,xmax,xmax,xpos,xpos/)
    ybox = (/ypos,ypos,ymax,ymax,ypos/)    
    str = unique_string("fill_lft")     
    plot@$str$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)
;Left column
    i = 0
    xpos = xmin
    do j = 1,ny-2
        fil_color = get_color_rgba(color, ranges, values(i,j))
        gnres@gsFillColor = fil_color;fil_color(i,j,:)
        ypos = ymin + (ysize/2.*(2*j-1))
        ;gnres@gsFillColor = fil_color(i,j,:)
        xbox = (/xpos,xpos+xsize/2.,xpos+xsize/2.,xpos,xpos/)
        ybox = (/ypos,ypos,ypos+ysize,ypos+ysize,ypos/)    
        str = unique_string("fill_lft")     
        plot@$str$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)  
    end do
; Right column
    i = nx-1
    xpos = xmax - xsize/2.
    do j = 1,ny-2
        fil_color = get_color_rgba(color, ranges, values(i,j))
        gnres@gsFillColor = fil_color;fil_color(i,j,:)
        ypos = ymin + (ysize/2.*(2*j-1))
        ;gnres@gsFillColor = fil_color(i,j,:)
        xbox = (/xpos,xmax,xmax,xpos,xpos/)
        ybox = (/ypos,ypos,ypos+ysize,ypos+ysize,ypos/)    
        str = unique_string("fill_lft")     
        plot@$str$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)  
    end do
;Middle values.
    do j=1,ny-2
        ypos = ymin + (ysize/2.*(2*j-1))
        do i=1,nx-2
            fil_color = get_color_rgba(color, ranges, values(i,j))
            gnres@gsFillColor = fil_color;fil_color(i,j,:)
            ;gnres@gsFillColor = fil_color(i,j,:)
            xpos = xmin + (xsize/2.*(2*i-1))
            xbox = (/xpos,xpos+xsize,xpos+xsize,xpos,xpos/)
            ybox = (/ypos,ypos,ypos+ysize,ypos+ysize,ypos/)    
            str = unique_string("fill_lft")     
            plot@$str$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)  
        end do
    end do

    return(plot)
end

;-----------------------------------------------------------------------------
;Procedure to attach a vertical lavelbar to the right of the plot.
;-----------------------------------------------------------------------------
undef("add_labelbar")
procedure add_labelbar(wks,plot,labels,color)
local vph,vpw,nboes,lbres,lbid,amres,annoid,color
begin
    getvalues plot
      "vpHeightF" :vph
      "vpWidthF"  :vpw
    end getvalues

    lbres = True

    lbres@vpWidthF           = 0.45*vpw
    lbres@vpHeightF          = 1.0*vph

    lbres@lbAutoManage = False
    lbres@lbFillColors = color
    lbres@lbMonoFillPattern  = True
    lbres@lbOrientation      = "Vertical"
    lbres@lbPerimOn          = False
    lbres@lbTitleOn          = False
    lbres@lbTitleFontHeightF = 0.02
    lbres@lbTitleDirection   = "Across"
    lbres@lbBoxMinorExtentF  = 0.2
    lbres@lbLabelOffsetF     = 0.02
    lbres@lbLabelAlignment   = "InteriorEdges"
    lbres@lbBoxLinesOn       = True
    lbres@lbFillOpacityF     = 1.0
    lbres@lbLabelStride      = 2
    lbres@lbBoxSeparatorLinesOn = False
    lbres@lbLabelFontHeightF = 0.01
    lbres@lbLabelAngleF      = -45
    lbid = gsn_create_labelbar(wks,dimsizes(color(:,0)),labels,lbres)
; attach labelbar to the plot.
    amres = True
    amres@amJust = "CenterLeft"
    amres@amParallelPosF   = 0.4
    amres@amOrthogonalPosF = 0.0
    annoid = gsn_add_annotation(plot,lbid,amres)
end 

;-----------------------------------------------------------------------------
; Main code.
;-----------------------------------------------------------------------------
;read in data and selected strong MJO days.

; jet indices plot.
f = addfile(datadir,"r")

time = f->time
TIME = cd_calendar(time, 0)
year = toint( TIME(:,0) )
mon  = toint( TIME(:,1) )
day  = toint( TIME(:,2) )

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

reamp   = new((/dimsizes(time)/),float)
rephase = new((/dimsizes(time)/),integer)
reamp   = 0.
rephase = 0

do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
    reamp(i)   = amp(nt)
    rephase(i) = phase(nt)
    delete(nt)
end do

jets1 = f->$var1$ ; jet speed
jetl1 = f->$var2$ ; jet latitude

jetl = ISO_Filter(jetl1,time,201,nw,nb)
jets = ISO_Filter(jets1,time,201,nw,nb)
;jetl  = deseason(jetl1,time)
;jets  = deseason(jets1,time)

;-----------------------------------------------------------------------------
;calculate phase-lag data. lag time 50 days
n = 30
plag   = new((/2,n+1,13/),float) ; (indices,lag days, phases)
plag   = 0.
plag1  = plag    ; amp = 1.5
plag2  = plag    ; amp = 1.0

;nd   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
;reamp    := reamp(nd)
;ndd  = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
;rephase  := rephase(ndd)
;in1  = ispan(0,dimsizes(time)-1,1)
;in2  = in1(nd)

;amp > 2.0
;am = 2.0
;do i = 1,8
;    nt  = ind(rephase.eq.i.and.reamp.gt.am)
;    nt1 = in2(nt)
;    do j = 0 ,n, 1
;        plag(0,j,i-1) = dim_avg(jets(nt1+j))
;        plag(1,j,i-1) = dim_avg(jetl(nt1+j))
;    end do
;    delete([/nt,nt1/])
;end do

;do i = 0,4
;    plag(:,:,i+8) = plag(:,:,i)
;end do

; amp > 1.5
;am = 1.5
;do i = 1,8
   ; nt  = ind(rephase.eq.i.and.reamp.gt.am)
  ;  nt1 = in2(nt)
 ;   do j = 0 ,n, 1
 ;       plag1(0,j,i-1) = dim_avg(jets(nt1+j))
 ;       plag1(1,j,i-1) = dim_avg(jetl(nt1+j))
;    end do
;    delete([/nt,nt1/])
;end do

;do i = 0,4
 ;   plag1(:,:,i+8) = plag1(:,:,i)
;end do

;amp>1.0
;am = 2.0
;do i = 1,8
 ;   nt  = ind(rephase.eq.i.and.reamp.gt.am)
 ;   nt1 = in2(nt)
;    do j = 0 ,n, 1
;        plag2(0,j,i-1) = dim_avg(jets(nt1+j))
;        plag2(1,j,i-1) = dim_avg(jetl(nt1+j))
;    end do
;    delete([/nt,nt1/])
;end do

;do i = 0,4
;    plag2(:,:,i+8) = plag2(:,:,i)
;end do

;first conver from 12 to 3.
n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.31)
reamp   := reamp(n1:n2)
rephase := rephase(n1:n2)
jetl    := jetl(n1:n2)
jets    := jets(n1:n2)
mon     := mon(n1:n2)
time    := time(n1:n2)

nd   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
reamp    := reamp(nd)
ndd  = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
rephase  := rephase(ndd)
in1  = ispan(0,dimsizes(time)-1,1)
in2  = in1(nd)

;amp > 2.0
am = 2.0
do i = 1,8
    nt  = ind(rephase.eq.i.and.reamp.gt.am)
    nt1 = in2(nt)
    do j = 0 ,n, 1
        plag(0,j,i-1) = dim_avg(jets(nt1+j))
        plag(1,j,i-1) = dim_avg(jetl(nt1+j))
    end do
    delete([/nt,nt1/])
end do

do i = 0,4
    plag(:,:,i+8) = plag(:,:,i)
end do

; amp > 1.5
am = 1.5
do i = 1,8
    nt  = ind(rephase.eq.i.and.reamp.gt.am)
    nt1 = in2(nt)
    do j = 0 ,n, 1
        plag1(0,j,i-1) = dim_avg(jets(nt1+j))
        plag1(1,j,i-1) = dim_avg(jetl(nt1+j))
    end do
    delete([/nt,nt1/])
end do

do i = 0,4
    plag1(:,:,i+8) = plag1(:,:,i)
end do

;amp>1.0
am = 1.0
do i = 1,8
    nt  = ind(rephase.eq.i.and.reamp.gt.am)
    nt1 = in2(nt)
    do j = 0 ,n, 1
        plag2(0,j,i-1) = dim_avg(jets(nt1+j))
        plag2(1,j,i-1) = dim_avg(jetl(nt1+j))
    end do
    delete([/nt,nt1/])
end do

do i = 0,4
    plag2(:,:,i+8) = plag2(:,:,i)
end do

print("calculate data done")
;Norm
; amp>2.0
;plagt =plag
;do i = 0,1
;    pmin = min(plag(i,:,:))
;    pmax = max(plag(i,:,:))
;    plagt(i,:,:) = 2*(plag(i,:,:)-pmin)/(pmax-pmin) - 1
;end do
;plag = plagt
;delete(plagt)
; amp>1.5
;plagt =plag1
;do i = 0,1
;    pmin = min(plag1(i,:,:))
;    pmax = max(plag1(i,:,:))
;    plagt(i,:,:) = 2*(plag1(i,:,:)-pmin)/(pmax-pmin) - 1
;end do
;plag1 = plagt
;delete(plagt)
;
; amp>1.0
;plagt =plag2
;do i = 0,1
;    pmin = min(plag2(i,:,:))
;    pmax = max(plag2(i,:,:))
;    plagt(i,:,:) = 2*(plag2(i,:,:)-pmin)/(pmax-pmin) - 1
;end do
;plag2 = plagt
;delete(plagt)

x = ispan(0,n,1)  ; x labels
y = (/1,2,3,4,5,6,7,8,1,2,3,4,5/)

color1 = "MPL_PRGn"
color2 = "MPL_PuOr"
;label1
be = floor(max(plag(0,:,:)))
en = ceil(min(plag(0,:,:)))
ll = (/be,abs(en)/)
a  = max(ll)
if (be.eq.abs(en).and.a.eq.0) then
    ncolor = toint((1 + 1)/0.1 + 2)
    cmap1  = span_color_rgba(color1, ncolor)
    label1 = fspan(-1., 1., ncolor-1)
else
    ncolor = toint((a + a)/0.2 + 2)
    cmap1  = span_color_rgba(color1, ncolor)
    label1 = fspan(-a, a, ncolor-1)
end if
;label2
be = floor(max(plag(1,:,:)))
en = ceil(min(plag(1,:,:)))
ll = (/be,abs(en)/)
a  = max(ll)
if (be.eq.abs(en).and.a.eq.0) then
    ncolor = toint((1 + 1)/0.1 + 2)
    cmap2  = span_color_rgba(color2, ncolor)
    label2 = fspan(-1., 1., ncolor-1)
else
    ncolor = toint((a + a)/0.2 + 2)
    cmap2  = span_color_rgba(color2, ncolor)
    label2 = fspan(-a, a, ncolor-1)
end if
;label3
be = floor(max(plag1(0,:,:)))
en = ceil(min(plag1(0,:,:)))
ll = (/be,abs(en)/)
a  = max(ll)
if (be.eq.abs(en).and.a.eq.0) then
    ncolor = toint((1 + 1)/0.1 + 2)
    cmap3  = span_color_rgba(color1, ncolor)
    label3 = fspan(-1., 1., ncolor-1)
else
    ncolor = toint((a + a)/0.2 + 2)
    cmap3  = span_color_rgba(color1, ncolor)
    label3 = fspan(-a, a, ncolor-1)
end if
;label4
be = floor(max(plag1(1,:,:)))
en = ceil(min(plag1(1,:,:)))
ll = (/be,abs(en)/)
a  = max(ll)
if (be.eq.abs(en).and.a.eq.0) then
    ncolor = toint((1 + 1)/0.1 + 2)
    cmap4  = span_color_rgba(color2, ncolor)
    label4 = fspan(-1., 1., ncolor-1)
else
    ncolor = toint((a + a)/0.2 + 2)
    cmap4  = span_color_rgba(color2, ncolor)
    label4 = fspan(-a, a, ncolor-1)
end if
;label5
be = floor(max(plag2(0,:,:)))
en = ceil(min(plag2(0,:,:)))
ll = (/be,abs(en)/)
a  = max(ll)
if (be.eq.abs(en).and.a.eq.0) then
    ncolor = toint((1 + 1)/0.1 + 2)
    cmap5  = span_color_rgba(color1, ncolor)
    label5 = fspan(-1., 1., ncolor-1)
else
    ncolor = toint((a + a)/0.2 + 2)
    cmap5  = span_color_rgba(color1, ncolor)
    label5 = fspan(-a, a, ncolor-1)
end if
;label6
be = floor(max(plag2(1,:,:)))
en = ceil(min(plag2(1,:,:)))
ll = (/be,abs(en)/)
a  = max(ll)
if (be.eq.abs(en).and.a.eq.0) then
    ncolor = toint((1 + 1)/0.1 + 2)
    cmap6  = span_color_rgba(color2, ncolor)
    label6 = fspan(-1., 1., ncolor-1)
else
    ncolor = toint((a + a)/0.2 + 2)
    cmap6  = span_color_rgba(color2, ncolor)
    label6 = fspan(-a, a, ncolor-1)
end if
print("calculate labels done")
;---Resources for blank plot  
bres                = True                    
bres@gsnDraw        = False
bres@gsnFrame       = False
bres@gsnMaximize    = False     
bres@gsnLeftStringFontHeightF  = 0.012
bres@gsnRightStringFontHeightF = 0.012

bres@vpWidthF       = 0.36
bres@vpHeightF      = 0.16

bres@trYMaxF  = 12
bres@trYMinF  = 0
bres@tmYLLabelFontHeightF = 0.012
bres@tmYLMinorPerMajor    = 0
bres@tmYLTickSpacingF     = 1
bres@tmYLMajorLengthF     = 0.01
bres@tmYLMode             = "Explicit"
bres@tmYLValues           = (/0,1,2,3,4,5,6,7,8,9,10,11,12/)
bres@tmYLLabels           = y 
bres@tiXAxisOffsetYF      = 0.01
bres@tiYAxisFontHeightF   = 0.012

bres@tmYROn = False
bres@tmXTOn = False

bres@trXMaxF  = n
bres@trXMinF  = 0

bres@tmXBMinorPerMajor     = 8
bres@tmXBLabelFontHeightF  = 0.012
bres@tmXBMajorLengthF      = 0.01
bres@tmXBTickSpacingF      = 5.
bres@tmXBMode              = "Manual"
bres@tiXAxisFontHeightF    = 0.012

; amp > 2.0
bres@gsnRightString = "~F21~Amp>2"
plot1 = new(2, graphic)
k = 0
bres@vpXF  = 0.07
bres@vpYF  = 0.69
bres@tiYAxisString         = "~F21~MJO phase"
bres@tiXAxisString         = "" 
bres@gsnLeftString = "~F22~b ~F21~Jet speed"
plot1(k) = gsn_csm_blank_plot(wks,bres)
plot1(k) = add_polygon(wks,plot1(k),plag(k,:,:),label1,color1)
add_labelbar(wks,plot1(k),label1,cmap1)
print("plot b finished")

k = 1
bres@vpXF  = 0.55
bres@vpYF  = 0.69
bres@tiYAxisString         = ""
bres@tiXAxisString         = "" 
bres@gsnLeftString = "~F22~e ~F21~Jet latitude"
plot1(k) = gsn_csm_blank_plot(wks,bres)
plot1(k) = add_polygon(wks,plot1(k),plag(k,:,:),label2,color2)
add_labelbar(wks,plot1(k),label2,cmap2)
print("plot e finished")

;amp >1.5
bres@gsnRightString = "~F21~Amp>1.5"
plot2 = new(2, graphic)
k = 0
bres@vpXF  = 0.07
bres@vpYF  = 0.46
bres@tiYAxisString         = "~F21~MJO phase"
bres@tiXAxisString         = "" 
bres@gsnLeftString = "~F22~c ~F21~Jet speed"
plot2(k) = gsn_csm_blank_plot(wks,bres)
plot2(k) = add_polygon(wks,plot2(k),plag1(k,:,:),label3,color1)
add_labelbar(wks,plot2(k),label3,cmap3)
print("plot c finished")

k = 1
bres@vpXF  = 0.55
bres@vpYF  = 0.46
bres@tiYAxisString         = ""
bres@tiXAxisString         = "" 
bres@gsnLeftString = "~F22~f ~F21~Jet latitude"
plot2(k) = gsn_csm_blank_plot(wks,bres)
plot2(k) = add_polygon(wks,plot2(k),plag1(k,:,:),label4,color2)
add_labelbar(wks,plot2(k),label4,cmap4)
print("plot f finished")

;amp >1.0
bres@gsnRightString = "~F21~Amp>1"
plot3 = new(2, graphic)
k = 0
bres@vpXF  = 0.07
bres@vpYF  = 0.23
bres@tiYAxisString         = "~F21~MJO phase"
bres@tiXAxisString         = "~F21~lag days" 
bres@gsnLeftString = "~F22~d ~F21~Jet speed"
plot3(k) = gsn_csm_blank_plot(wks,bres)
plot3(k) = add_polygon(wks,plot3(k),plag2(k,:,:),label5,color1)
add_labelbar(wks,plot3(k),label5,cmap5)
print("plot d finished")

k = 1
bres@vpXF  = 0.55
bres@vpYF  = 0.23
bres@tiYAxisString         = ""
bres@tiXAxisString         = "~F21~lag days" 
bres@gsnLeftString = "~F22~g ~F21~Jet latitude"
plot3(k) = gsn_csm_blank_plot(wks,bres)
plot3(k) = add_polygon(wks,plot3(k),plag2(k,:,:),label6,color2)
add_labelbar(wks,plot3(k),label6,cmap6)
print("plot g finished")

do k = 0, 1
    draw(plot1(k))
    draw(plot2(k))
    draw(plot3(k))
end do
