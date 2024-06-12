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
var3    = getenv("var3")

rmmdir  = getenv("rmmdir")
plotdir = getenv("plotdir")

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
    nrng    = dimsizes(ranges)
    ncolor  = dimsizes(color(:,0))
;---create colors.
   ; fil_color = new((/nx,ny,4/),float)
   ; fil_color = 0.
   ; fil_color(:,:,0) = where(values.le.ranges(1),color(0,0),fil_color(:,:,0))
   ; fil_color(:,:,1) = where(values.le.ranges(1),color(0,1),fil_color(:,:,1))
   ; fil_color(:,:,2) = where(values.le.ranges(1),color(0,2),fil_color(:,:,2))
  ;  fil_color(:,:,3) = where(values.le.ranges(1),color(0,3),fil_color(:,:,3))

  ;  fil_color(:,:,0) = where(values.gt.ranges(nrng-2),color(ncolor-1,0),fil_color(:,:,0))
  ;  fil_color(:,:,1) = where(values.gt.ranges(nrng-2),color(ncolor-1,1),fil_color(:,:,1))
   ; fil_color(:,:,2) = where(values.gt.ranges(nrng-2),color(ncolor-1,2),fil_color(:,:,2))
  ;  fil_color(:,:,3) = where(values.gt.ranges(nrng-2),color(ncolor-1,3),fil_color(:,:,3))

  ;  do i=1,nrng-3
  ;      rng_beg = ranges(i)
  ;      rng_end = ranges(i+1)
 ;       fil_color(:,:,0) = where(rng_beg.le.values.and.values.lt.rng_end,color(i,0),fil_color(:,:,0))
 ;       fil_color(:,:,1) = where(rng_beg.le.values.and.values.lt.rng_end,color(i,1),fil_color(:,:,1))
 ;       fil_color(:,:,2) = where(rng_beg.le.values.and.values.lt.rng_end,color(i,2),fil_color(:,:,2))
;        fil_color(:,:,3) = where(rng_beg.le.values.and.values.lt.rng_end,color(i,3),fil_color(:,:,3))
;    end do

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
    dum = new(dimsizes(values),graphic)
; First row.
    j = 0
    ypos = ymin
    do i=1,nx-2
        fil_color = get_color_rgba("MPL_PRGn", ranges, values(i,j))
        gnres@gsFillColor = fil_color;fil_color(i,j,:)
        xpos = xmin + (xsize/2. * (2*i-1))
        xbox = (/xpos,xpos+xsize,xpos+xsize,xpos,xpos/)
        ybox = (/ypos,ypos,ypos+ysize/2.,ypos+ysize/2.,ypos/)    
        str = unique_string("fill_lft")     
        plot@$str$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)   
    end do
; Left
    i = 0
    fil_color = get_color_rgba("MPL_PRGn", ranges, values(i,j))
    gnres@gsFillColor = fil_color;fil_color(i,j,:)
    ;gnres@gsFillColor = fil_color(i,j,:)
    xpos = xmin 
    xbox = (/xpos,xpos+xsize/2.,xpos+xsize/2.,xpos,xpos/)
    ybox = (/ypos,ypos,ypos+ysize/2.,ypos+ysize/2.,ypos/)
    str = unique_string("fill_lft")     
    plot@$str$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)
; The rightmost.
    i = nx-1
    fil_color = get_color_rgba("MPL_PRGn", ranges, values(i,j))
    gnres@gsFillColor = fil_color;fil_color(i,j,:)
    ;gnres@gsFillColor = fil_color(i,j,:)
    xpos = xmin + (xsize/2. * (2*i-1))
    xbox = (/xpos,xmax,xmax,xpos,xpos/)
    ybox = (/ypos,ypos,ypos+ysize/2.,ypos+ysize/2.,ypos/)    
    str = unique_string("fill_lft")     
    plot@$str$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)

; Top
    j = 0
    ypos = ymax - ysize/2.
    do i=1,nx-2
        fil_color = get_color_rgba("MPL_PRGn", ranges, values(i,j))
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
    fil_color = get_color_rgba("MPL_PRGn", ranges, values(i,j))
    gnres@gsFillColor = fil_color;fil_color(i,j,:)
    ;gnres@gsFillColor = fil_color(i,j,:)
    xpos = xmin 
    xbox = (/xpos,xpos+xsize/2.,xpos+xsize/2.,xpos,xpos/)
    ybox = (/ypos,ypos,ymax,ymax,ypos/)    
    str = unique_string("fill_lft")     
    plot@$str$ = gsn_add_polygon(wks,plot,xbox,ybox,gnres)
; The rightmost.
    i = nx-1
    fil_color = get_color_rgba("MPL_PRGn", ranges, values(i,j))
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
        fil_color = get_color_rgba("MPL_PRGn", ranges, values(i,j))
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
        fil_color = get_color_rgba("MPL_PRGn", ranges, values(i,j))
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
            fil_color = get_color_rgba("MPL_PRGn", ranges, values(i,j))
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

    lbres@vpWidthF           = 0.15*vpw
    lbres@vpHeightF          = 0.95*vph

    lbres@lbAutoManage = False
    lbres@lbFillColors = color
    lbres@lbMonoFillPattern  = True
    lbres@lbLabelFontHeightF = 0.02
    lbres@lbOrientation      = "Vertical"
    lbres@lbPerimOn          = False
    lbres@lbTitleOn          = False
    lbres@lbTitleFontHeightF = 0.02
    lbres@lbTitleDirection   = "Across"
    lbres@lbBoxMinorExtentF  = 0.2
    lbres@lbTitlePosition    = "Right"
    lbres@lbLabelOffsetF     = 0.05
    lbres@lbLabelAlignment   = "InteriorEdges"
    lbres@lbBoxLinesOn       = True
    lbres@lbFillOpacityF     = 1.0
    lbres@lbTitleString      = labels
    lbid = gsn_create_labelbar(wks,13,labels,lbres)
; attach labelbar to the plot.
    amres = True
    amres@amJust = "CenterLeft"
    amres@amParallelPosF   = 0.55
    amres@amOrthogonalPosF = 0.0
    annoid = gsn_add_annotation(plot,lbid,amres)
end 

;-----------------------------------------------------------------------------
; Main code.
;-----------------------------------------------------------------------------
;read in data and selected strong MJO days.
f = addfile(datadir,"r")

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
jetz1 = f->$var3$ ; jet zonality

jets2 = new((/dimsizes(jets1),1,1/),float)
jetl2 = new((/dimsizes(jets1),1,1/),float)
jetz2 = new((/dimsizes(jets1),1,1/),float)
jets2 = 0.
jetl2 = 0.
jetz2 = 0.
jets2(:,0,0) = jets1
jetl2(:,0,0) = jetl1
jetz2(:,0,0) = jetz1
copy_VarCoords(jets1, jets2(:,0,0))
copy_VarCoords(jets1, jetl2(:,0,0))
copy_VarCoords(jets1, jetz2(:,0,0))

jets3 = iso_filter(jets2,201,100,20,0)
jetl3 = iso_filter(jetl2,201,100,20,0)
jetz3 = iso_filter(jetz2,201,100,20,0)

jets  = jets3(:,0,0)
jetl  = jetl3(:,0,0)
jetz  = jetz3(:,0,0)

;-----------------------------------------------------------------------------
;calculate phase-lag data. lag time 50 days
n = 30
plag   = new((/3,n+1,16/),float) ; (indices,lag days, phases)
plag!0 = "indices"
plag!1 = "lagdays"
plag!2 = "phase"
plag&lagdays = ispan(0, n, 1)
plag&phase   = ispan(1, 16, 1)

nd   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
reamp    := reamp(nd)
ndd  = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
rephase  := rephase(ndd)

;indices
in1  = ispan(0,dimsizes(time)-1,1)
in2  = in1(nd)

;phase 1-8
do i = 1,8
    nt  = ind(rephase.eq.i.and.reamp.gt.1.)
    nt1 = in2(nt)
    print(dimsizes(nt))
    do j = 0 ,n, 1
        plag(0,j,i-1) = dim_avg(jets(nt1+j))
        plag(1,j,i-1) = dim_avg(jetl(nt1+j))
        plag(2,j,i-1) = dim_avg(jetz(nt1+j))
    end do
    delete([/nt,nt1/])
end do

do i = 1,8
    nt  = ind(rephase.eq.i.and.reamp.gt.1.)
    nt1 = in2(nt)
    print(dimsizes(nt))
    do j = 0 ,n, 1
        plag(0,j,i+7) = dim_avg(jets(nt1+j))
        plag(1,j,i+7) = dim_avg(jetl(nt1+j))
        plag(2,j,i+7) = dim_avg(jetz(nt1+j))
    end do
    delete([/nt,nt1/])
end do

plag := dim_standardize_n_Wrap(plag,1,(/1,2/))

;plot
label = new((/3,13/),float)
label = 0.
label(0,:) = fspan(-3.,3., 13)
label(1,:) = fspan(-3.,3., 13)
label(2,:) = fspan(-3.,3., 13)

;label(0,:) = fspan(-3.,3., 13)
;label(1,:) = fspan(-3.,3., 13)
;label(2,:) = fspan(-3.,3., 13)
Label = (/"Jet speed","Jet latitude","Jet zonality"/)
Title = (/"(a) jet speed","(b) jet latitude","(c) jet zonality"/)
x = ispan(0,n,1)  ; x labels
y = (/1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8/)


cmap  = span_color_rgba("MPL_PRGn", 13)


;---Resources for blank plot  
bres                = True                    
bres@gsnDraw        = False
bres@gsnFrame       = False
bres@gsnMaximize    = False     

bres@vpWidthF       = 0.7
bres@vpHeightF      = 0.3

bres@trYMaxF  = 15
bres@trYMinF  = 0
bres@tmYLLabelFontHeightF = 0.015
bres@tmYLMinorPerMajor    = 0
bres@tmYLTickSpacingF     = 1
bres@tmYLMajorLengthF     = 0.01
bres@tmYLMode             = "Explicit"
bres@tmYLValues           = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
bres@tmYLLabels           = y
bres@tiYAxisString         = "MJO phase"
bres@tiXAxisString         = "lag days" 
bres@tiXAxisOffsetYF       = 0.03
bres@tiYAxisFontHeightF    = 0.015

bres@trXMaxF  = n
bres@trXMinF  = 0

bres@tmXBMinorPerMajor     = 8
bres@tmXBLabelFontHeightF  = 0.015
bres@tmXBMajorLengthF      = 0.01
bres@tmXBTickSpacingF      = 10.
bres@tmXBMode              = "Manual"
bres@tiXAxisFontHeightF    = 0.015

wks = gsn_open_wks("x11",plotdir)
plot = new(3, graphic)

do k = 0,2,1  
    bres@gsnLeftString = Title(k)
    plot(k) = gsn_csm_blank_plot(wks,bres)
    plot(k) = add_polygon(wks,plot(k),plag(k,:,:),label(k,:),cmap)
    add_labelbar(wks,plot(k),label(k,:),cmap)
end do

; Panel these 3 plots.
resP                  = True
resP@gsnMaximize      = False
resP@gsnFrame         = False
resP@gsnDraw          = True
resP@gsnPanelLabelBar = False
resP@gsnPanelTop      = 0.98
resP@gsnPanelBottom   = 0.02
resP@gsnPanelRowSpec  = True
gsn_panel(wks,plot,(/1,1,1/),resP)

i = 1
do while (i .eq.1)
    i = 1
end do