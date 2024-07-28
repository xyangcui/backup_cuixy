
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
datadir = getenv("indir")
mdir    = getenv("outdir")
var1    = getenv("var1")
var2    = getenv("var2")
var3    = getenv("var3")

rmmdir  = getenv("rmmdir")
plotdir = getenv("plotdir")

nwaa = tofloat(getenv("n1"))
nwbb = tofloat(getenv("n2"))
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
; procedure to attach dot.
; mkres@gsMarkerIndex = 16
;-----------------------------------------------------------------------------
undef("add_sigpoint")
procedure add_sigpoint(wks:graphic,plot:graphic,plagt,mvalue)
local plot,plagt,xmin,xmax,ymin,ymax,nx,ny,mvalue
begin
    dims    = dimsizes(plagt)
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
    gnres@gsMarkerIndex = 16 
    gnres@gsMarkerSizeF = 0.003
; First row.
    j = 0
    ypos = ymin
    do i=1,nx-2
        if(plagt(i,j).ge.mvalue(1,i,j).or.plagt(i,j).le.mvalue(0,i,j))then
          xpos = xmin + (xsize/2. * (2*i-1)) +0.5   
          str = unique_string("fill_dot")     
          plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
        end if   
    end do
; Left
    i = 0
    xpos = xmin +0.5
    if(plagt(i,j).ge.mvalue(1,i,j).or.plagt(i,j).le.mvalue(0,i,j))then
      str = unique_string("fill_dot")     
      plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
    end if
; The rightmost.
    i = nx-1
    xpos = xmin + (xsize/2. * (2*i-1))+0.5
    if(plagt(i,j).ge.mvalue(1,i,j).or.plagt(i,j).le.mvalue(0,i,j))then
      str = unique_string("fill_dot")     
      plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
    end if

; Top
    j = ny - 1
    ypos = ymax - ysize/2.+0.5
    do i=1,nx-2
      xpos = xmin + (xsize/2. * (2*i-1))+0.5
      if(plagt(i,j).ge.mvalue(1,i,j).or.plagt(i,j).le.mvalue(0,i,j))then   
        str = unique_string("fill_dot")     
        plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
      end if      
    end do
; Left
    i = 0
    xpos = xmin +0.5
    if(plagt(i,j).ge.mvalue(1,i,j).or.plagt(i,j).le.mvalue(0,i,j))then
      str = unique_string("fill_dot")     
      plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
    end if
; The rightmost.
    i = nx-1
    xpos = xmin + (xsize/2. * (2*i-1))+0.5
    if(plagt(i,j).ge.mvalue(1,i,j).or.plagt(i,j).le.mvalue(0,i,j))then
      str = unique_string("fill_dot")     
      plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
    end if
;Left column
    i = 0
    xpos = xmin+0.5
    do j = 1,ny-2
        ypos = ymin + (ysize/2.*(2*j-1))+0.5
        if(plagt(i,j).ge.mvalue(1,i,j).or.plagt(i,j).le.mvalue(0,i,j))then
            str = unique_string("fill_dot")     
            plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
        end if  
    end do
; Right column
    i = nx-1
    xpos = xmax - xsize/2.+0.5
    do j = 1,ny-2
        ypos = ymin + (ysize/2.*(2*j-1))+0.5
        if(plagt(i,j).ge.mvalue(1,i,j).or.plagt(i,j).le.mvalue(0,i,j))then
            str = unique_string("fill_dot")     
            plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
        end if  
    end do
;Middle values.
    do j=1,ny-2
        ypos = ymin + (ysize/2.*(2*j-1))+0.5
        do i=1,nx-2
            xpos = xmin + (xsize/2.*(2*i-1))+0.5
            if(plagt(i,j).ge.mvalue(1,i,j).or.plagt(i,j).le.mvalue(0,i,j))then
                str = unique_string("fill_dot")     
                plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
            end if  
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
    lbres@lbLabelStride      = 4
    lbres@lbLabelFontHeightF = 0.01
    lbres@lbLabelAngleF      = -45
    lbres@lbBoxSeparatorLinesOn = False
    lbid = gsn_create_labelbar(wks,dimsizes(color(:,0)),labels,lbres)
; attach labelbar to the plot.
    amres = True
    amres@amJust = "CenterLeft"
    amres@amParallelPosF   = 0.42
    amres@amOrthogonalPosF = 0.0
    annoid = gsn_add_annotation(plot,lbid,amres)
end 

;-----------------------------------------------------------------------------
; Main code.
;-----------------------------------------------------------------------------
;read in data and selected strong MJO days.
;; m_value.
fm = addfile(mdir,"r")
m_value1 = fm->m_value  ;(times,lei,:,:)
m_value  = new((/2,3,31,13/),float)  ;(0.25/97.5,lei,:,:)
m_value  = 0.
n25  = toint(dimsizes(m_value1(:,0,0,0))*0.05) - 1
n975 = toint(dimsizes(m_value1(:,0,0,0))*0.95) - 1
m_value(0,:,:,:) = m_value1(n25,:,:,:)
m_value(1,:,:,:) = m_value1(n975,:,:,:)

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
jetz1 = f->$var3$ ; jet zonality

jetl = ISO_Filter(jetl1,time,201,nwaa,nwbb)
jets = ISO_Filter(jets1,time,201,nwaa,nwbb)
jetz = ISO_Filter(jetz1,time,201,nwaa,nwbb)
;-----------------------------------------------------------------------------
;calculate phase-lag data. lag time 50 days
n = 30
plag   = new((/3,n+1,13/),float) ; (indices,lag days, phases)
plag   = 0.
plag1  = plag    ; amp = 1.5
plag2  = plag    ; amp = 1.0

;first conver from 12 to 3.
n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.31)
reamp   := reamp(n1:n2)
rephase := rephase(n1:n2)
jetl    := jetl(n1:n2)
jets    := jets(n1:n2)
jetz    := jetz(n1:n2)
mon     := mon(n1:n2)
time    := time(n1:n2)

nd   = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
reamp    := reamp(nd)
ndd  = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
rephase  := rephase(ndd)
in1  = ispan(0,dimsizes(time)-1,1)
in2  = in1(nd)

;amp > 2.0
am = 1.0
do i = 1,8
    nt  = ind(rephase.eq.i.and.reamp.gt.am)
    nt1 = in2(nt)
    do j = 0 ,n, 1
        plag(0,j,i-1) = dim_avg(jets(nt1+j))
        plag(1,j,i-1) = dim_avg(jetl(nt1+j))
        plag(2,j,i-1) = dim_avg(jetz(nt1+j))
    end do
    delete([/nt,nt1/])
end do

do i = 0,4
    plag(:,:,i+8) = plag(:,:,i)
end do


x = ispan(0,n,1)  ; x labels
y = (/1,2,3,4,5,6,7,8,1,2,3,4,5/)
color1 = "MPL_PRGn"
color2 = "MPL_PuOr"
color3 = "MPL_PiYG"
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
be = floor(max(plag(2,:,:)))
en = ceil(min(plag(2,:,:)))
ll = (/be,abs(en)/)
a  = max(ll)
if (be.eq.abs(en).and.a.eq.0) then
    ncolor = toint((1 + 1)/0.1 + 2)
    cmap3  = span_color_rgba(color3, ncolor)
    label3 = fspan(-1., 1., ncolor-1)
else
    ncolor = toint((a + a)/0.2 + 2)
    cmap3  = span_color_rgba(color3, ncolor)
    label3 = fspan(-a, a, ncolor-1)
end if

;---Resources for blank plot  
bres                = True                    
bres@gsnDraw        = False
bres@gsnFrame       = False
bres@gsnMaximize    = False     
bres@gsnLeftStringFontHeightF  = 0.012
bres@gsnRightStringFontHeightF = 0.012

bres@vpWidthF       = 0.23
bres@vpHeightF      = 0.18

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


wks = gsn_open_wks("x11", plotdir)
; amp > 2.0
bres@gsnRightString = "~F21~Amp>2"
plot1 = new(3, graphic)
k = 0
bres@vpXF  = 0.07
bres@vpYF  = 0.95
bres@tiYAxisString         = "~F21~MJO phase"
bres@tiXAxisString         = "" 
bres@gsnLeftString = "~F22~a ~F21~Jet speed"
plot1(k) = gsn_csm_blank_plot(wks,bres)
plot1(k) = add_polygon(wks,plot1(k),plag(k,:,:),label1,color1)
add_sigpoint(wks,plot1(k),plag(k,:,:),m_value(:,k,:,:))
add_labelbar(wks,plot1(k),label1,cmap1)
print("plot a finished")

k = 1
bres@vpXF  = 0.385
bres@vpYF  = 0.95
bres@tiYAxisString         = ""
bres@tiXAxisString         = "" 
bres@gsnLeftString = "~F22~d ~F21~Jet latitude"
plot1(k) = gsn_csm_blank_plot(wks,bres)
plot1(k) = add_polygon(wks,plot1(k),plag(k,:,:),label2,color2)
add_sigpoint(wks,plot1(k),plag(k,:,:),m_value(:,k,:,:))
add_labelbar(wks,plot1(k),label2,cmap2)
print("plot d finished")

k = 2
bres@vpXF  = 0.7
bres@vpYF  = 0.95
bres@tiYAxisString         = ""
bres@tiXAxisString         = "" 
bres@gsnLeftString = "~F22~g ~F21~Jet zonality"
plot1(k) = gsn_csm_blank_plot(wks,bres)
plot1(k) = add_polygon(wks,plot1(k),plag(k,:,:),label3,color3)
add_sigpoint(wks,plot1(k),plag(k,:,:),m_value(:,k,:,:))
add_labelbar(wks,plot1(k),label3,cmap3)
print("plot g finished")

do k = 0, 2
    draw(plot1(k))
end do

i = 0
do while(i.eq.0)
    i = 0
end do
