load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
;-----------------------------------------------------------------------------
;
;-----------------------------------------------------------------------------
ddir = "/home/sunming/data5/cuixy/ENSO-MJO-Jet/data"
udir = ddir+"/u_DJF_mean_1979-2022.nc"
;; jet indices.
eadir = ddir+"/jet_indices_EAJS.nc"
medir = ddir+"/jet_indices_MEJS.nc"
nadir = ddir+"/jet_indices_NAJS.nc"

sigeadir = ddir+"/monte_carlo/jet_indices_EAJS.nc"
sigmedir = ddir+"/monte_carlo/jet_indices_MEJS.nc"
signadir = ddir+"/monte_carlo/jet_indices_NAJS.nc"

rmmdir  = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
ninodir = "/home/sunming/data5/cuixy/data/clim_indices/Nino_3.4_index.csv"

pdir = "/home/sunming/data5/cuixy/ENSO-MJO-Jet/figure/uclim_jet_indices"

nwaa = 90.
nwbb = 10.

wks = gsn_open_wks("x11",pdir)
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
  xClimsm = smthClmDayTLL(xClim, 3)
  xAnom = calcDayAnomTLL(xtmp, yyyyddd, xClimsm)
  xr    = xAnom(:,0,0)
  xr := wgt_runave_Wrap(xr, wts, 0)

  return(xr)

end
;-----------------------------------------------------------------------------
; Function to draw polygons in a viewpoint space.
;-----------------------------------------------------------------------------
undef("add_polygon")
function add_polygon(wks:graphic,plot:graphic,values,ranges,color)
local plot,values,xmin,xmax,ymin,ymax,nx,ny,fil_color
begin
  dims = dimsizes(values)
  ny   = dims(1)
  nx   = dims(0)
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
procedure add_sigpoint(wks:graphic,plot:graphic,plagt,mvalue,sigcolor)
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
  gnres@gsMarkerColor = sigcolor
; First row.
  j = 0
  ypos = ymin
  do i=1,nx-2
    if(plagt(i,j).gt.mvalue(1,i,j).or.plagt(i,j).lt.mvalue(0,i,j))then
      xpos = xmin + (xsize/2. * (2*i-1)) +0.5   
      str = unique_string("fill_dot")     
      plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
    end if   
  end do
; Left
  i = 0
  xpos = xmin +0.5
  if(plagt(i,j).gt.mvalue(1,i,j).or.plagt(i,j).lt.mvalue(0,i,j))then
    str = unique_string("fill_dot")     
    plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
  end if
; The rightmost.
  i = nx-1
  xpos = xmin + (xsize/2. * (2*i-1))+0.5
  if(plagt(i,j).gt.mvalue(1,i,j).or.plagt(i,j).lt.mvalue(0,i,j))then
    str = unique_string("fill_dot")     
    plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
  end if

; Top
  j = ny - 1
  ypos = ymax - ysize/2.+0.5
  do i=1,nx-2
    xpos = xmin + (xsize/2. * (2*i-1))+0.5
    if(plagt(i,j).gt.mvalue(1,i,j).or.plagt(i,j).lt.mvalue(0,i,j))then   
      str = unique_string("fill_dot")     
      plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
    end if      
  end do
; Left
  i = 0
  xpos = xmin +0.5
  if(plagt(i,j).gt.mvalue(1,i,j).or.plagt(i,j).lt.mvalue(0,i,j))then
    str = unique_string("fill_dot")     
    plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
  end if
; The rightmost.
  i = nx-1
  xpos = xmin + (xsize/2. * (2*i-1))+0.5
  if(plagt(i,j).gt.mvalue(1,i,j).or.plagt(i,j).lt.mvalue(0,i,j))then
    str = unique_string("fill_dot")     
    plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
  end if
;Left column
  i = 0
  xpos = xmin+0.5
  do j = 1,ny-2
    ypos = ymin + (ysize/2.*(2*j-1))+0.5
    if(plagt(i,j).gt.mvalue(1,i,j).or.plagt(i,j).lt.mvalue(0,i,j))then
      str = unique_string("fill_dot")     
      plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
    end if  
  end do
; Right column
  i = nx-1
  xpos = xmax - xsize/2.+0.5
  do j = 1,ny-2
    ypos = ymin + (ysize/2.*(2*j-1))+0.5
    if(plagt(i,j).gt.mvalue(1,i,j).or.plagt(i,j).lt.mvalue(0,i,j))then
      str = unique_string("fill_dot")     
      plot@$str$ = gsn_add_polymarker(wks,plot,xpos,ypos,gnres)
    end if  
  end do
;Middle values.
  do j=1,ny-2
    ypos = ymin + (ysize/2.*(2*j-1))+0.5
    do i=1,nx-2
      xpos = xmin + (xsize/2.*(2*i-1))+0.5
      if(plagt(i,j).gt.mvalue(1,i,j).or.plagt(i,j).lt.mvalue(0,i,j))then
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

  lbres@vpWidthF   = 0.45*vpw
  lbres@vpHeightF  = 1.0*vph

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
  lbres@lbLabelFont        = 21
  lbres@lbBoxLinesOn       = True
  lbres@lbFillOpacityF     = 1.0
  lbres@lbLabelStride      = 1.0
  lbres@lbLabelFontHeightF = 0.01
  lbres@lbLabelAngleF      = 0.
  lbres@lbBoxSeparatorLinesOn = True
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
;----------------------------------------------------------------
; jet indices plot.
  fea = addfile(eadir,"r")
  fme = addfile(medir,"r")
  fna = addfile(nadir,"r")

  time = fea->time
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
  delete([/phase,amp/])
;--------------------------------------------------------------
  var1 = "jets"
  var2 = "jetl"
;; NAJS.
  jetsna1 = fna->$var1$ ; jet speed
  jetlna1 = fna->$var2$ ; jet latitude

  jetlna = ISO_Filter(jetlna1,time,201,nwaa,nwbb)
  jetsna = ISO_Filter(jetsna1,time,201,nwaa,nwbb)

;; EAJS.
  jetsea1 = fea->$var1$ ; jet speed
  jetlea1 = fea->$var2$ ; jet latitude

  jetlea = ISO_Filter(jetlea1,time,201,nwaa,nwbb)
  jetsea = ISO_Filter(jetsea1,time,201,nwaa,nwbb)
;; MEJS.
  jetsme1 = fme->$var1$ ; jet speed
  jetlme1 = fme->$var2$ ; jet latitude

  jetlme = ISO_Filter(jetlme1,time,201,nwaa,nwbb)
  jetsme = ISO_Filter(jetsme1,time,201,nwaa,nwbb)

;first conver from 12 to 3.
;; time: 19791201--20220331
  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.31)
  reamp   := reamp(n1:n2)
  rephase := rephase(n1:n2)
  jetlna  := jetlna(n1:n2)
  jetsna  := jetsna(n1:n2)

  jetlea  := jetlea(n1:n2)
  jetsea  := jetsea(n1:n2)

  jetlme  := jetlme(n1:n2)
  jetsme  := jetsme(n1:n2)

  year    := year(n1:n2)
  mon     := mon(n1:n2)
  time    := time(n1:n2)

;; select for EI, LA or NE. Dec year + 1.

  delete(iyear)
  delim  = ","
  ninof  = asciiread(ninodir,-1,"string")
  iyear  = str_get_field(ninof,1,delim)
  ievent = str_get_field(ninof,6,delim)

  iyear  := toint(iyear(1::))
  ievent := tofloat(ievent(1::))

;;EI
  nei = ind(ievent.eq.1)
  yei = iyear(nei)
;;LA
  nla = ind(ievent.eq.-1)
  yla = iyear(nla)
;;NE
  nne = ind(ievent.eq.0)
  yne = iyear(nne)

  nt = ind(mon.eq.12)
  year(nt) = year(nt)+1
;-----------------------------------------------------------------------------
;calculate phase-lag data. lag time 50 days
  n = 30
  plag = new((/3,2,n+1,13/),float) ; (indices,lag days, phases)
  plag = 0.

  jetl = jetlna
  jets = jetsna

;; significant level.
  sigl = 0.05
  sigh = 0.95
  fm = addfile(signadir,"r")
  mvalue1 = fm->mvalue  ;(times,lei,:,:)
  mvalue  = new((/3,2,2,31,13/),float)  ;(0.25/97.5,lei,:,:)
  mvalue  = 0.
  n05  = toint(dimsizes(mvalue1(0,:,0,0,0))*sigl) - 1
  n95  = toint(dimsizes(mvalue1(0,:,0,0,0))*sigh) - 1
  mvalue(:,0,:,:,:) = mvalue1(:,n05,:,:,:)
  mvalue(:,1,:,:,:) = mvalue1(:,n95,:,:,:)

  sigl = 0.1
  sigh = 0.9
  mvaluet  = new((/3,2,2,31,13/),float)  ;(0.25/97.5,lei,:,:)
  mvaluet  = 0.
  n10  = toint(dimsizes(mvalue1(0,:,0,0,0))*sigl) - 1
  n90  = toint(dimsizes(mvalue1(0,:,0,0,0))*sigh) - 1
  mvaluet(:,0,:,:,:) = mvalue1(:,n10,:,:,:)
  mvaluet(:,1,:,:,:) = mvalue1(:,n90,:,:,:)


;; here for EI.
  ny = yei
  nc = ind(year.eq.ny(0))
  do i = 1, dimsizes(ny)-1
    ntemp = ind(year.eq.ny(i))
    nc := array_append_record(nc, ntemp, 0)
    delete(ntemp)
  end do

  reampt = reamp(nc)
  rephat = rephase(nc)
  jetlt  = jetl(nc)
  jetst  = jets(nc)

  yeart = year(nc)
  mont  = mon(nc)
  timet = time(nc) 

  nd   = ind(mont.eq.12.or.mont.eq.1.or.mont.eq.2)
  amp   = reampt(nd)
  phase = rephat(nd)
  in1  = ispan(0,dimsizes(timet)-1,1)
  in2  = in1(nd)

  delete(nt)
  am = 1.0
  do i = 1,8
    nt  = ind(phase.eq.i.and.amp.gt.am)
    nt1 = in2(nt)
    do j = 0 ,n, 1
      plag(0,0,j,i-1) = dim_avg(jetst(nt1+j))
      plag(0,1,j,i-1) = dim_avg(jetlt(nt1+j))
    end do
    delete([/nt,nt1/])
  end do

  do i = 0,4
    plag(0,:,:,i+8) = plag(0,:,:,i)
  end do

  delete([/ny,nc,reampt,rephat,jetlt,jetst,yeart,mont,timet,nd,\
         amp,phase,in1,in2/])

;; for ne.
  ny = yne
  nc = ind(year.eq.ny(0))
  do i = 1, dimsizes(ny)-1
    ntemp = ind(year.eq.ny(i))
    nc := array_append_record(nc, ntemp, 0)
    delete(ntemp)
  end do

  reampt = reamp(nc)
  rephat = rephase(nc)
  jetlt  = jetl(nc)
  jetst  = jets(nc)

  yeart = year(nc)
  mont  = mon(nc)
  timet = time(nc) 

  nd   = ind(mont.eq.12.or.mont.eq.1.or.mont.eq.2)
  amp   = reampt(nd)
  phase = rephat(nd)
  in1  = ispan(0,dimsizes(timet)-1,1)
  in2  = in1(nd)

  delete(nt)
  am = 1.0
  do i = 1,8
    nt  = ind(phase.eq.i.and.amp.gt.am)
    nt1 = in2(nt)
    do j = 0 ,n, 1
      plag(1,0,j,i-1) = dim_avg(jetst(nt1+j))
      plag(1,1,j,i-1) = dim_avg(jetlt(nt1+j))
    end do
    delete([/nt,nt1/])
  end do

  do i = 0,4
    plag(1,:,:,i+8) = plag(1,:,:,i)
  end do

  delete([/ny,nc,reampt,rephat,jetlt,jetst,yeart,mont,timet,nd,\
         amp,phase,in1,in2/])

;; for LA.
;; here for EI.
  ny = yla
  nc = ind(year.eq.ny(0))
  do i = 1, dimsizes(ny)-1
    ntemp = ind(year.eq.ny(i))
    nc := array_append_record(nc, ntemp, 0)
    delete(ntemp)
  end do

  reampt = reamp(nc)
  rephat = rephase(nc)
  jetlt  = jetl(nc)
  jetst  = jets(nc)

  yeart = year(nc)
  mont  = mon(nc)
  timet = time(nc) 

  nd   = ind(mont.eq.12.or.mont.eq.1.or.mont.eq.2)
  amp   = reampt(nd)
  phase = rephat(nd)
  in1  = ispan(0,dimsizes(timet)-1,1)
  in2  = in1(nd)

  delete(nt)
  am = 1.0
  do i = 1,8
    nt  = ind(phase.eq.i.and.amp.gt.am)
    nt1 = in2(nt)
    do j = 0 ,n, 1
      plag(2,0,j,i-1) = dim_avg(jetst(nt1+j))
      plag(2,1,j,i-1) = dim_avg(jetlt(nt1+j))
    end do
    delete([/nt,nt1/])
  end do

  do i = 0,4
    plag(2,:,:,i+8) = plag(2,:,:,i)
  end do

  delete([/ny,nc,reampt,rephat,jetlt,jetst,yeart,mont,timet,nd,\
         amp,phase,in1,in2/])

  x = ispan(0,n,1)  
  y = (/1,2,3,4,5,6,7,8,1,2,3,4,5/)
  color1 = "MPL_BrBG"
  color2 = "MPL_PuOr"
  undef("cmap")
  cmap = read_colormap_file("MPL_BrBG")

;; MEJS js:3 jl:1 EAJS js:4 jl:1 NAJS js:2 jl:1.4
  jsbe =  2
  jsen = -2
  jlbe =  1.4
  jlen = -1.4
;label1  EI jets.
  ncolor = toint((jsbe + jsbe)/0.2 + 2)
  cmap1  = span_color_rgba(color1, ncolor)
  label1 = fspan(jsen, jsbe, ncolor-1)

;label2 EI jetl.
  ncolor = toint((jlbe + jlbe)/0.1 + 2)
  cmap2  = span_color_rgba(color2, ncolor)
  label2 = fspan(jlen, jlbe, ncolor-1)

;label3 NE jets
  ncolor = toint((jsbe + jsbe)/0.2 + 2)
  cmap3  = span_color_rgba(color1, ncolor)
  label3 = fspan(jsen, jsbe, ncolor-1)

;label4 NE jetl
  ncolor = toint((jlbe + jlbe)/0.1 + 2)
  cmap4  = span_color_rgba(color2, ncolor)
  label4 = fspan(jlen, jlbe, ncolor-1)

;label5 LA jets
  ncolor = toint((jsbe + jsbe)/0.2 + 2)
  cmap5  = span_color_rgba(color1, ncolor)
  label5 = fspan(jsen, jsbe, ncolor-1)

;label6 LA jetl
  ncolor = toint((jlbe + jlbe)/0.1 + 2)
  cmap6  = span_color_rgba(color2, ncolor)
  label6 = fspan(jlen, jlbe, ncolor-1)

  ;i = 0
  ;do while(i.eq.0)
  ;  i = 0
  ;end do

;---Resources for blank plot  
  bres                = True                    
  bres@gsnDraw        = False
  bres@gsnFrame       = False
  bres@gsnMaximize    = False     
  bres@gsnLeftStringFontHeightF  = 0.012
  bres@gsnRightStringFontHeightF = 0.012

  bres@vpWidthF       = 0.35
  bres@vpHeightF      = 0.15

  bres@trYMaxF  = 12
  bres@trYMinF  = 0
  bres@tmYLLabelFontHeightF = 0.010
  bres@tmYLMinorPerMajor    = 0
  bres@tmYLTickSpacingF     = 1
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
  bres@tmXBMajorLengthF      = 0.005
  bres@tmXBTickSpacingF      = 5.
  bres@tmXBMode              = "Manual"
  bres@tiXAxisFontHeightF    = 0.012

  bres@tmXBLabelFont = 21
  bres@tmYRLabelFont = 21

  bres@vpWidthF  = 0.23
  bres@vpHeightF = 0.18

; EI
  plot1 = new(2, graphic)

  k = 0
  bres@vpXF  = 0.05
  bres@vpYF  = 0.7

  bres@tiYAxisString         = "~F21~MJO phase"
  bres@tiXAxisString         = "" 
  bres@gsnLeftString = "~F21~(a) ~F21~JSI"
  bres@gsnRightString = "~F21~EI"
  plot1(k) = gsn_csm_blank_plot(wks,bres)
  plot1(k) = add_polygon(wks,plot1(k),plag(0,k,:,:),label1,cmap)
  ;add_sigpoint(wks,plot1(k),plag(0,k,:,:),mvaluet(0,:,k,:,:),"black")
  add_sigpoint(wks,plot1(k),plag(0,k,:,:),mvalue(0,:,k,:,:),"red")
  add_labelbar(wks,plot1(k),label1,cmap1)
  print("plot d finished")

  k = 1
  bres@vpXF  = 0.05
  bres@vpYF  = 0.45

  bres@tiYAxisString         = ""
  bres@tiXAxisString         = "~F21~Lag days" 
  bres@gsnLeftString = "~F21~(b) ~F21~JLI"
  plot1(k) = gsn_csm_blank_plot(wks,bres)
  plot1(k) = add_polygon(wks,plot1(k),plag(0,k,:,:),label2,color2)
  ;add_sigpoint(wks,plot1(k),plag(0,k,:,:),mvaluet(0,:,k,:,:),"black")
  add_sigpoint(wks,plot1(k),plag(0,k,:,:),mvalue(0,:,k,:,:),"red")
  add_labelbar(wks,plot1(k),label2,cmap2)
  print("plot g finished")

;; NE
  plot2 = new(2, graphic)

  k = 0
  bres@vpXF  = 0.38
  bres@vpYF  = 0.7

  bres@tiYAxisString         = ""
  bres@tiXAxisString         = "" 
  bres@gsnLeftString = "~F21~(c) ~F21~JSI"
  bres@gsnRightString = "~F21~NE"
  plot2(k) = gsn_csm_blank_plot(wks,bres)
  plot2(k) = add_polygon(wks,plot2(k),plag(1,k,:,:),label3,cmap)
  ;add_sigpoint(wks,plot2(k),plag(1,k,:,:),mvaluet(1,:,k,:,:),"black")
  add_sigpoint(wks,plot2(k),plag(1,k,:,:),mvalue(1,:,k,:,:),"red")
  add_labelbar(wks,plot2(k),label3,cmap3)
  print("plot c finished")

  k = 1
  bres@vpXF  = 0.38
  bres@vpYF  = 0.45

  bres@tiYAxisString = ""
  bres@tiXAxisString = "~F21~Lag days" 
  bres@gsnLeftString = "~F21~(d) ~F21~JLI"
  plot2(k) = gsn_csm_blank_plot(wks,bres)
  plot2(k) = add_polygon(wks,plot2(k),plag(1,k,:,:),label4,color2)
  ;add_sigpoint(wks,plot2(k),plag(1,k,:,:),mvaluet(1,:,k,:,:),"black")
  add_sigpoint(wks,plot2(k),plag(1,k,:,:),mvalue(1,:,k,:,:),"red")
  add_labelbar(wks,plot2(k),label4,cmap4)
  print("plot f finished")

;; LA
  plot3 = new(3, graphic)

  k = 0
  bres@vpXF  = 0.7
  bres@vpYF  = 0.7

  bres@tiYAxisString = ""
  bres@tiXAxisString = "" 
  bres@gsnLeftString = "~F21~(e) ~F21~JSI"
  bres@gsnRightString = "~F21~LA"
  plot3(k) = gsn_csm_blank_plot(wks,bres)
  plot3(k) = add_polygon(wks,plot3(k),plag(2,k,:,:),label5,cmap)
  ;add_sigpoint(wks,plot3(k),plag(2,k,:,:),mvaluet(2,:,k,:,:),"black")
  add_sigpoint(wks,plot3(k),plag(2,k,:,:),mvalue(2,:,k,:,:),"red")
  add_labelbar(wks,plot3(k),label5,cmap5)
  print("plot b finished")

  k = 1
  bres@vpXF  = 0.7
  bres@vpYF  = 0.45

  bres@tiYAxisString = ""
  bres@tiXAxisString = "~F21~Lag days" 
  bres@gsnLeftString = "~F21~(f) ~F21~JLI"
  plot3(k) = gsn_csm_blank_plot(wks,bres)
  plot3(k) = add_polygon(wks,plot3(k),plag(2,k,:,:),label6,color2)
  ;add_sigpoint(wks,plot3(k),plag(2,k,:,:),mvaluet(2,:,k,:,:),"black")
  add_sigpoint(wks,plot3(k),plag(2,k,:,:),mvalue(2,:,k,:,:),"red")
  add_labelbar(wks,plot3(k),label6,cmap6)
  print("plot e finished")

  draw(plot1(0))
  draw(plot1(1))
  draw(plot2(0))
  draw(plot2(1))
  draw(plot3(0))
  draw(plot3(1))

  i = 0
  do while(i.eq.0)
    i = 0
  end do