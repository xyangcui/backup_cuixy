;-----------------------------------------------------------------
;; draw general budget terms. 
;; MEJS box (10-70E,10-40N)
;; EAJS box (110-180E,15-50N)
;; time: 2023-12-13
;-----------------------------------------------------------------
  rad   = 4.*atan(1.)/180. 
  tday  = 24.*60.*60.
  mlatn = 30
  mlats = 25
  mlonw = 70
  mlone = 100

  ddir = "/home/sunming/data5/cuixy/global_jets/data/temp_budget"
  tdir = ddir+"/"+"anom_dtdt_r2.5.nc"
  udir = ddir+"/"+"anom_udtdx_r2.5.nc"
  vdir = ddir+"/"+"anom_vdtdy_r2.5.nc"
  wdir = ddir+"/"+"anom_wdtdp_r2.5.nc"

  t1dir = "/home/sunming/data5/cuixy/global_jets/data/ERA5_daily_anom_tvert_1979-2022_2.5x2.5.nc"
  pdir  = "/home/sunming/data5/cuixy/global_jets/plot/Fig4_bar_MEJS"
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
  ft  = addfile(tdir,"r")
  fu  = addfile(udir,"r")
  fv  = addfile(vdir,"r")
  fw  = addfile(wdir,"r")

  ft1 = addfile(t1dir,"r")
  t1  = ft1->t

  t = ft->dtdt
  t = t
  u = fu->udtdx
  u = u*tday
  v = fv->vdtdy
  v = v*tday
  w = fw->wdtdp
  w = w*tday

;residual
  r = t + (u + v + w) ;- r1)
  copy_VarCoords(t, r)
  u = -1*u
  v = -1*v
  w = -1*w;+r1
  lat = t&lat

  adir = ddir+"/anom_udtdx_seperate_r2.5.nc"

  fa  = addfile(adir,"r")

  ultl = fa->ultl
  ulti = fa->ulti
  ulth = fa->ulth
  uitl = fa->uitl
  uiti = fa->uiti
  uith = fa->uith
  uhtl = fa->uhtl
  uhti = fa->uhti
  uhth = fa->uhth

;---------------------------------------------------------------=
;; MJO phases.
  time = ft->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
  delim  = ","
  rmmf   = asciiread(rmmdir,-1,"string")
  iyear  = str_get_field(rmmf,1,delim)
  imon   = str_get_field(rmmf,2,delim)
  iday   = str_get_field(rmmf,3,delim)
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

;; store lead-lag values for plotting. 
;;MEJS.
  llmejs = new((/8,41,6/),float)
  llmejs = 0.
  tm = t(:,{mlats:mlatn},{mlonw:mlone})
  clatjs = cos(tm&lat*rad)

  t1js = wgt_areaave_Wrap(t1(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  tjs = wgt_areaave_Wrap(t(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  ujs = wgt_areaave_Wrap(u(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  vjs = wgt_areaave_Wrap(v(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  wjs = wgt_areaave_Wrap(w(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  rjs = wgt_areaave_Wrap(r(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)

  do i = 0, 7
    ip = i + 1
    nt = ind(phase.eq.ip.and.amp.gt.1.)
    nind = ind_DJF(nt)
    do j = 0, 40
      ntt   = -15 + j
      nindd = nind + ntt
      llmejs(i,j,0) = dim_avg_n(tjs(nindd), 0)
      llmejs(i,j,1) = dim_avg_n(ujs(nindd), 0)
      llmejs(i,j,2) = dim_avg_n(vjs(nindd), 0)
      llmejs(i,j,3) = dim_avg_n(wjs(nindd), 0)
      llmejs(i,j,4) = dim_avg_n(rjs(nindd), 0)
      llmejs(i,j,5) = dim_avg_n(t1js(nindd), 0)
      delete([/ntt,nindd/])
    end do
    delete([/ip,nt,nind/])
  end do

  delete([/clatjs,tjs,ujs,vjs,wjs,rjs/])
  print("calc MEJS data done. "+ systemfunc("date"))

;; line for dtdt
  line = new((/2,41/),float)
  line(0,:) = llmejs(2,:,0)
  line(1,:) = llmejs(2,:,5)
;; bar for other value. EAJS:27-32 day: 12-17 MEJS: day 2-7 17-22 ;;
  bar = new((/5/),float)
  bar = dim_avg_n(llmejs(2,6:15,0:4),0)

;; seperate
;; store lead-lag values for plotting. 
;;MEJS.
  lsmejs = new((/10/),float)
  lsmejs = 0.

  tm = t(:,{mlats:mlatn},{mlonw:mlone})
  clatjs = cos(tm&lat*rad)
  
  ujs = wgt_areaave_Wrap(u(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  ultljs = wgt_areaave_Wrap(ultl(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  ultijs = wgt_areaave_Wrap(ulti(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  ulthjs = wgt_areaave_Wrap(ulth(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  uitljs = wgt_areaave_Wrap(uitl(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  uitijs = wgt_areaave_Wrap(uiti(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  uithjs = wgt_areaave_Wrap(uith(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  uhtljs = wgt_areaave_Wrap(uhtl(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  uhtijs = wgt_areaave_Wrap(uhti(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)
  uhthjs = wgt_areaave_Wrap(uhth(:,{mlats:mlatn},{mlonw:mlone})\
                         , clatjs, 1.0, 1)

  nt = ind(phase.eq.3.and.amp.gt.1.)
  nind  = ind_DJF(nt)
  ntt   = -9
  nin = nind + ntt
  do j = 1, 9
    ntt   = -9 + j
    nindd = nind + ntt
    nin  := array_append_record(nin, nindd, 0)
  end do
  lsmejs(0) = dim_avg_n(ujs(nin), 0)
  lsmejs(1) = dim_avg_n(ultljs(nin), 0)*tday*(-1)
  lsmejs(2) = dim_avg_n(ultijs(nin), 0)*tday*(-1)
  lsmejs(3) = dim_avg_n(ulthjs(nin), 0)*tday*(-1)
  lsmejs(4) = dim_avg_n(uitljs(nin), 0)*tday*(-1)
  lsmejs(5) = dim_avg_n(uitijs(nin), 0)*tday*(-1)
  lsmejs(6) = dim_avg_n(uithjs(nin), 0)*tday*(-1)
  lsmejs(7) = dim_avg_n(uhtljs(nin), 0)*tday*(-1)
  lsmejs(8) = dim_avg_n(uhtijs(nin), 0)*tday*(-1)
  lsmejs(9) = dim_avg_n(uhthjs(nin), 0)*tday*(-1)

  delete([/clatjs,ujs/])
  print("calc MEJS data done. "+ systemfunc("date"))



;-----------------------------------------------------------------------------
; Function to draw polygons in a viewpoint space.
;-----------------------------------------------------------------------------
undef("add_polygon")
function add_polygon(wks:graphic,plot:graphic,values,ranges,color)
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
;; m_value.
;; EASJ.
sigl = 0.05
sigh = 0.95
fmea = addfile(meadir,"r")
mea_value1 = fmea->m_value  ;(times,lei,:,:)
mea_value  = new((/2,2,31,13/),float)  ;(0.25/97.5,lei,:,:)
mea_value  = 0.
n10  = toint(dimsizes(mea_value1(:,0,0,0))*sigl) - 1
n90  = toint(dimsizes(mea_value1(:,0,0,0))*sigh) - 1
mea_value(0,:,:,:) = mea_value1(n10,:,:,:)
mea_value(1,:,:,:) = mea_value1(n90,:,:,:)
;; MESJ.
fmme = addfile(mmedir,"r")
mme_value1 = fmme->m_value  ;(times,lei,:,:)
mme_value  = new((/2,2,31,13/),float)  ;(0.25/97.5,lei,:,:)
mme_value  = 0.
n10  = toint(dimsizes(mme_value1(:,0,0,0))*sigl) - 1
n90  = toint(dimsizes(mme_value1(:,0,0,0))*sigh) - 1
mme_value(0,:,:,:) = mme_value1(n10,:,:,:)
mme_value(1,:,:,:) = mme_value1(n90,:,:,:)
;; NAJS.
fmna = addfile(mnadir,"r")
mna_value1 = fmna->m_value  ;(times,lei,:,:)
mna_value  = new((/2,2,31,13/),float)  ;(0.25/97.5,lei,:,:)
mna_value  = 0.
n10  = toint(dimsizes(mna_value1(:,0,0,0))*sigl) - 1
n90  = toint(dimsizes(mna_value1(:,0,0,0))*sigh) - 1
mna_value(0,:,:,:) = mna_value1(n10,:,:,:)
mna_value(1,:,:,:) = mna_value1(n90,:,:,:)

;; EASJ.
sigl = 0.1
sigh = 0.9

mea_valuet  = new((/2,2,31,13/),float)  ;(0.25/97.5,lei,:,:)
mea_valuet  = 0.
n10  = toint(dimsizes(mea_value1(:,0,0,0))*sigl) - 1
n90  = toint(dimsizes(mea_value1(:,0,0,0))*sigh) - 1
mea_valuet(0,:,:,:) = mea_value1(n10,:,:,:)
mea_valuet(1,:,:,:) = mea_value1(n90,:,:,:)
;; MESJ.
mme_valuet  = new((/2,2,31,13/),float)  ;(0.25/97.5,lei,:,:)
mme_valuet  = 0.
n10  = toint(dimsizes(mme_value1(:,0,0,0))*sigl) - 1
n90  = toint(dimsizes(mme_value1(:,0,0,0))*sigh) - 1
mme_valuet(0,:,:,:) = mme_value1(n10,:,:,:)
mme_valuet(1,:,:,:) = mme_value1(n90,:,:,:)
;; NAJS.
mna_valuet  = new((/2,2,31,13/),float)  ;(0.25/97.5,lei,:,:)
mna_valuet  = 0.
n10  = toint(dimsizes(mna_value1(:,0,0,0))*sigl) - 1
n90  = toint(dimsizes(mna_value1(:,0,0,0))*sigh) - 1
mna_valuet(0,:,:,:) = mna_value1(n10,:,:,:)
mna_valuet(1,:,:,:) = mna_value1(n90,:,:,:)
;-----------------------------------------------------------

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
;--------------------------------------------------------------
var1 = "jets"
var2 = "jetl"
var3 = "jetz"
;; NAJS.
jetsna1 = fna->$var1$ ; jet speed
jetlna1 = fna->$var2$ ; jet latitude
;jetzna1 = fna->$var3$ ; jet zonality

jetlna = ISO_Filter(jetlna1,time,201,nwaa,nwbb)
jetsna = ISO_Filter(jetsna1,time,201,nwaa,nwbb)
;jetzna = ISO_Filter(jetzna1,time,201,nwaa,nwbb)

;; EAJS.
jetsea1 = fea->$var1$ ; jet speed
jetlea1 = fea->$var2$ ; jet latitude

jetlea = ISO_Filter(jetlea1,time,201,nwaa,nwbb)
jetsea = ISO_Filter(jetsea1,time,201,nwaa,nwbb)
;; EAJS.
jetsme1 = fme->$var1$ ; jet speed
jetlme1 = fme->$var2$ ; jet latitude

jetlme = ISO_Filter(jetlme1,time,201,nwaa,nwbb)
jetsme = ISO_Filter(jetsme1,time,201,nwaa,nwbb)
;-----------------------------------------------------------------------------
;calculate phase-lag data. lag time 50 days
n = 30
plagna   = new((/2,n+1,13/),float) ; (indices,lag days, phases)
plagna   = 0.

plagea   = new((/2,n+1,13/),float) ; (indices,lag days, phases)
plagea   = 0.

plagme   = new((/2,n+1,13/),float) ; (indices,lag days, phases)
plagme   = 0.

;first conver from 12 to 3.
n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.31)
reamp   := reamp(n1:n2)
rephase := rephase(n1:n2)
jetlna    := jetlna(n1:n2)
jetsna    := jetsna(n1:n2)
;jetzna    := jetzna(n1:n2)

jetlea    := jetlea(n1:n2)
jetsea    := jetsea(n1:n2)

jetlme    := jetlme(n1:n2)
jetsme    := jetsme(n1:n2)

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
        plagna(0,j,i-1) = dim_avg(jetsna(nt1+j))
        plagna(1,j,i-1) = dim_avg(jetlna(nt1+j))
        ;plagna(2,j,i-1) = dim_avg(jetzna(nt1+j))
        plagea(0,j,i-1) = dim_avg(jetsea(nt1+j))
        plagea(1,j,i-1) = dim_avg(jetlea(nt1+j))
        plagme(0,j,i-1) = dim_avg(jetsme(nt1+j))
        plagme(1,j,i-1) = dim_avg(jetlme(nt1+j))
    end do
    delete([/nt,nt1/])
end do

do i = 0,4
    plagna(:,:,i+8) = plagna(:,:,i)
    plagea(:,:,i+8) = plagea(:,:,i)
    plagme(:,:,i+8) = plagme(:,:,i)
end do


x = ispan(0,n,1)  ; x labels
y = (/1,2,3,4,5,6,7,8,1,2,3,4,5/)
color1 = "MPL_BrBG";"MPL_RdYlBu"
color2 = "MPL_PuOr"
color3 = "MPL_PRGn"
undef("cmap")
cmap =  read_colormap_file("MPL_BrBG")
;cmap := cmap(::-1,:)
;label1
be = floor(max(plagna(0,:,:)))
en = ceil(min(plagna(0,:,:)))
ll = (/be,abs(en)/)
a  = max(ll)
if (be.eq.abs(en).and.a.eq.0) then
    ncolor = toint((1 + 1)/0.1 + 2)
    cmap1  = span_color_rgba(color1, ncolor)
    ;cmap1 := cmap1(::-1,:)
    label1 = fspan(-1., 1., ncolor-1)
else
    ncolor = toint((a + a)/0.2 + 2)
    cmap1  = span_color_rgba(color1, ncolor)
    ;cmap1 := cmap1(::-1,:)
    label1 = fspan(-a, a, ncolor-1)
end if
;label2
be = floor(max(plagna(1,:,:)))
en = ceil(min(plagna(1,:,:)))
ll = (/be,abs(en)/)
a  = 1.4
;if (be.eq.abs(en).and.a.eq.0) then
;    ncolor = toint((1 + 1)/0.1 + 2)
;    cmap2  = span_color_rgba(color2, ncolor)
;    label2 = fspan(-1., 1., ncolor-1)
;else
    ncolor = toint((a + a)/0.1 + 2)
    cmap2  = span_color_rgba(color2, ncolor)
    label2 = fspan(-a, a, ncolor-1)
;end if
;label3
;be = floor(max(plagna(2,:,:)))
;en = ceil(min(plagna(2,:,:)))
;ll = (/be,abs(en)/)
;a  = max(ll)
;if (be.eq.abs(en).and.a.eq.0) then
;    ncolor = toint((1 + 1)/0.1 + 2)
;    cmap3  = span_color_rgba(color3, ncolor)
;    label3 = fspan(-1., 1., ncolor-1)
;else
;    ncolor = toint((a + a)/0.2 + 2)
;    cmap3  = span_color_rgba(color3, ncolor)
;    label3 = fspan(-a, a, ncolor-1)
;end if
;; EASJ.
;label1
be = floor(max(plagea(0,:,:)))
en = ceil(min(plagea(0,:,:)))
ll = (/be,abs(en)/)
a  = max(ll)
if (be.eq.abs(en).and.a.eq.0) then
    ncolor = toint((1 + 1)/0.1 + 2)
    cmap4  = span_color_rgba(color1, ncolor)
    ;cmap4 := cmap4(::-1,:)
    label4 = fspan(-1., 1., ncolor-1)
else
    ncolor = toint((a + a)/0.2 + 2)
    cmap4  = span_color_rgba(color1, ncolor)
    ;cmap4 := cmap4(::-1,:)
    label4 = fspan(-a, a, ncolor-1)
end if
;label2
be = floor(max(plagea(1,:,:)))
en = ceil(min(plagea(1,:,:)))
ll = (/be,abs(en)/)
a  = max(ll)
if (be.eq.abs(en).and.a.eq.0) then
    ncolor = toint((1 + 1)/0.1 + 2)
    cmap5  = span_color_rgba(color2, ncolor)
    label5 = fspan(-1., 1., ncolor-1)
else
    ncolor = toint((a + a)/0.2 + 2)
    cmap5  = span_color_rgba(color2, ncolor)
    label5 = fspan(-a, a, ncolor-1)
end if
;; MEJS.
;label1
be = floor(max(plagme(0,:,:)))
en = ceil(min(plagme(0,:,:)))
ll = (/be,abs(en)/)
a  = max(ll)
if (be.eq.abs(en).and.a.eq.0) then
    ncolor = toint((1 + 1)/0.1 + 2)
    cmap6  = span_color_rgba(color1, ncolor)
    ;cmap6 := cmap6(::-1,:)
    label6 = fspan(-1., 1., ncolor-1)
else
    ncolor = toint((3 + 3)/0.2 + 2)
    cmap6  = span_color_rgba(color1, ncolor)
    ;cmap6 := cmap6(::-1,:)
    label6 = fspan(-3, 3, ncolor-1)
end if
;label2
be = floor(max(plagme(1,:,:)))
en = ceil(min(plagme(1,:,:)))
ll = (/be,abs(en)/)
a  = max(ll)
if (be.eq.abs(en).and.a.eq.0) then
    ncolor = toint((1 + 1)/0.1 + 2)
    cmap7  = span_color_rgba(color2, ncolor)
    label7 = fspan(-1., 1., ncolor-1)
else
    ncolor = toint((a + a)/0.1 + 2)
    cmap7  = span_color_rgba(color2, ncolor)
    label7 = fspan(-a, a, ncolor-1)
end if

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
;bres@tmYLMajorLengthF     = 0.01
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

; EAJS
plot1 = new(2, graphic)

k = 0

bres@vpXF  = 0.38
bres@vpYF  = 0.7
bres@tiYAxisString         = ""
bres@tiXAxisString         = "" 
bres@gsnLeftString = "~F21~(c) ~F21~JSI"
bres@gsnRightString = "~F21~EAJS"
plot1(k) = gsn_csm_blank_plot(wks,bres)
plot1(k) = add_polygon(wks,plot1(k),plagea(k,:,:),label4,cmap)
add_sigpoint(wks,plot1(k),plagea(k,:,:),mea_valuet(:,k,:,:),"black")
add_sigpoint(wks,plot1(k),plagea(k,:,:),mea_value(:,k,:,:),"red")
add_labelbar(wks,plot1(k),label4,cmap4)
print("plot d finished")

k = 1
bres@vpXF  = 0.38
bres@vpYF  = 0.45
bres@tiYAxisString         = ""
bres@tiXAxisString         = "~F21~Lag days" 
bres@gsnLeftString = "~F21~(f) ~F21~JLI"
plot1(k) = gsn_csm_blank_plot(wks,bres)
plot1(k) = add_polygon(wks,plot1(k),plagea(k,:,:),label5,color2)
add_sigpoint(wks,plot1(k),plagea(k,:,:),mea_valuet(:,k,:,:),"black")
add_sigpoint(wks,plot1(k),plagea(k,:,:),mea_value(:,k,:,:),"red")
add_labelbar(wks,plot1(k),label5,cmap5)
print("plot g finished")

;; MEJS
plot2 = new(2, graphic)

k = 0
bres@vpXF  = 0.05
bres@vpYF  = 0.7
bres@tiYAxisString         = "~F21~MJO phase"
bres@tiXAxisString         = "" 
bres@gsnLeftString = "~F21~(b) ~F21~JSI"
bres@gsnRightString = "~F21~MEJS"
plot2(k) = gsn_csm_blank_plot(wks,bres)
plot2(k) = add_polygon(wks,plot2(k),plagme(k,:,:),label6,cmap)
add_sigpoint(wks,plot2(k),plagme(k,:,:),mme_valuet(:,k,:,:),"black")
add_sigpoint(wks,plot2(k),plagme(k,:,:),mme_value(:,k,:,:),"red")
add_labelbar(wks,plot2(k),label6,cmap6)
print("plot c finished")

k = 1
bres@vpXF  = 0.05
bres@vpYF  = 0.45
bres@tiYAxisString         = "~F21~MJO phase"
bres@tiXAxisString         = "~F21~Lag days" 
bres@gsnLeftString = "~F21~(e) ~F21~JLI"
plot2(k) = gsn_csm_blank_plot(wks,bres)
plot2(k) = add_polygon(wks,plot2(k),plagme(k,:,:),label7,color2)
add_sigpoint(wks,plot2(k),plagme(k,:,:),mme_valuet(:,k,:,:),"black")
add_sigpoint(wks,plot2(k),plagme(k,:,:),mme_value(:,k,:,:),"red")
add_labelbar(wks,plot2(k),label7,cmap7)
print("plot f finished")

;; NAJS
plot3 = new(3, graphic)

k = 0
bres@vpXF  = 0.7
bres@vpYF  = 0.7
bres@tiYAxisString         = ""
bres@tiXAxisString         = "" 
bres@gsnLeftString = "~F21~(g) ~F21~JSI"
bres@gsnRightString = "~F21~NAJS"
plot3(k) = gsn_csm_blank_plot(wks,bres)
plot3(k) = add_polygon(wks,plot3(k),plagna(k,:,:),label1,cmap)
add_sigpoint(wks,plot3(k),plagna(k,:,:),mna_valuet(:,k,:,:),"black")
add_sigpoint(wks,plot3(k),plagna(k,:,:),mna_value(:,k,:,:),"red")
add_labelbar(wks,plot3(k),label1,cmap1)
print("plot b finished")

k = 1
bres@vpXF  = 0.7
bres@vpYF  = 0.45
bres@tiYAxisString         = ""
bres@tiXAxisString         = "~F21~Lag days" 
bres@gsnLeftString = "~F21~(e) ~F21~JLI"
plot3(k) = gsn_csm_blank_plot(wks,bres)
plot3(k) = add_polygon(wks,plot3(k),plagna(k,:,:),label2,color2)
add_sigpoint(wks,plot3(k),plagna(k,:,:),mna_valuet(:,k,:,:),"black")
add_sigpoint(wks,plot3(k),plagna(k,:,:),mna_value(:,k,:,:),"red")
add_labelbar(wks,plot3(k),label2,cmap2)
print("plot e finished")

;k = 2
;bres@vpXF  = 0.05
;bres@vpYF  = 0.23
;bres@tiYAxisString         = "~F21~MJO phase"
;bres@tiXAxisString         = "~F21~Lag days" 
;bres@gsnLeftString = "~F21~(h) ~F21~JZI"
;plot3(k) = gsn_csm_blank_plot(wks,bres)
;plot3(k) = add_polygon(wks,plot3(k),plagna(k,:,:),label3,color3)
;add_sigpoint(wks,plot3(k),plagna(k,:,:),mna_value(:,k,:,:))
;add_labelbar(wks,plot3(k),label3,cmap3)
;print("plot h finished")

  draw(plot1(0))
  draw(plot1(1))
  draw(plot2(0))
  draw(plot2(1))
  draw(plot3(0))
  draw(plot3(1))
;draw(plot3(2))

  i = 0
  do while(i.eq.0)
    i = 0
  end do