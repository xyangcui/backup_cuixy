load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/calc_rws.ncl"

;------------------------------------------------------------
procedure create_nc_file(x1,x2,x3,TIME,LAT,LON,DIR:string,VAR1:string,\
                         VAR2,VAR3)
local x1,LON,LAT,TIME,DIR,VAR,fout
begin
    system("/bin/rm -f "+DIR)

  fout = addfile(DIR,"c")

  dimNames = (/"time", "lat", "lon"/)  
  dimSizes = (/ -1   ,  dimsizes(LAT),  dimsizes(LON)/) 
  dimUnlim = (/ True , False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout, "time" ,typeof(TIME),getvardims(TIME))                           
  filevardef(fout, "lat"  ,typeof(LAT),getvardims(LAT))                          
  filevardef(fout, "lon"  ,typeof(LON),getvardims(LON))                          
  filevardef(fout, VAR1    ,typeof(x1)  ,getvardims(x1))    
  filevardef(fout, VAR2    ,typeof(x2)  ,getvardims(x2)) 
  filevardef(fout, VAR3    ,typeof(x3)  ,getvardims(x3)) 
                                                                                             
  filevarattdef(fout,"time" ,TIME)                    
  filevarattdef(fout,"lat"  ,LAT)                   
  filevarattdef(fout,"lon"  ,LON)                     
  filevarattdef(fout,VAR1   ,x1)   
  filevarattdef(fout,VAR2   ,x2)  
  filevarattdef(fout,VAR3   ,x3)  

  setfileoption(fout,"DefineMode",False)

  fout->time  = (/TIME/)     
  fout->lat   = (/LAT/)
  fout->lon   = (/LON/) 
  fout->$VAR1$ = (/x1/)
  fout->$VAR2$ = (/x2/)
  fout->$VAR3$ = (/x3/)
  
end

;-------------------------------------------------------
;; draw rws.
;-------------------------------------------------------

  ddir  = "/home/sunming/data5/cuixy/global_jets/data"
  pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  udir  = ddir+"/"+"ERA5_daily_divu_1979-2022_r2.5.nc"
  vdir  = ddir+"/"+"ERA5_daily_divv_1979-2022_r2.5.nc"
  vodir = ddir+"/"+"ERA5_daily_vor_1979-2022_r2.5.nc"

  idir  = ddir+"/"+"ERA5_daily_anom_rws_1979-2022.nc"
  ldir  = ddir+"/"+"ERA5_daily_low_rws_1979-2022.nc"

  latmin = -20
  latmax =  70
  lonmin =  0
  lonmax =  220
;-------------------------------------------------------
;; main code.
;-------------------------------------------------------
  print("Prepare data begin. "+systemfunc("date"))
  fu  = addfile(udir,"r")
  fv  = addfile(vdir,"r")
  fvo = addfile(vodir,"r")

  u  = fu->u(:,{latmin:latmax},{lonmin:lonmax})
  v  = fv->v(:,{latmin:latmax},{lonmin:lonmax})
  vor= fvo->vor(:,{latmin:latmax},{lonmin:lonmax})

;; ISO filter.
  fi  = addfile(idir,"r")

  ui  = fi->u(:,{latmin:latmax},{lonmin:lonmax})
  vi  = fi->v(:,{latmin:latmax},{lonmin:lonmax})
  vori= fi->vor(:,{latmin:latmax},{lonmin:lonmax})

  ;ui   = iso_filter(u,201,100.,20.,0)
  ;vi   = iso_filter(v,201,100.,20.,0)
  ;vori = iso_filter(vor,201,100.,20.,0)
  ;copy_VarCoords(u, ui)
  ;copy_VarCoords(v, vi)
  ;copy_VarCoords(vor, vori)
  ;create_nc_file(ui,vi,vori,time,lat,lon,\
  ;               ddir+"/ERA5_daily_anom_rws_1979-2022.nc","u","v","vor")

;; Low filter.
  fl  = addfile(ldir,"r")

  ;ul  = fl->u(:,{latmin:latmax},{lonmin:lonmax})
  ;vl  = fl->v(:,{latmin:latmax},{lonmin:lonmax})
  ;vorl= fl->vor(:,{latmin:latmax},{lonmin:lonmax})

  ;wgt = filwgts_lancos(201, 0, 1/100., -999., 1) 
  ;ul  = wgt_runave_n_Wrap(u, wgt, 0, 0)
  ;vl  = wgt_runave_n_Wrap(v, wgt, 0, 0)
  ;vorl= wgt_runave_n_Wrap(vor, wgt, 0, 0)

  ;create_nc_file(ul,vl,vorl,time,lat,lon,\
  ;              ddir+"/ERA5_daily_low_rws_1979-2022.nc","u","v","vor")

  time = ui&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  lat  = ui&lat
  lon  = ui&lon
  print("Prepare data done. "+systemfunc("date"))

;--MJO date and DJF indices.
  print("sort MJO dates begin. "+systemfunc("date"))
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

  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.2.and.day.eq.28)

  amp   := amp(n1:n2)
  phase := phase(n1:n2)
  mon   := mon(n1:n2)
  time  := time(n1:n2)
  u    := u(n1:n2,:,:)
  v    := v(n1:n2,:,:)
  vor  := vor(n1:n2,:,:)

  nt     = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  amp   := amp(nt)
  phase := phase(nt)
  ul    = dim_avg_n_Wrap(u(nt,:,:),0)
  vl    = dim_avg_n_Wrap(v(nt,:,:),0)
  vorl  = dim_avg_n_Wrap(vor(nt,:,:),0)

  delete(nt)
  ind_DJF = amp&timeind   ; indice of DJF between 1979-2022 DJF.

  print("sort MJO dates done. "+systemfunc("date"))

;--calc rws and rws seperate.
  print("calc rws begin. "+systemfunc("date"))

  rws = new((/8,5,dimsizes(lat),dimsizes(lon)/),typeof(ui))
  rws = 0.
  rwse = new((/8,2,dimsizes(lat),dimsizes(lon)/),typeof(ui))
  rwse = 0.
  nt = ind(phase.eq.1.and.amp.gt.1.)
  nind = ind_DJF(nt)
  do j = 0, 7
    ntt   = -15 + j*5
    nindd = nind + ntt

    llui = dim_avg_n_Wrap(ui(nindd,:,:), 0)
    llvi = dim_avg_n_Wrap(vi(nindd,:,:), 0)
    llvoi= dim_avg_n_Wrap(vori(nindd,:,:), 0)

    llum = ul;dim_avg_n_Wrap(ul(nindd,:,:), 0)
    llvm = vl;dim_avg_n_Wrap(vl(nindd,:,:), 0)
    llvom= vorl;dim_avg_n_Wrap(vorl(nindd,:,:), 0)

    rws(j,:,:,:)  = calc_rws(llui,llum,llvi,llvm,llvoi,llvom)
    rwse(j,:,:,:) = calc_rwssep(llui,llvi,llvom)
    delete([/ntt,nindd,llui,llum,llvi,llvm,llvoi,llvom/])
  end do
  copy_VarCoords(ui(0,:,:), rws(0,0,:,:))
  copy_VarCoords(ui(0,:,:), rwse(0,0,:,:))
  printVarSummary(ui)
  printVarSummary(rws)
  print("calc rws done. "+systemfunc("date"))

;------------res------------------------------------------------

;shading T plot
  cmap = read_colormap_file("MPL_RdBu")

  res = True
  res@gsnFrame           = False
  res@gsnDraw            = False
  res@gsnMaximize        = False
  res@gsnAddCyclic       = False
  res@gsnMajorLatSpacing = 15
  res@gsnMajorLonSpacing = 30
  res@gsnLeftStringFontHeightF = 0.025 

  res@tmYLLabelFontHeightF  = 0.015
  res@tmXBLabelFontHeightF  = 0.015

  res@mpFillOn        = False
  res@mpFillDrawOrder = "Postdraw"
  res@mpMinLatF   = latmin
  res@mpMaxLatF   = latmax
  res@mpMinLonF   = lonmin
  ;res@mpCenterLonF= 180
  res@mpMaxLonF   = lonmax
  res@mpLandFillColor = "white"

  res@cnFillOn               = True
  res@cnFillPalette          = cmap(::-1,:)
  res@cnLinesOn              = False
  res@cnLineLabelsOn         = False
  res@cnLevelSelectionMode   = "Manuallevels"
  res@cnMinLevelValF         = -1.0
  res@cnMaxLevelValF         =  1.0
  res@cnLevelSpacingF        =  0.2
  ;res@cnMinLevelValF         = -0.5
  ;res@cnMaxLevelValF         =  0.5
  ;res@cnLevelSpacingF        =  0.1
  res@cnLabelBarEndStyle     = "ExcludeOuterBoxes"
  res@lbLabelBarOn           = False

  resP = True
  resP@gsnMaximize      = True
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02
  resP@gsnPanelYWhiteSpacePercent = 5

  resP@pmLabelBarWidthF         = 0.85
  resP@pmLabelBarOrthogonalPosF = -0.01

  resP@lbTitleString      = "unit: 10~S~-10~N~ s~S~-2~N~"
  resP@lbTitlePosition    = "Bottom"
  resP@lbTitleFontHeightF = 0.015
  resP@lbBoxLinesOn       = True
  resP@lbLabelFontHeightF = 0.012
  resP@gsnPanelRowSpec    = True

  respoly = True
  respoly@gsLineThicknessF = 3.0
  respoly@gsLineColor      = "black"

  print("prepare res done. "+ systemfunc("date"))

;---------------draw Fig.-------------------------------------

  wks = gsn_open_wks("x11", pdir+"rws_analysis.")

  ;zeta  = "~F33~z~F21~_~H-40~"
  zeta  = "~F21~_~H-18~~V-30~~F33~z~H23~"
  po    = "~F33~c"
  pi    = "~F34~6"
  term1 = "~V-30~~F21~-"+zeta+"~F34~Q~F34~w~F22~v'~B~"+po+"~N~"
  label = (/"~F21~(a) RWS:total","~F21~(b) RWS:term1","~F21~(c) RWS:term2",\
           "~F21~(d) RWS:term3","~F21~(e) RWS:term4"/)
  label2 = (/"~F21~(f) uadv",\
            "~F21~(g) vadv"/)
  plot1 = new(5,graphic)
  do i = 0,4
    res@gsnLeftString = label(i)
    plot1(i) = gsn_csm_contour_map(wks, rws(5,i,:,:), res)
    dum1 = gsn_add_polyline(wks, plot1(i),(/60,90,90,60,60/),\
     (/20,20,30,30,20/), respoly)
    dum2 = gsn_add_polyline(wks, plot1(i),(/140,170,170,140,140/),\
     (/25,25,35,35,25/), respoly)
    overlay(plot1(i),dum1)
    overlay(plot1(i),dum2)
  end do
  ;gsn_panel(wks, (/plot(0),plot(1),plot(2),plot(3),plot(4)/), (/1,2,2/), resP)

  plot2 = new(2,graphic)
  do i = 0,1
    res@gsnLeftString = label2(i)
    plot2(i) = gsn_csm_contour_map(wks, rwse(5,i,:,:), res)
    dum1 = gsn_add_polyline(wks, plot2(i),(/60,90,90,60,60/),\
     (/20,20,30,30,20/), respoly)
    dum2 = gsn_add_polyline(wks, plot2(i),(/140,170,170,140,140/),\
     (/25,25,35,35,25/), respoly)
    overlay(plot2(i),dum1)
    overlay(plot2(i),dum2)
  end do
  gsn_panel(wks, (/plot1(0),plot1(1),plot1(2),plot1(3),plot1(4),\
                  plot2(0),plot2(1)/), (/1,2,2,2/), resP)

  i = 0
  do while (i.eq.0)
    i = 0
  end do
