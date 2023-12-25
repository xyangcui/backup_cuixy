;----------------------------------------------------
;; draw vertical-latitude plot of temp.
;; author: cui xiangyang time: 2023-12-03
;----------------------------------------------------

ddir   = "/home/sunming/data5/cuixy/global_jets/data"
pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

phidir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022.nc"
idir   = ddir + "/ERA5_daily_anom_tlevel"
files  = systemfunc("ls " + idir + "/day_t_*.nc")

;---------------------------------------------------
; Main code.
;---------------------------------------------------
begin

  fphi = addfile(phidir,"r")
  ft   = addfiles(files,"r")
  ListSetType (ft, "cat")

  t   = dim_avg_n_Wrap(ft[:]->t(:,:,:,\
                                {lonmin:lonmax}),3)

  print("complete data selling "+ systemfunc("date"))

  time = fphi->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  t&time = time

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

;;-store lead-lag values for plotting. (-15 -10 ---- 15 25)

  llt = new((/8,dimsizes(t&level),dimsizes(t&latitude)/),float)
  copy_VarCoords(t(0,:,:), llt(0,:,:))
  llt = 0.

  do i = 0, 7
      ip = i + 1
      nt = ind(phase.eq.ip.and.amp.gt.1.)
      nind = ind_DJF(nt)
      llt(i,:,:) = dim_avg_n(t(nind,:,:), 0)
      delete([/ip,nt,nind/])
  end do

  print("complete calc leadlag data. "+ systemfunc("date"))

;; resources

;shading T plot
  cmap = read_colormap_file("MPL_RdBu")

  res = True
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@gsnMaximize = False
  res@gsnRightString     = "Avg: "+lonmin+":"+lonmax
  res@gsnRightStringFontHeightF = 0.01
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnMajorLatSpacing = 10

  res@vpHeightF   = 0.2
  res@vpWidthF    = 0.3
  res@tmYRMode    = "Automatic"

  res@cnLevelSelectionMode = "ManualLevels"       
  res@cnLevelSpacingF      = 0.2                  
  res@cnMinLevelValF       = -1.2               
  res@cnMaxLevelValF       = 1.2                
  res@cnLineLabelsOn       = False              
  res@cnFillOn             = True                
  res@cnFillPalette        = cmap(::-1,:)

  res@lbLabelBarOn   = False

  resP = True
  resP@gsnMaximize      = False
  resP@gsnFrame         = False
  resP@gsnDraw          = True
  resP@gsnPanelLabelBar = True
  resP@gsnPanelTop      = 0.98
  resP@gsnPanelBottom   = 0.02


  resP@pmLabelBarOrthogonalPosF = 0.0
  resP@pmLabelBarWidthF         = 0.15

  resP@lbOrientation      = "Vertical"
  resP@lbTitleString      = "Temp anom (K)"
  resP@lbTitlePosition    = "Bottom"
  ;resP@lbTitleDirection   = "Across"
  resP@lbTitleFontHeightF = 0.01
  resP@lbBoxLinesOn       = True
   ;resP@lbLabelStride     = 2
   ;resP@lbFillOpacityF    = 0.6
  resP@gsnPanelRowSpec    = True
  resP@lbLabelFontHeightF = 0.012
 
  print("complete res. "+ systemfunc("date"))

;; plot.

  wks = gsn_open_wks("eps", pdir+"/"+"t_ver_lat_"+region)

  label = (/"a ","b ","c ",\
           "d ","e ","f ","g ",\
           "h "/) + region + " phase" + ispan(1, 8, 1)


  resp = True
  resp@gsLineColor = "Green4"
  resp@gsLineDashPattern = 1
  resp@gsLineThicknessF  = 4.0
  resp@tfPolyDrawOrder   = "PostDraw"

  i = 0
  res@gsnLeftString = label(i)
  plot1 = gsn_csm_pres_hgt(wks, llt(i,:,:), res)
  str1  = unique_string("fill_lft")     
  plot1@$str1$ = gsn_add_polyline(wks,plot1,(/latmax,latmax/),(/1000,200/),resp)
  str2 = unique_string("fill_lft")     
  plot1@$str2$ = gsn_add_polyline(wks,plot1,(/latmin,latmin/),(/1000,200/),resp)

  i = 1
  res@gsnLeftString = label(i)
  plot2 = gsn_csm_pres_hgt(wks, llt(i,:,:), res)
  str1  = unique_string("fill_lft")     
  plot2@$str1$ = gsn_add_polyline(wks,plot2,(/latmax,latmax/),(/1000,200/),resp)
  str2 = unique_string("fill_lft")     
  plot2@$str2$ = gsn_add_polyline(wks,plot2,(/latmin,latmin/),(/1000,200/),resp)

  i = 2
  res@gsnLeftString = label(i)
  plot3 = gsn_csm_pres_hgt(wks, llt(i,:,:), res)
  str1 = unique_string("fill_lft")     
  plot3@$str1$ = gsn_add_polyline(wks,plot3,(/latmax,latmax/),(/1000,200/),resp)
  str2 = unique_string("fill_lft")     
  plot3@$str2$ = gsn_add_polyline(wks,plot3,(/latmin,latmin/),(/1000,200/),resp)

  i = 3
  res@gsnLeftString = label(i)
  plot4 = gsn_csm_pres_hgt(wks, llt(i,:,:), res)
  str1 = unique_string("fill_lft")     
  plot4@$str1$ = gsn_add_polyline(wks,plot4,(/latmax,latmax/),(/1000,200/),resp)
  str2 = unique_string("fill_lft")     
  plot4@$str2$ = gsn_add_polyline(wks,plot4,(/latmin,latmin/),(/1000,200/),resp)

  i = 4
  res@gsnLeftString = label(i)
  plot5 = gsn_csm_pres_hgt(wks, llt(i,:,:), res)
  str1 = unique_string("fill_lft")     
  plot5@$str1$ = gsn_add_polyline(wks,plot5,(/latmax,latmax/),(/1000,200/),resp)
  str2 = unique_string("fill_lft")     
  plot5@$str2$ = gsn_add_polyline(wks,plot5,(/latmin,latmin/),(/1000,200/),resp)

  i = 5
  res@gsnLeftString = label(i)
  plot6 = gsn_csm_pres_hgt(wks, llt(i,:,:), res)
  str1 = unique_string("fill_lft")     
  plot6@$str1$ = gsn_add_polyline(wks,plot6,(/latmax,latmax/),(/1000,200/),resp)
  str2 = unique_string("fill_lft")     
  plot6@$str2$ = gsn_add_polyline(wks,plot6,(/latmin,latmin/),(/1000,200/),resp)

  i = 6
  res@gsnLeftString = label(i)
  plot7 = gsn_csm_pres_hgt(wks, llt(i,:,:), res)
  str1 = unique_string("fill_lft")     
  plot7@$str1$ = gsn_add_polyline(wks,plot7,(/latmax,latmax/),(/1000,200/),resp)
  str2 = unique_string("fill_lft")     
  plot7@$str2$ = gsn_add_polyline(wks,plot7,(/latmin,latmin/),(/1000,200/),resp)

  i = 7
  res@gsnLeftString = label(i)
  plot8 = gsn_csm_pres_hgt(wks, llt(i,:,:), res)
  str1 = unique_string("fill_lft")     
  plot8@$str1$ = gsn_add_polyline(wks,plot8,(/latmax,latmax/),(/1000,200/),resp)
  str2 = unique_string("fill_lft")     
  plot8@$str2$ = gsn_add_polyline(wks,plot8,(/latmin,latmin/),(/1000,200/),resp)


  gsn_panel(wks, (/plot1,plot5,plot2,plot6,\
                plot3,plot7,plot4,plot8/), (/2,2,2,2/), resP)

  print("Draw done. "+ systemfunc("date"))

  ;i = 1
  ;do while (i.eq.1)
  ;  i = 1
  ;end do

end