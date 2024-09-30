load "/home/sunming/data5/cuixy/script/ncl_func/calc_rws.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"

function stats_frequency(timt,pha,ampt,is,ie,x,y,n)
local is,ie,x,y
begin
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ;event all days.
  nda  := ndtooned(nda)
  igood = ind(.not.ismissing(nda))
  nda  := nda(igood)  ;index for z.
  pt = pha(nda)
  at = ampt(nda)
  pt  := pt(ind(at.gt.1.))
  nda := nda(ind(at.gt.1.))
  delete(igood)
  ;select amp >1. and phase.
  rws_f = new((/8/),float)
  do m = 0,7
    np = m+1
    print("phase "+ np +" begin.")
    ip = ind(pt.eq.np)
    ip@_FillValue = -999
    if (any(ip.eq.ip@_FillValue).eq."Missing")
      continue
    else
      out  = new((/n,dimsizes(ip)/),float)
      do k = 0,dimsizes(ip)-1
        do j = 1,n
          r1 = x(nda(ip(k))+j)
          r2 = y(nda(ip(k))+j)
          if (r1.gt.0.and.r2.lt.0) then
            out(j-1,k) = 1
          else
            out(j-1,k) = 0
          end if
        end do
      end do
      out := dim_sum(out)
      rws_f(m) = 100.*dim_sum(out)/(dimsizes(ip)*n)
      delete(out)
    end if 
    delete(ip)  
    print("phase "+ np +" done.")  
  end do

  return rws_f
end
;----------------------------------------------------------------------
  lats = 10
  latn = 60
  lonw = 80
  lone = 240

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  zdir = ddir+"/"+"rws.nc"
;;jet events
  ndir = ddir+"/"+"IOnP_date_omi.csv"
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"
  idir = ddir+"/"+"IOiP_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
;-----------------------------------------------------------
; main code.
;-----------------------------------------------------------
  f = addfile(zdir,"r")
  z = f->rws(:,{lats:latn},{lonw:lone},0)
  printVarSummary(z)
  time = z&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.10.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.5.and.day.eq.31)

  mon := mon(n1:n2)
  z  := z(n1:n2,:,:)
  nt = ind(mon.eq.10.or.mon.eq.11.or.mon.eq.12.or.mon.eq.1\
           .or.mon.eq.2.or.mon.eq.3.or.mon.eq.4.or.mon.eq.5)
  z   := z(nt,:,:)
  delete([/TIME,time,year,mon,day,nt/])
;; ONDJFMAM
  time = z&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  timi = cd_calendar(time,-2)
;;OMI
  delim = ","
  rmmf  = asciiread(omidir,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
  phaset = str_get_field(rmmf,4,delim)
  ampt = str_get_field(rmmf,7,delim)

  iyear := toint(iyear(1::))
  imon  := toint(imon(1::))
  iday  := toint(iday(1::))
  phaset := tointeger(phaset(1::))
  ampt := tofloat(ampt(1::))

  phase = new((/dimsizes(time)/),integer)
  amp   = new((/dimsizes(time)/),float)
  phase = 0
  amp   = 0.

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i)\
             .and.iday.eq.day(i))
    phase(i) = phaset(nt)
    amp(i) = ampt(nt)
    delete(nt)
  end do
  delete([/phaset,ampt,iyear,imon,iday/])
  print("complete OMI selling "+ systemfunc("date"))
;----------------------------------------------------------------------
;;jet events (year,mon,day) (events,days,lat,lon)
  nf  = asciiread(idir,-1,"string")
  is  = str_get_field(nf,1,delim)
  ie  = str_get_field(nf,2,delim)
  ist = toint(is(1::))
  iet = toint(ie(1::))

  nf := asciiread(ndir,-1,"string")
  is := str_get_field(nf,1,delim)
  ie := str_get_field(nf,2,delim)

  ist := array_append_record(ist, toint(is(1::)), 0)
  iet := array_append_record(iet, toint(ie(1::)), 0)

  nf := asciiread(fdir,-1,"string")
  is := str_get_field(nf,1,delim)
  ie := str_get_field(nf,2,delim)

  isf = toint(is(1::))
  ief = toint(ie(1::))

  ist := array_append_record(ist, toint(is(1::)), 0)
  iet := array_append_record(iet, toint(ie(1::)), 0)

  nf := asciiread(sdir,-1,"string")
  is := str_get_field(nf,1,delim)
  ie := str_get_field(nf,2,delim) 

  iss = toint(is(1::))
  ies = toint(ie(1::))

  ist := array_append_record(ist, toint(is(1::)), 0)
  iet := array_append_record(iet, toint(ie(1::)), 0)

  delete([/ie,is/])

; rws averaged.
  latn = 30
  lats = 40
  lonw = 105
  lone = 120

  ;latn = 35
  ;lats = 43
  ;lonw = 105
  ;lone = 130
  rad  = 4.0*atan(1.0)/180.0
  z1 = z(:,{lats:latn},{lonw:lone})
  lat1 = z1&lat
  clat1 = cos(lat1*rad)
  rws1  = wgt_areaave_Wrap(z1, clat1, 1.0, 1)

  latn = 25;33
  lats = 35;43
  lonw = 135;140
  lone = 150;160

  ;latn = 33
  ;lats = 43
  ;lonw = 140
  ;lone = 160
  z2 = z(:,{lats:latn},{lonw:lone})
  lat2 = z2&lat
  clat2 = cos(lat2*rad)
  rws2  = wgt_areaave_Wrap(z2, clat2, 1.0, 1)  

;rws1>0 rws2<0.
  mt = 5
  rft = stats_frequency(timi,phase,amp,ist,iet,rws1,rws2,mt)
  rff = stats_frequency(timi,phase,amp,isf,ief,rws1,rws2,mt)
  rfs = stats_frequency(timi,phase,amp,iss,ies,rws1,rws2,mt)

  print(rft)
  print(rff)
  print(rfs)

  wks = gsn_open_wks("x11", "bar")

  res             = True  
  res@gsnDraw     = False
  res@gsnMaximize = False
  res@gsnFrame    = False
  res@gsnRightStringFontHeightF = 0.015
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnRightStringOrthogonalPosF = 0.0
  res@gsnLeftStringOrthogonalPosF  = 0.0

  res@tiXAxisString = "~F21~MJO phase"
  res@tiYAxisString = "~F21~[%]"
  res@tiYAxisFontHeightF = 0.015
  res@tiXAxisFontHeightF = 0.015

  res@tmYROn = False
  res@tmXTOn = False
  res@tmYLLabelFontHeightF = 0.012
  res@tmXBLabelFontHeightF = 0.012
  res@tmXBMode = "Explicit" 
  res@tmXBValues = (/1,2,3,4,5,6,7,8/)
  res@tmXBLabels = (/"1","2","3","4","5","6","7","8"/)

  ;res@vpXF = 0.18
  ;res@vpYF = 0.95   
  ;res@vpWidthF  = 0.25
  ;res@vpHeightF = 0.35

  res@trXMaxF = 9
  res@trXMinF = 0
  res@trYMinF = 10.
  res@trYMaxF = 80.
  res@gsnYRefLine  = 10.


  res@gsnXYBarChart = True
  res@gsnXYBarChartBarWidth = 0.25

  res@gsnLeftString = "(a)"
  res@gsnXYBarChartColors = (/"red3"/)  
  plot1 = gsn_csm_xy(wks,fspan(.875,7.875,8),rff,res)   
  res@gsnXYBarChartColors = (/"blue3"/)         
  plot2 = gsn_csm_xy(wks,fspan(1.125,8.125,8),rfs,res)  

  overlay(plot1,plot2)
  draw(plot1)


  i = 0 
  do while(i.eq.0)
    i = 0
  end do