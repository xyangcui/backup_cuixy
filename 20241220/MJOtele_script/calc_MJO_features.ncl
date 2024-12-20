
;--------------------------------------------------------------------------
;function: 返回位相3对应的日期 格式：YYYYMMDD
;--------------------------------------------------------------------------
function return_days(timt:integer,xdir)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,pt,at,ip
begin
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nf = asciiread(xdir,-1,"string")
  is = str_get_field(nf,1,delim)
  ie = str_get_field(nf,2,delim)
  is := toint(is(1::))
  ie := toint(ie(1::))
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ;event all days.
  nda := ndtooned(nda)
  nda := nda(ind(.not.ismissing(nda)))  ;index for z.
  return nda

end

;--------------------------------------------------------------------------
;procedure: 输出对应的天数：取MJO AMp大于1的天计数。
;--------------------------------------------------------------------------
procedure return_p3(timt:integer,xdir,amt:float,name:string)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,pt,at
begin
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nf = asciiread(xdir,-1,"string")
  is = str_get_field(nf,1,delim)
  ie = str_get_field(nf,2,delim)
  is := toint(is(1::))
  ie := toint(ie(1::))
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ni = dimsizes(is)
  ;event all days.
  nda := ndtooned(nda)
  nda := nda(ind(.not.ismissing(nda)))  ;index for z.
  at = amt(nda)
  print("event "+name+" :"+dimsizes(at)/ni)
end

;--------------------------------------------------------------------------
;function: 计算每个事件 每个位相的频次。
;--------------------------------------------------------------------------
function calc_pro(timt:integer,xdir,pht:integer,amt:float)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,igood,pt,at,ip,out,pct1
begin
  pct1 = new((/8/),float) 
  pct1 = 0.
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nf = asciiread(xdir,-1,"string")
  is = str_get_field(nf,1,delim)
  ie = str_get_field(nf,2,delim)
  is := toint(is(1::))
  ie := toint(ie(1::))
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ;event all days.
  nda := ndtooned(nda)
  nda := nda(ind(.not.ismissing(nda)))  ;index for z.
  pt = pht(nda)
  at = amt(nda)
  ;pt  := pt(ind(at.gt.1.))
  nie  = tofloat(dimsizes(ie))
  ;select amp >1. and phase.
  do m = 0, 7
    ip = ind(pt.eq.m+1)
    ip@_FillValue = -999
    if (any(ip.eq.ip@_FillValue).eq."Missing")
      pct1(m) = 0.
    else
      pct1(m) = dimsizes(ip)/nie
    end if
    delete(ip)
  end do
  return pct1
end
;-----------------------------------------------------------------
  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  odir = ddir+"/"+"olr_daily_anom_1979-2022.nc"
;;jet events
  ndir = ddir+"/"+"IOnP_date_omi.csv"
  fdir = ddir+"/"+"IOfP_date_omi.csv"
  sdir = ddir+"/"+"IOsP_date_omi.csv"
  idir = ddir+"/"+"IOiP_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
;-----------------------------------------------------------------
; main code.
;-----------------------------------------------------------------
  f = addfile(odir,"r")
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  timi = cd_calendar(time, -2) 
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

;;MJO events.
  freq = new((/4,8/),float)
  freq = 0.
;;jet events (year,mon,day) (events,days,lat,lon)
  print("frequency begin. "+ systemfunc("date"))
  freq(0,:) = calc_pro(timi,ndir,phase,amp) ;IOnP
  freq(1,:) = calc_pro(timi,idir,phase,amp) ;IOiP
  freq(2,:) = calc_pro(timi,fdir,phase,amp) ;IOfP
  freq(3,:) = calc_pro(timi,sdir,phase,amp) ;IOsP
  print("frequency done. "+ systemfunc("date"))

  IOnP = return_days(timi,ndir) ;IOnP
  IOiP = return_days(timi,idir) ;IOiP
  IOfP = return_days(timi,fdir) ;IOfP
  IOsP = return_days(timi,sdir) ;IOsP
;;jet events (year,mon,day) (events,days,lat,lon)
  return_p3(timi,ndir,amp,"IOnP") ;IOnP
  return_p3(timi,idir,amp,"IOiP") ;IOiP
  return_p3(timi,fdir,amp,"IOfP") ;IOfP
  return_p3(timi,sdir,amp,"IOsP") ;IOsP

  pn = amp(IOnP)
  pi = amp(IOiP)
  pf = amp(IOfP)
  ps = amp(IOsP)

  opt         = True
  opt@bin_min = 0.
  opt@bin_max = 4.

  pdfn = pdfx(pn,10,opt)
  pdfi = pdfx(pi,10,opt)
  pdff = pdfx(pf,10,opt)
  pdfs = pdfx(ps,10,opt)
  pdf = new((/4,dimsizes(pdfn)/),double)
  pdf(0,:) = pdfn
  pdf(1,:) = pdfi
  pdf(2,:) = pdff
  pdf(3,:) = pdfs
  pdf@bin_center = pdfn@bin_center


  wks  = gsn_open_wks ("x11","PDFX")
  res  = True
  res@gsnCenterString = "default 25 bins"
  res@xyLineColors = (/"red4","blue4","yellow4","black"/)
  res@xyLineThicknesses = (/4,4,4,4/)
  ;plot = gsn_csm_xy (wks, pdf@bin_center, pdf, res)

  resc  = True
  resc@gsnCenterString = "default 25 bins"
  resc@xyLineColors = (/"red4","blue4","yellow4","black"/)
  resc@xyLineThicknesses = (/4,4,4,4/)
  ;resc@tmXBFormat = 
  plot = gsn_csm_xy (wks, ispan(1,8,1), freq, resc)

  draw(plot)

  i = 0
  do while(i.eq.0)
    i = 0
  end do
