load "/home/sunming/data5/cuixy/script/ncl_func/monte_caro.ncl"
function stats_frequency(timt,pha,ampt,is,ie,x,y)
local is,ie,x,y,iw
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
  ip = ind(pt.eq.7.or.pt.eq.6)
  ip@_FillValue = -999
  if (any(ip.eq.ip@_FillValue).eq."Missing")
    continue
  else
    out = new(dimsizes(ip),float)
    do k = 0,dimsizes(ip)-1
      r1 = x(nda(ip(k)))
      r2 = y(nda(ip(k)))
      if (r1.gt.0.and.r2.lt.0) then
        out(k) = 1
      else
        out(k) = 0
      end if
    end do
    iw = generate_sample_indices(dimsizes(out), 1)
    rws_f = dim_sum(out(iw))/dimsizes(ip)
  end if 

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

  nf := asciiread(fdir,-1,"string")
  is := str_get_field(nf,1,delim)
  ie := str_get_field(nf,2,delim)

  isf = toint(is(1::))
  ief = toint(ie(1::))

  nf := asciiread(sdir,-1,"string")
  is := str_get_field(nf,1,delim)
  ie := str_get_field(nf,2,delim) 

  iss = toint(is(1::))
  ies = toint(ie(1::))

  delete([/ie,is/])
; rws averaged.
  latn = 27;30
  lats = 37;40
  lonw = 105
  lone = 120

  rad  = 4.0*atan(1.0)/180.0
  z1 = z(:,{lats:latn},{lonw:lone})
  lat1 = z1&lat
  clat1 = cos(lat1*rad)
  rws1  = wgt_areaave_Wrap(z1, clat1, 1.0, 1)

  latn = 25
  lats = 35
  lonw = 130;135
  lone = 145;150

  z2 = z(:,{lats:latn},{lonw:lone})
  lat2 = z2&lat
  clat2 = cos(lat2*rad)
  rws2  = wgt_areaave_Wrap(z2, clat2, 1.0, 1)  

;rws1>0 rws2<0.
;----------------------------------------------------------------------
;计算PCC
  ind_time = ispan(0,dimsizes(time)-1,1) 
  nrand = 1000
  rf = new((/2,nrand/),float) ;(events,phases,days)
  rf = 0.
  print("random calc begin. "+ systemfunc("date"))
  rseed1 = toint(systemfunc(" date +%s"))
  rseed2 = toint((12345l*rseed1)%2147483398l)
  random_setallseed(rseed1, rseed2)
  do k = 0, nrand - 1
    ;random seed.
    rseed1 = toint(systemfunc(" date +%s"))
    rseed2 = toint((12345l*rseed1)%2147483398l)+toint(random_normal(100, 10, 1))
    random_setallseed(rseed1, rseed2)

    rf(0,k) = stats_frequency(timi,phase,amp,isf,ief,rws1,rws2) ;IOfP
    rf(1,k) = stats_frequency(timi,phase,amp,iss,ies,rws1,rws2) ;IOsP
  end do
  print("random calc done. "+ systemfunc("date"))
  print(rf(0,10:20))
;;store.
  print("store begin. "+ systemfunc("date"))

  odir = ddir+"/"+"boostrap_rws_omi_p67_region2.nc"
  system("/bin/rm -f "+odir)   
  fout = addfile(odir,"c")  

  fAtt = True           
  fAtt@creation_date = systemfunc ("date")

  fileattdef( fout, fAtt )            
  fout->rf = rf   

  print("store done. "+ systemfunc("date"))  