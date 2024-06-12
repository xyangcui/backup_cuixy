load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
load "/home/sunming/data5/cuixy/script/ncl_func/monte_caro.ncl"
;--------------------------------------------------------
; monte+carlo test method 1.
; randomly shifting original dateset in recycle form.
;--------------------------------------------------------
datadir = "/home/sunming/data5/cuixy/global_jets/data/jet_indices_NAJS_1.nc"
outdir  = "/home/sunming/data5/cuixy/global_jets/data/monte_carlo/jet_indices_NAJS_1.nc"
var1    = "jets"
var2    = "jetl"
;var3   = "jetz"

rmmdir  = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

nw  = 100
nb  = 20

;--------------------------------------------------------
; Function to filter.
;--------------------------------------------------------
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
;--------------------------------------------------------

;--------------------------------------------------------
; function calc_lagvalue.
;--------------------------------------------------------
function calc_lagvalue(in22,phaset,ampt,jl,js,jz)
local am,nt,nt1,phaset,ampt,jl,js,jz,in22
begin
  am = 1.0
  plag = new((/3,31,13/),float)
  plag = 0.
  do i = 1,8
    nt  = ind(phaset.eq.i.and.ampt.gt.am)
    nt1 = in22(nt)
    do j = 0 ,30, 1
      plag(0,j,i-1) = dim_avg(js(nt1+j))
      plag(1,j,i-1) = dim_avg(jl(nt1+j))
      plag(2,j,i-1) = dim_avg(jz(nt1+j))
    end do
    delete([/nt,nt1/])
  end do
  do i = 0,4
    plag(:,:,i+8) = plag(:,:,i)
  end do

  return plag
end
;--------------------------------------------------------
;--------------------------------------------------------
; function calc_lagvalue.
;--------------------------------------------------------
function calc_lagvalue2(in22,phaset,ampt,jl,js)
local am,nt,nt1,phaset,ampt,jl,js,in22
begin
  am = 1.0
  plag = new((/2,31,13/),float)
  plag = 0.
  do i = 1,8
    nt  = ind(phaset.eq.i.and.ampt.gt.am)
    nt1 = in22(nt)
    do j = 0 ,30, 1
      plag(0,j,i-1) = dim_avg(js(nt1+j))
      plag(1,j,i-1) = dim_avg(jl(nt1+j))
    end do
    delete([/nt,nt1/])
  end do
  do i = 0,4
    plag(:,:,i+8) = plag(:,:,i)
  end do

  return plag
end

;--------------------------------------------------------
; Main code.
;--------------------------------------------------------
;read in data and selected strong MJO days.

; jet indices plot.
  f = addfile(datadir,"r")

  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint( TIME(:,0) )
  mon  = toint( TIME(:,1) )
  day  = toint( TIME(:,2) )

  jets1 = f->$var1$ ; jet speed
  jetl1 = f->$var2$ ; jet latitude
  ;jetz1 = f->$var3$ ; jet zonality.

  jetl = ISO_Filter(jetl1,time,201,nw,nb)
  jets = ISO_Filter(jets1,time,201,nw,nb)
  ;jetz = ISO_Filter(jetz1,time,201,nw,nb)

; MJO phase and amp.
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

  delete([/amp,phase/])
  print("sort MJO done. "+systemfunc("date"))

;------------------------------------------------------------
; sort NDJFM.
  n1 = ind(year.eq.1979.and.mon.eq.10.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.4.and.day.eq.1)
  mont = mon(n1:n2)
  jetl  := jetl(n1:n2)  ; NDJFM can replace by variable.
  jets  := jets(n1:n2)  ; NDJFM can replace by variable.
  ;jetz  := jetz(n1:n2)  ; NDJFM can replace by variable.
  reamp := reamp(n1:n2)
  rephase := rephase(n1:n2)
  n_ndjfm = ind(mont.eq.11.or.mont.eq.12.or.mont.eq.1.or.mont.eq.2.or.mont.eq.3)
  jetl  := jetl(n_ndjfm)  ; NDJFM can replace by variable.
  jets  := jets(n_ndjfm)  ; NDJFM can replace by variable.
  ;jetz  := jetz(n_ndjfm)  ; NDJFM can replace by variable.
  reamp := reamp(n_ndjfm)
  rephase := rephase(n_ndjfm)

  delete([/time,TIME,year,mon,day/])

  time = jetl&time
  TIME = cd_calendar(time, 0)
  year = toint( TIME(:,0) )
  mon  = toint( TIME(:,1) )
  day  = toint( TIME(:,2) )

;-----------------------------------------------------------------------------
;;calculate phase-lag data. lag time 31 days. NDJFM. for shuttle.
;;first sel data DJF.  

  n_djf = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  in1   = ispan(0,dimsizes(time)-1,1)   ; NDJFM.
  in2   = in1(n_djf)  ; DJF indices
  amp   = reamp(n_djf)
  phase = rephase(n_djf)

  ind_time = ispan(0,dimsizes(time)-1,1)
  n = 30
  nrandom = 10000
  m_value = new((/nrandom,2,n+1,13/),float)
  m_value = 0.

;amp>1.0
  do k = 0, nrandom-1
    rd_time = monte_carlo1(ind_time)
    jetl_tmp = jetl(rd_time)
    jets_tmp = jets(rd_time)
    ;jetz_tmp = jetz(rd_time)
    ;m_value(k,:,:,:) = calc_lagvalue(in2,phase,amp,\
    ;                             jetl_tmp,jets_tmp,jetz_tmp)
    m_value(k,:,:,:) = calc_lagvalue2(in2,phase,amp,jetl_tmp,jets_tmp)
    delete([/rd_time,jetl_tmp,jets_tmp/]) ;,jetz_tmp
  end do
  print("random calc done. "+ systemfunc("date"))

; sort values
  m_value1 = dim_pqsort_n(m_value, 2, 0)

  system("/bin/rm -f "+outdir)
  ncdf = addfile(outdir,"c")
  ncdf->m_value = m_value
  print("sort done. "+ systemfunc("date"))
