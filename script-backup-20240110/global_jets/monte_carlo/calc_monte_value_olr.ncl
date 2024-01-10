load "/home/sunming/data5/cuixy/script/ncl_func/monte_caro.ncl"
;--------------------------------------------------------
; draw temperature anomaly in shading, streamfunction
; anomaly in contour and OLR anomaly in contour.
; time: 2023-12-17
;--------------------------------------------------------
;-----------------------------------------------------------
ddir   = "/home/sunming/data5/cuixy/global_jets/data"
rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

udir   = ddir+"/"+"olr_daily_anom_1979-2022.nc"
outdir = ddir+"/"+"monte_carlo/olr_sig_phase15.nc"
latmin = -30
latmax =  50
;--------------------------------------------------------
; function calc_lagvalue.
;--------------------------------------------------------
function calc_lagvalue(in22,phaset,ampt,value)
local in22,phaset,ampt,value,DIMS,dim1,dim2,llvalue,\
ip,nt,ntt,nt1,nindd
begin
  DIMS = dimsizes(value)
  dim1 = DIMS(1)
  dim2 = DIMS(2)
  llvalue = new((/2,8,dim1,dim2/),float)
  llvalue = 0.
  iphase  = (/1,5/)
  do i = 0, 1
    ip = iphase(i)
    nt = ind(phaset.eq.ip.and.ampt.gt.1.)
    nt1 = in22(nt)
    do j = 0, 7
      ntt   = -15 + j*5
      nindd = nt1 + ntt
      llvalue(i,j,:,:) = dim_avg_n(value(nindd,:,:), 0)
      delete([/ntt,nindd/])
    end do
    delete([/ip,nt,nt1/])
  end do

  return llvalue
end
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------
fu   = addfile(udir,"r")
u    = fu->olr(:,{latmin:latmax},:)
print("complete data selling "+ systemfunc("date"))

time = u&time
TIME = cd_calendar(time, 0)
year = toint(TIME(:,0))
mon  = toint(TIME(:,1))
day  = toint(TIME(:,2))

;-------------------------------------------------------
;; MJO phase and amp.
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
u   := u(n1:n2,:,:)
reamp   := reamp(n1:n2)
rephase := rephase(n1:n2)
n_ndjfm  = ind(mont.eq.11.or.mont.eq.12.or.mont.eq.1.or.mont.eq.2.or.mont.eq.3)
u       := u(n_ndjfm,:,:)    ; NDJFM can replace by variable.
reamp   := reamp(n_ndjfm)
rephase := rephase(n_ndjfm)

delete([/time,TIME,year,mon,day/])

time = u&time
TIME = cd_calendar(time, 0)
year = toint( TIME(:,0) )
mon  = toint( TIME(:,1) )
day  = toint( TIME(:,2) )

print("NDJFM done. "+ systemfunc("date"))
;-----------------------------------------------------------------------------
;;calculate phase-lag data. lag time 31 days. NDJFM. for shuttle.
;;first sel data DJF.  

n_djf = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
in1   = ispan(0,dimsizes(time)-1,1)   ; NDJFM.
in2   = in1(n_djf)  ; DJF indices
amp   = reamp(n_djf)
phase = rephase(n_djf)

ind_time = ispan(0,dimsizes(time)-1,1) 
print("ind_time done. "+ systemfunc("date"))
;; store lead-lag values for plotting. (-15 -10 ---- 15 20)
;; (in22,phaset,ampt,value)
nrandom = 10000

m_value_2 = new((/nrandom,2,8,dimsizes(u&lat),dimsizes(u&lon)/),float)
m_value_2 = 0.
print("random calc u begin. "+ systemfunc("date"))
do k = 0, nrandom - 1
  rd_time = monte_carlo1(ind_time)
  ut   = u(rd_time,:,:)
  m_value_2(k,:,:,:,:) = calc_lagvalue(in2,phase,amp,ut)
  delete([/ut,rd_time/])
end do
print("random calc done. "+ systemfunc("date"))

; sort values
m_value2 = dim_pqsort_n(m_value_2, 2, 0)
copy_VarCoords(u(0,:,:), m_value_2(0,0,0,:,:))
system("/bin/rm -f "+outdir)
ncdf = addfile(outdir,"c")
ncdf->m_value = m_value_2
print("sort done. "+ systemfunc("date"))