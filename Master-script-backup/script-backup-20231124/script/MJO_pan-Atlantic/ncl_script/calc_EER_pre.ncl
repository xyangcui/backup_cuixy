;----------------------------------------------------------
; calculate extreme event ratio.
; f0: reference ratio. 0.05
;----------------------------------------------------------
VAR        = "total_precipitation"
ftdir      = "/home/sunming/data5/cuixy/data/MJO_pan-Atlantic/"\
+"ERA5_s10/ERA5_s10_1979-2022.nc"
datadir    = "/home/sunming/data5/cuixy/data/MJO_pan-Atlantic"
mjodatedir = datadir+"/selected_MJO_dates.csv"
vardir     = datadir+"/ERA5_"+VAR+"/selected_precip_dates.nc"
outdir     = datadir+"/ERA5_"+VAR+"/precip_EER.nc"
f0         = 0.05
;----------------------------------------------------------

ft = addfile(ftdir,"r")
timet = ft->time
TIME  = cd_calendar(timet, 0)
year  = toint(TIME(:,0))
mon   = toint(TIME(:,1))
day   = toint(TIME(:,2))

fmjo  = asciiread(mjodatedir, -1, "string")
delim = ","
iyear = str_get_field(fmjo,1,delim)
imon  = str_get_field(fmjo,2,delim)
iday  = str_get_field(fmjo,3,delim)

iyear := toint(iyear(1::))
imon  := toint(imon(1::))
iday  := toint(iday(1::))

f  = addfile(vardir,"r")

cdate = f->cdate
lat   = cdate&lat
lon   = cdate&lon

;--MJO phases. central is phase 3.
phase_tmp = new(dimsizes(timet),integer)
phase_tmp = 0

do i = 0,dimsizes(iyear)-1
    nt = ind(year.eq.iyear(i).and.mon.eq.imon(i).and.day.eq.iday(i))
    phase_tmp(nt-2:nt+2)       = 3
    phase_tmp(nt-2+5:nt+2+5)   = 4
    phase_tmp(nt-2+10:nt+2+10) = 5
    phase_tmp(nt-2+15:nt+2+15) = 6
    phase_tmp(nt-2+20:nt+2+20) = 7
    phase_tmp(nt-2-15:nt+2-15) = 8
    phase_tmp(nt-2-10:nt+2-10) = 1
    phase_tmp(nt-2-5:nt+2-5)   = 2
    delete(nt)
end do

nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
phase = phase_tmp(nt)

delete([/nt,phase_tmp/])
;--cal-ratio
Rp = new((/8,dimsizes(lat),dimsizes(lon)/),float)
Rp = 0.
do i = 0,7
    ip = i + 1
    nt = ind(phase.eq.ip)
    fp = dim_sum_n(cdate(nt,:,:),0)*1./dimsizes(nt)
    Rp(i,:,:) = fp/f0
    delete([/nt,fp/])
end do
print(max(Rp))

Rp!0     = "phase"
Rp!1     = "lat"
Rp!2     = "lon"
Rp&lat   = lat
Rp&lon   = lon
Rp&phase = ispan(1,8,1)

system("rm -f "+outdir)

f1 = addfile(outdir,"c")

f1->Rp = Rp
