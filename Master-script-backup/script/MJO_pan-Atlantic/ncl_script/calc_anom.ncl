;----------------------------------------------------------
; calculate t2m anomalies in NA. 
; using a 30-year running window, function of detrending. 
; Information
; author: cui xiangyang time: 2023-10-29 release:v1
;----------------------------------------------------------
VAR = getenv("VAR")
datadir = "/home/sunming/data5/cuixy/data/MJO_pan-Atlantic/ERA5_"+VAR+"/"
inname  = "ERA5_"+VAR+"_1979-2022.nc"
ouname  = "ERA5_"+VAR+"_anom_1979-2022.nc"

var = getenv("var")
;----------------------------------------------------------
; Main code
;----------------------------------------------------------

f = addfile(datadir+inname,"r")

time = f->time
TIME = cd_calendar(time, 0)
year = toint(TIME(:,0))
mon  = toint(TIME(:,1))
day  = toint(TIME(:,2))
ddd      = day_of_year(year, mon, day)
yyyyddd  = year*1000 + ddd

x = f->$var$

xClimDay = clmDayTLL(x, yyyyddd)
; 15 day running.
xClimDay := runave_n_Wrap(xClimDay, 15, -1, 0)

xAnom = calcDayAnomTLL(x, yyyyddd, xClimDay)

copy_VarCoords(x, xAnom)


system("rm -f "+datadir+ouname)
f1  = addfile(datadir+ouname,"c")
f1->$var$ = xAnom
