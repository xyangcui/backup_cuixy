load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"
;-------------------------------------------------------
; calculate OLR anomalies by subtracting daily mean 
; and the first 3 harmonics of climatology annual cycle.
; 
; ------------------------------------------------------

f = addfile("/home/sunming/data5/cuixy/DATA/NOAA-OLR/daily/raw_data/olr.day.mean_19740601-20221231.nc","r")

time = f->time
TIME = cd_calendar(time, 0)
year = toint(TIME(:,0))
nt   = ind(year.ge.1979)

olr  = f->olr(nt,:,:)

olranom = iso_filter(olr,201,70,20,1)

system ("rm -r "+"/home/sunming/data5/cuixy/data/MJO_pan-Atlantic/olr.daily.anom.1979-2022.nc")

fout = addfile("/home/sunming/data5/cuixy/data/MJO_pan-Atlantic/olr.daily.anom.1979-2022.nc","c")
fout ->olr = olranom
