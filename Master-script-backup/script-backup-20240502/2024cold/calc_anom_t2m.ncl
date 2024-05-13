;-----------------------------------------------------------------
;; this script is just used to filter values by removing seasonal
;; cycle. Lanczos filter.
;; the clim used a smoothing window to reduce the warming trend.
;; time: 2024-03-31  author: cui xiangyang
;-----------------------------------------------------------------
load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_3D.ncl"
;----------------------------------------------------------
; Function: single climatology
;----------------------------------------------------------
undef("single_clim")
function single_clim(nt, x,yr,mo,da, ymS, ymL)
local nt,x,ddd,yyyyddd,xClimDay,xAnom,year1,mon1,day1,xtmp,Xtmp
begin
  xtmp     = x(nt,:,:)
  year1    = yr(nt)
  mon1     = mo(nt)
  day1     = da(nt)
  ddd      = day_of_year(year1, mon1, day1)
  yyyyddd  = year1*1000 + ddd
  xClimDay = clmDayTLL(xtmp, yyyyddd)

  nd = ind(yr.ge.ymS.and.yr.le.ymL)
  year2     = yr(nd)
  mon2      = mo(nd)
  day2      = da(nd)
  ddd2      = day_of_year(year2, mon2, day2)
  yyyyddd2  = year2*1000 + ddd2
  Xtmp  = x(nd,:,:)
  xAnom = calcDayAnomTLL(Xtmp, yyyyddd2, xClimDay)

  return(xAnom)

end
;-----------------------------------------------------------------
;; Main code.
;-----------------------------------------------------------------
  var  = "t"

  nwa  = 90.
  nwb  = 10.

  yrStrt = 19800301
  yrLast = 20240326

  lats = 0
  latn = 60
  lonw = 60
  lone = 200

  ddir = "/home/sunming/data5/cuixy/2024cold/data"

  fin  = addfile(ddir+"/ncep_t2m_mean_daily_1980-2024.nc","r")

  ymd = cd_calendar(fin->time, -2)
  ny  = ind(ymd.ge.yrStrt.and.ymd.le.yrLast)

  x   = fin->$var$(ny,{lats:latn},{lonw:lone})

  time = x&time
  lat  = x&lat
  lon  = x&lon

  delete([/ymd,ny,yrStrt,yrLast/])
;--5-day running mean
  x := runave_n_Wrap(x, 5, -1, 0)

  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

;--New martrix contains anomalies.
  nt       = ind(mon.eq.1.or.mon.eq.2)
  year(nt) = year(nt) - 1
  xA = x
  delete(nt)
;--period 1. (1980-1995) (80/81-10/11)
  ymStrt1 = year(0)
  ymLast1 = ymStrt1 + 15
  nt11    = ind(year.ge.ymStrt1.and.year.le.(ymStrt1+29))    ;clim
  nt12    = ind(year.ge.ymStrt1.and.year.le.ymLast1)
  t2mAnom = single_clim(nt11,x,year,mon,day,ymStrt1,ymLast1)

  xA(nt12,:,:) = t2mAnom
  delete([/t2mAnom/])

;--period 2. (2008-2023) (94/95-23/24)
  ymLast2 = year(dimsizes(year)-1)
  ymStrt2 = ymLast2 - 15
  nt21    = ind(year.ge.(ymLast2-29).and.year.le.ymLast2)    ;clim
  nt22    = ind(year.ge.ymStrt2.and.year.le.ymLast2)    
  t2mAnom = single_clim(nt21,x,year,mon,day,ymStrt2,ymLast2)

  xA(nt22,:,:) = t2mAnom
  delete([/t2mAnom/])

;--period 3. (1996-2007)
  ymStrt = ymStrt1 
  ymLast = ymStrt + 29
  do i = ymLast1+1,ymStrt2-1
    ymStrt  = ymStrt + 1
    ymLast  = ymLast + 1

    nt31 = ind(year.ge.ymStrt.and.year.le.ymLast) ;clim
    t2mAnom  = single_clim(nt31,x,year,mon,day,i,i)
    
    nt32 = ind(year.eq.i)
    xA(nt32,:,:) = t2mAnom
    delete([/t2mAnom,nt31,nt32/])
  end do

  wts = filwgts_lanczos (41, 2, 1./nwa,1./nwb , 1.0)
  xa  = wgt_runave_n_Wrap(xA, wts, 0, 0)

  print("calc done. "+ systemfunc("date"))

  create_nc_file(xa,time,lat,lon,ddir+"/ncep_daily_filtered_anom_"+var+"_1980-2024.nc",var)
  create_nc_file(xA,time,lat,lon,ddir+"/ncep_daily_anom_"+var+"_1980-2024.nc",var)