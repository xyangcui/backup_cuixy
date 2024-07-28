;----------------------------------------
;; select Jet events by PC1.
;; "candidate": consecutive 1 pentad 1.0.
;; Pre 12 days must grow and decay.
;; remove each year pre 12 days and late.
;----------------------------------------

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

  ddir = "/home/sunming/data5/cuixy/global_jets/data/NAJ"
  indir = ddir+"/Ulevel_PC.csv"
  oudir = ddir+"/centered_date.csv"

;; function calc deria.
function deriv(pctt)
begin
  rv = new((/4/),float)
  rv = 0.
  do k = 1,4
    rv(k-1) = pctt(k)-pctt(k-1)
  end do
  return rv
end
; time table.
  idir = ddir + "/ERA5_daily_anom_ulevel_1979-2022_r2.5.nc"
  f = addfile(idir,"r")
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.4.and.day.eq.1)
  year := year(n1:n2)
  mon  := mon(n1:n2)
  day  := day(n1:n2)
  time := time(n1:n2)

  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  year := year(nt)
  mon  := mon(nt)
  day  := day(nt)
  time := time(nt)
  delete([/n1,n2,nt/])
;-----------------------------------------
;Main code.
;-----------------------------------------
  delim = ","
  rmmf  = asciiread(indir,-1,"string")
  pc  = str_get_field(rmmf,1,delim)
  pc := tofloat(pc)
  pc!0 = "time"
  pc&time = time

  re = new((/dimsizes(time),4/),float)
  re = 0.

; Dec year+1.
  ny = ind(mon.eq.12)
  yeart = year
  yeart(ny) = year(ny)+1
  delete(ny)

  j = 0
  do iy = 1980,2022
    print("year "+iy+" begin.")
    ny  = ind(yeart.eq.iy)  ; extract iy year.
    nb  = 5
    nl  = dimsizes(ny)-5
    pct = pc(ny)
    i = nb
    do while (i.lt.nl) 

      if (pct(i).gt.1.) then ;; positive event.
        tmp  = pct(i-2:i+2)  ; A pentad.
        tmin = tmp(minind(tmp)) ; min value.
        if (tmin.gt.1.) then
          imax = maxind(tmp)   ; index of max.
          if (imax.eq.2) then
            yearc = tofloat(year(ny))
            monc  = tofloat(mon(ny))
            dayc  = tofloat(day(ny))
            re(j,0) = yearc(i)
            re(j,1) = monc(i)
            re(j,2) = dayc(i)
            re(j,3) = pct(i)
            i = i+6
            j = j+1
            delete([/yearc,monc,dayc/])
          else
            i = i+1
          end if
        else
          i = i+1
        end if

      else if (pct(i).lt.-1.) then ;negative event.
        tmp  = pct(i-2:i+2)  ; A pentad.
        tmax = tmp(maxind(tmp)) ; min value.
        if (tmax.lt.-1.) then
          imin = minind(tmp)   ; index of max.
          if (imin.eq.2) then
            yearc = tofloat(year(ny))
            monc  = tofloat(mon(ny))
            dayc  = tofloat(day(ny))
            re(j,0) = yearc(i)
            re(j,1) = monc(i)
            re(j,2) = dayc(i)
            re(j,3) = pct(i)
            i = i+6
            j = j+1
            delete([/yearc,monc,dayc/])
          else
            i = i+1
          end if
        else
          i = i+1
        end if

      else
        i = i+1
      end if

      end if

    end do   ; do while.
  delete([/pct,ny/])
  print("year "+iy+" end.")
  end do     ; do.

;; write PC1 in CSV.
  system("rm -rf " + oudir)
  write_table(oudir,"w",[/re(:,0),re(:,1),re(:,2),re(:,3)/]\
              ,"%5.2f,%5.2f,%5.2f,%5.2f")    