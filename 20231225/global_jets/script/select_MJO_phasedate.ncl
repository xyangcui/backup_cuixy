;----------------------------------------------------
;; select MJO phases date.
;; every first date of phases is sorted.
;; return file: MJO_phase.csv
;; including MJO phases and corelated dates.
;; Phase(nphase,ndates)
;; author: cui xiangyang time: 2023-12-07
;----------------------------------------------------

ddir   = "/home/sunming/data5/cuixy/global_jets/data"
pdir   = "/home/sunming/data5/cuixy/global_jets/plot"
rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

phidir = ddir+"/"+"ERA5_daily_anom_phi250_1979-2022.nc"

;---------------------------------------------------
; Main code.
;---------------------------------------------------
begin

  fphi = addfile(phidir,"r")

  time = fphi->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

;-------------------------------------------------------
;; selling MJO.

  delim = ","
  rmmf  = asciiread(rmmdir,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
  phaset = str_get_field(rmmf,6,delim)
  ampt   = str_get_field(rmmf,7,delim)

  iyear  := toint(iyear(1::))
  imon   := toint(imon(1::))
  iday   := toint(iday(1::))
  phaset := toint(phaset(1::))
  ampt   := tofloat(ampt(1::))

  amp   = new((/dimsizes(time)/),float)
  phase = new((/dimsizes(time)/),integer)
  amp   = 0.
  phase = 0

  do i = 0, dimsizes(time)-1
      nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
      amp(i)   = ampt(nt)
      phase(i) = phaset(nt)
      delete(nt)
  end do

  nn = ispan(0, dimsizes(time)-1, 1)
  amp!0   = "timeind"
  phase!0 = "timeind"
  amp&timeind   = nn
  phase&timeind = nn
  delete([/nn,ampt,phaset,iyear,imon,iday/])

  print("complete RMM selling "+ systemfunc("date"))

;-------------------------------------------------------------

  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.2.and.day.eq.28)

  amp   := amp(n1:n2)
  phase := phase(n1:n2)
  mon   := mon(n1:n2)
  time  := time(n1:n2)

  nt     = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  amp   := amp(nt)
  phase := phase(nt)
  time  := time(nt)
  year  := year(nt)
  mon   := mon(nt)
  delete(nt)
  nt  = ind(amp.gt.1)
  print(dimsizes(nt))
  phase := phase(nt)
  time  := time(nt)
  year  := year(nt)
  mon   := mon(nt)
  delete(nt)
;-------------------------------------------------------------
; 8 phases. ()
  i = 0
  j = 1
  t = 0
  ndiff = 0.
  phase1 = new((/300/),integer,-999)
  phase_tmp = phase(ind(phase.eq.j))
  tind   = phase_tmp&timeind
  dtind  = dimsizes(tind)-2
  dn  = 0  ; length of each phase event.
  do while (i.le.dtind)
    ndiff = tind(i+1) - tind(i)
    if(ndiff.eq.1) then
      i  = i + 1
      dn = dn + 1
    else
      t  = t + 1       ; number of phase event.
      i  = i + 1
      phase1(t) = dn + 1
      dn = 0
    end if
  end do
  if (ndiff.gt.1) then
    t = t + 1
    phase1(t) = 1
  end if

  igood   = ind(.not.ismissing(phase1))
  phase1t = phase1(igood)
  delete([/phase_tmp,tind,dtind,igood/])
print("Phase1 finish. "+systemfunc("date"))

  i = 0
  j = 2
  t = 0
  ndiff = 0.
  phase2 = new((/300/),integer,-999)
  phase_tmp = phase(ind(phase.eq.j))
  tind   = phase_tmp&timeind
  dtind  = dimsizes(tind)-2
  dn  = 0  ; length of each phase event.
  do while (i.le.dtind)
    ndiff = tind(i+1) - tind(i)
    if(ndiff.eq.1) then
      i  = i + 1
      dn = dn + 1
    else
      t  = t + 1       ; number of phase event.
      i  = i + 1
      phase2(t) = dn + 1
      dn = 0
    end if
  end do
  if (ndiff.gt.1) then
    t = t + 1
    phase2(t) = 1
  end if

  igood   = ind(.not.ismissing(phase2))
  phase2t = phase2(igood)

  delete([/phase_tmp,tind,dtind,igood/])
print("Phase2 finish. "+systemfunc("date"))

  i = 0
  j = 3
  t = 0
  ndiff = 0.
  phase3 = new((/300/),integer,-999)
  phase_tmp = phase(ind(phase.eq.j))
  tind   = phase_tmp&timeind
  dtind  = dimsizes(tind)-2
  dn  = 0  ; length of each phase event.
  do while (i.le.dtind)
    ndiff = tind(i+1) - tind(i)
    if(ndiff.eq.1) then
      i  = i + 1
      dn = dn + 1
    else
      t  = t + 1       ; number of phase event.
      i  = i + 1
      phase3(t) = dn + 1
      dn = 0
    end if
  end do

  if (ndiff.gt.1) then
    t = t + 1
    phase3(t) = 1
  end if

  igood   = ind(.not.ismissing(phase3))
  phase3t = phase3(igood)

  delete([/phase_tmp,tind,dtind,igood/])
print("Phase3 finish. "+systemfunc("date"))

  i = 0
  j = 4
  t = 0
  ndiff = 0.
  phase4 = new((/300/),integer,-999)
  phase_tmp = phase(ind(phase.eq.j))
  tind   = phase_tmp&timeind
  dtind  = dimsizes(tind)-2
  dn  = 0  ; length of each phase event.
  do while (i.le.dtind)
    ndiff = tind(i+1) - tind(i)
    if(ndiff.eq.1) then
      i  = i + 1
      dn = dn + 1
    else
      t  = t + 1       ; number of phase event.
      i  = i + 1
      phase4(t) = dn + 1
      dn = 0
    end if
  end do

  if (ndiff.gt.1) then
    t = t + 1
    phase4(t) = 1
  end if

  igood   = ind(.not.ismissing(phase4))
  phase4t = phase4(igood)

  delete([/phase_tmp,tind,dtind,igood/])
print("Phase4 finish. "+systemfunc("date"))

  i = 0
  j = 5
  t = 0
  ndiff = 0.
  phase5 = new((/300/),integer,-999)
  phase_tmp = phase(ind(phase.eq.j))
  tind   = phase_tmp&timeind
  dtind  = dimsizes(tind)-2
  dn  = 0  ; length of each phase event.
  do while (i.le.dtind)
    ndiff = tind(i+1) - tind(i)
    if(ndiff.eq.1) then
      i  = i + 1
      dn = dn + 1
    else
      t  = t + 1       ; number of phase event.
      i  = i + 1
      phase5(t) = dn + 1
      dn = 0
    end if
  end do

  if (ndiff.gt.1) then
    t = t + 1
    phase5(t) = 1
  end if

  igood   = ind(.not.ismissing(phase5))
  phase5t = phase5(igood)

  delete([/phase_tmp,tind,dtind,igood/])
print("Phase5 finish. "+systemfunc("date"))

  i = 0
  j = 6
  t = 0
  ndiff = 0.
  phase6 = new((/300/),integer,-999)
  phase_tmp = phase(ind(phase.eq.j))
  tind   = phase_tmp&timeind
  dtind  = dimsizes(tind)-2
  dn  = 0  ; length of each phase event.
  do while (i.le.dtind)
    ndiff = tind(i+1) - tind(i)
    if(ndiff.eq.1) then
      i  = i + 1
      dn = dn + 1
    else
      t  = t + 1       ; number of phase event.
      i  = i + 1
      phase6(t) = dn + 1
      dn = 0
    end if
  end do

  if (ndiff.gt.1) then
    t = t + 1
    phase6(t) = 1
  end if

  igood   = ind(.not.ismissing(phase6))
  phase6t = phase6(igood)

  delete([/phase_tmp,tind,dtind,igood/])
print("Phase6 finish. "+systemfunc("date"))

  i = 0
  j = 7
  t = 0
  ndiff = 0.
  phase7 = new((/300/),integer,-999)
  phase_tmp = phase(ind(phase.eq.j))
  tind   = phase_tmp&timeind
  dtind  = dimsizes(tind)-2
  dn  = 0  ; length of each phase event.
  do while (i.le.dtind)
    ndiff = tind(i+1) - tind(i)
    if(ndiff.eq.1) then
      i  = i + 1
      dn = dn + 1
    else
      t  = t + 1       ; number of phase event.
      i  = i + 1
      phase7(t) = dn + 1
      dn = 0
    end if
  end do

  if (ndiff.gt.1) then
    t = t + 1
    phase1(t) = 1
  end if
 
  igood   = ind(.not.ismissing(phase7))
  phase7t = phase7(igood)

  delete([/phase_tmp,tind,dtind,igood/])
print("Phase7 finish. "+systemfunc("date"))

  i = 0
  j = 8
  t = 0
  ndiff = 0.
  phase8 = new((/300/),integer,-999)
  phase_tmp = phase(ind(phase.eq.j))
  tind   = phase_tmp&timeind
  dtind  = dimsizes(tind)-2
  dn  = 0  ; length of each phase event.
  do while (i.le.dtind)
    ndiff = tind(i+1) - tind(i)
    if(ndiff.eq.1) then
      i  = i + 1
      dn = dn + 1
    else
      t  = t + 1       ; number of phase event.
      i  = i + 1
      phase8(t) = dn + 1
      dn = 0
    end if
  end do

  if (ndiff.gt.1) then
    t = t + 1
    phase8(t) = 1
  end if
 
  igood   = ind(.not.ismissing(phase8))
  phase8t = phase8(igood)

  delete([/phase_tmp,tind,dtind,igood/])
print("Phase8 finish. "+systemfunc("date"))
;---------------------------------------------------------
print("---------------------------------------------------")
print("write .csv file begin. "+systemfunc("date"))

  csv_filename = ddir+"/mjo_phase_date.csv"
  dq = str_get_dq()

  fields = (/"phase1","phase2","phase3","phase4","phase5",\
            "phase6","phase7","phase8"/)
  fields = dq + fields + dq                    
  header = [/str_join(fields,",")/]              
                                                
  alist  = [/phase1t,phase2t,phase3t,phase4t,phase5t,phase6t,phase7t,phase8t/]  

  format = "%d,%d,%d,%d,%d,%d,%d,%d"          

  system("rm -rf " + csv_filename)               

  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format)
print("write .csv file done. "+ systemfunc("date"))

end


      
