;--------------------------------------------------------------
; 具体步骤
; 1 确定待选事件： P13的第一天和最后一天，之后判断是否至少8天大于1.
; 2 确定类型： 同样的操作挑选P68的端点，之后判断是否至少8天大于1.
; 是，则为传播事件，否则任何情况为不传播事件。
; 3 传播类型： 判断P13最后一天和P68第一天的间隔，快或者慢事件。
; 4 存储： P13的第一天和P68的最后一天。
; 注意： 这些事件总体是MJO事件，还有其他特征可能。
;--------------------------------------------------------------
  amp_cri  = 1.0 
  ds_cri   = 10
  ds_cri2  = 5
  f_cri   = 10
  s_cri   = 20
  nyear   = 43
  ndays1  = 220 ;; all days.
  ndays2  = 181 ;; NDJFM.
  ndays3  = 200 ;; NDJFMA.

  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy"
  udir = ddir+"/data/ERA5_z500_daily_mean_1979-2022_r2.5.nc"

  f    = addfile(udir,"r")
  time = f->time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
;-------------------------------------------------------
;; selling MJO.
  delim  = ","
  rmmf   = asciiread(rmmdir,-1,"string")
  iyear  = str_get_field(rmmf,1,delim)
  imon   = str_get_field(rmmf,2,delim)
  iday   = str_get_field(rmmf,3,delim)
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
  delete([/ampt,phaset,iyear,imon,iday/])
  print("complete RMM selling "+ systemfunc("date"))

;------------------------------------------------------------
  n1 = ind(year.eq.1979.and.mon.eq.9.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.5.and.day.eq.1)

  amp   := amp(n1:n2)
  phase := phase(n1:n2)
  mon   := mon(n1:n2)
  day   := day(n1:n2)
  year  := year(n1:n2)
  time  := time(n1:n2)

  nt     = ind(mon.eq.10.or.mon.eq.11.or.mon.eq.12\
               .or.mon.eq.1.or.mon.eq.2.or.mon.eq.3.or.mon.eq.4.or.mon.eq.5)
  amp   := amp(nt)
  phase := phase(nt)
  mon   := mon(nt)
  day   := day(nt)
  year  := year(nt)
  time  := time(nt)

  delete(nt)
  nt = ind(mon.ge.1.and.mon.le.5)
  year(nt) = year(nt) - 1 ;; year 1979-2021.
  delete(nt)

  year_rc = new((/nyear,ndays1/),integer)
  ;year_rc = 0
  mon_rc  = year_rc
  day_rc  = year_rc
;;P6-8
  year_ne = new((/nyear,ndays1/),integer)
  ;year_rc = 0
  mon_ne  = year_ne
  day_ne  = year_ne
  print("select candidate event begin.")
;;挑出来待选择事件
  do i = 0, nyear-1
    ;print("year "+ i+" begin.")
    pt = phase( ind(year.eq.i+1979) )
    at = amp( ind(year.eq.i+1979) )
    yt = year( ind(year.eq.i+1979) )
    mt = mon( ind(year.eq.i+1979) )
    dt = day( ind(year.eq.i+1979) )
    j = 0
    do while (j.le.ndays2) 
      if (pt(j).le.3) then ;在1-3之中
        k = 1
        pte = pt(j+k)
        do while (pte.le.3)  
          k = k + 1 
          pte = pt(j+k)
        end do
        ;判断是否有n天为MJO
        ate = at(j:j+k-1)
        nn = dimsizes(ind(ate.gt.amp_cri))
        if (nn.ge.ds_cri) then
          year_rc(i,j) = yt(j)
          mon_rc(i,j)  = mt(j)
          day_rc(i,j)  = dt(j)
          year_rc(i,j+1) = yt(j+k-1)
          mon_rc(i,j+1)  = mt(j+k-1)
          day_rc(i,j+1)  = dt(j+k-1)
        end if 
        j = j+k-1 + 1
        delete(ate)
      else ;不在P1-3之中，找下一个。
        j = j + 1
      end if
    end do
    ;delete([/pt,at,yt,mt,dt/])

    j = 0
    do while (j.le.ndays3) 
          ;;P6-8
      if (pt(j).ge.6) then ;在1-3之中
        k = 1
        pte = pt(j+k)
        do while (pte.ge.6)  
          k = k + 1 
          pte = pt(j+k)
        end do
        ;判断是否有n天为MJO
        ate = at(j:j+k-1)
        nn = dimsizes(ind(ate.gt.amp_cri))
        if (nn.ge.ds_cri2) then
          year_ne(i,j) = yt(j)
          mon_ne(i,j)  = mt(j)
          day_ne(i,j)  = dt(j)
          year_ne(i,j+1) = yt(j+k-1)
          mon_ne(i,j+1)  = mt(j+k-1)
          day_ne(i,j+1)  = dt(j+k-1)
        end if 
        j = j+k-1 + 1
        delete(ate)
      else ;不在P6-8之中，找下一个。
        j = j + 1
      end if
    end do
    delete([/pt,at,yt,mt,dt/])
    ;print("year "+ i+" done.")
  end do
;; reshape year, mon and day rc to a vector.
  year_rc := ndtooned(year_rc)
  mon_rc  := ndtooned(mon_rc)
  day_rc  := ndtooned(day_rc)

  igood = ind(.not.ismissing(year_rc))
  year_rc := year_rc(igood)
  mon_rc  := mon_rc(igood)
  day_rc  := day_rc(igood)

  nt = ind(mon_rc.ge.1.and.mon_rc.le.5)
  year_rc(nt) = year_rc(nt) + 1 ;; year 1979-2021.
  delete([/igood, nt/])

  year_ne := ndtooned(year_ne)
  mon_ne  := ndtooned(mon_ne)
  day_ne  := ndtooned(day_ne)

  igood = ind(.not.ismissing(year_ne))
  year_ne := year_ne(igood)
  mon_ne  := mon_ne(igood)
  day_ne  := day_ne(igood)

  nt = ind(mon_ne.ge.1.and.mon_ne.le.5)
  year_ne(nt) = year_ne(nt) + 1 ;; year 1979-2021.
  delete(igood)
  print("select candidate event done.")
  timep = year_rc*10000+mon_rc*100+day_rc
  timen = year_ne*10000+mon_ne*100+day_ne

;;convert to UTC date.
  yyyy  = year_rc              
  mm    = mon_rc
  dd    = day_rc                        
  hh    = dd                                 
  mn    = dd
  sc    = dd
  hh    = 0                                 
  mn    = 0
  sc    = 0                                          
  units = "days after 1979-01-01 00:00:00"                                        
  utcp  = cd_inv_calendar(yyyy,mm,dd,hh,mn,sc,units, 0)
  delete([/yyyy,mm,dd,hh,mn,sc/])

  yyyy  = year_ne              
  mm    = mon_ne
  dd    = day_ne                      
  hh    = dd                                 
  mn    = dd
  sc    = dd
  hh    = 0                                 
  mn    = 0
  sc    = 0                                          
  units = "days after 1979-01-01 00:00:00"                                         
  utcn  = cd_inv_calendar(yyyy,mm,dd,hh,mn,sc,units, 0)
  delete([/yyyy,mm,dd,hh,mn,sc/])
  delete(nt)
;; pairs.
  pend = timep(1::2)
  nsrt = timen(0::2)
  psrt = timep(0::2)
  nend = timen(1::2)

  pendu = tointeger(utcp(1::2))
  nsrtu = tointeger(utcn(0::2))
  psrtu = tointeger(utcp(0::2))
  nendu = tointeger(utcn(1::2))

  IOnP = new((/200,2/),integer)
  IOfP = IOnP
  IOsP = IOnP
  IOiP = IOnP
  nt = new((/1/),integer)
  nn = nt
  nt = 0
  do i = 0, dimsizes(pend)-1
    do j = 0, dimsizes(nsrt)-1
      pnn = nsrtu(j)-pendu(i)-1
      if (pnn.gt.0.and.pnn.le.f_cri) then
        IOfP(i,0) = psrt(i)
        IOfP(i,1) = nend(j)
        nn  = i
        nt := array_append_record(nt, nn, 0)
        break
      else if (pnn.ge.s_cri.and.pnn.lt.40) then
        IOsP(i,0) = psrt(i)
        IOsP(i,1) = nend(j)
        nn  = i
        nt := array_append_record(nt, nn, 0)
        break
      else if (pnn.lt.s_cri.and.pnn.gt.f_cri) then
        IOiP(i,0) = psrt(i)
        IOiP(i,1) = nend(j)
        nn  = i
        nt := array_append_record(nt, nn, 0) 
        break
      end if
      end if
      end if
    end do
  end do

  nt := nt(1:)
  ne  = ispan(0, dimsizes(pend)-1, 1)
  ;print(ne)
  ne@_FillValue = -999
  do i = 0, dimsizes(nt)-1
    ne = where(ne.eq.nt(i), ne@_FillValue, ne)
  end do
  ;print(ne)
  delete(nt)
  igood = ind(.not.ismissing(ne))
  ne := ne(igood)     ;;non-propagation.
;;
  nt = ind(mon.ge.1.and.mon.le.5)
  year(nt) = year(nt) + 1 ;; year 1979-2021.
  delete(nt)
  pend_ne = pend(ne)

  year_ne := pend_ne/10000
  mon_ne  := (pend_ne-year_ne*10000)/100
  day_ne  := (pend_ne-year_ne*10000)-mon_ne*100
  delete(igood)
  ;i = 0
  do i = 0, dimsizes(ne)-1
    nt = ind(year.eq.year_ne(i).and.mon.eq.mon_ne(i).and.day.eq.day_ne(i))
    ip = phase(nt+1)
    do while (ip.gt.3)
      nt = nt + 1
      ip = phase(nt)
    end do
    tss = year(nt-1)*10000+mon(nt-1)*100+day(nt-1)
    IOnP(i,0) = psrt(ne(i))
    IOnP(i,1) = tss
    delete([/nt,ip,tss/])
  end do

  IOnP := ndtooned(IOnP)
  IOsP := ndtooned(IOsP)
  IOfP := ndtooned(IOfP)
  IOiP := ndtooned(IOiP)

  ;print(dimsizes(pend))

  igood = ind(.not.ismissing(IOnP))
  IOnP := IOnP(igood)
  delete(igood)
  igood = ind(.not.ismissing(IOsP))
  IOsP := IOsP(igood)
  delete(igood)
  igood = ind(.not.ismissing(IOfP))
  IOfP := IOfP(igood)
  delete(igood)
  igood = ind(.not.ismissing(IOiP))
  IOiP := IOiP(igood)
  delete(igood)

  print(IOnP)
  print(IOfP)
  print(IOsP)
  print(IOiP)


  dq = str_get_dq()                
  fields = (/"start_date", "end_date"/)
  fields = dq + fields + dq                      
  header = [/str_join(fields,",")/]              
  format = "%d,%d" 

  head = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
;IOnP
  csv_filename = head+"/IOnP_date.csv"    
  alist  = [/IOnP(0::2),IOnP(1::2)/]       
  system("rm -rf " + csv_filename)               
  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format)  
;IOfP
  csv_filename = head+"/IOfP_date.csv"    
  alist  = [/IOfP(0::2),IOfP(1::2)/]       
  system("rm -rf " + csv_filename)               
  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format)  
;IOsP
  csv_filename = head+"/IOsP_date.csv"    
  alist  = [/IOsP(0::2),IOsP(1::2)/]       
  system("rm -rf " + csv_filename)               
  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format)  
;IOiP
  csv_filename = head+"/IOiP_date.csv"    
  alist  = [/IOiP(0::2),IOiP(1::2)/]       
  system("rm -rf " + csv_filename)               
  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format)  