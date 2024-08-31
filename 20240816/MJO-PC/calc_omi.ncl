;-----------------------------------------------------
;read in RMM indices.
;-----------------------------------------------------
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/OMI.csv"

  delim = ","
  rmmf  = asciiread(rmmdir,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
  omi1t = str_get_field(rmmf,5,delim)
  omi2t = str_get_field(rmmf,6,delim)
  ampt  = str_get_field(rmmf,7,delim)

  iyear := toint(iyear(1::))
  imon  := toint(imon(1::))
  iday  := toint(iday(1::))
  omi1  := tofloat(omi1t(1::))
  omi2  := tofloat(omi2t(1::))
  amp  := tofloat(ampt(1::))

  rmm1 = omi2
  rmm2 = -1*omi1
  print("complete RMM selling "+ systemfunc("date"))
;; calculate coamplitude and cophase.
  coPha = new(dimsizes(rmm1),integer)
  coPha = 0 
  cop = atan2(omi2, omi1)
  r2d = 45.0/atan(1.0)
  cop = cop *r2d
  pi  = atan(1.0)*4
;; convert to MJO phase.
  do i = 0,dimsizes(coPha)-1
      xt = cop(i)
      if (xt.ge.0.and.xt.lt.45.) then
        coPha(i) = 3
      else if (xt.ge.45..and.xt.lt.90.) then
        coPha(i) = 4
      else if (xt.ge.90..and.xt.lt.135.) then
        coPha(i) = 5
      else if (xt.ge.135..and.xt.lt.180.) then
        coPha(i) = 6
      else if (xt.ge.-180..and.xt.lt.-135.) then
        coPha(i) = 7
      else if (xt.ge.-135..and.xt.lt.-90.) then
        coPha(i) = 8
      else if (xt.ge.-90..and.xt.lt.-45.) then
        coPha(i) = 1
      else
        coPha(i) = 2
      end if
      end if
      end if
      end if
      end if
      end if
      end if
  end do
;;write csv.
  dq = str_get_dq()                
  fields = (/"year", "mon", "day", "phase", "OMI1", "OMI2", "Amp"/)
  fields = dq + fields + dq                      
  header = [/str_join(fields,",")/]              
  format = "%d,%d,%d,%d,%f,%f,%f" 

  head = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
;IOnP
  csv_filename = head    
  alist  = [/iyear,imon,iday,coPha,omi1,omi2,amp/]       
  system("rm -rf " + csv_filename)               
  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format) 