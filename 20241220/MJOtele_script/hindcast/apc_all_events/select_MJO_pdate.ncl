;--------------------------------------------------------------------------
;function: 返回位相对应的日期 格式：YYYYMMDD
;--------------------------------------------------------------------------
function return_p3(timt:integer,xdir,pht:integer,nphase:integer,ampt:float)
local delim,nint,nf,is,ie,nda,its,ite,ide,nda,pt,at,ip
begin
  delim = ","
  nint = ispan(0,dimsizes(timt),1) ; time index.
  nf = asciiread(xdir,-1,"string")
  is = str_get_field(nf,1,delim)
  ie = str_get_field(nf,2,delim)
  is := toint(is(1::))
  ie := toint(ie(1::))
  nda = new((/dimsizes(is),100/),integer)
  do i = 0,dimsizes(is)-1
    its = ind(timt.eq.is(i))
    ite = ind(timt.eq.ie(i))
    ide = ite-its
    nda(i,0:ide) = nint(its:ite)
  end do
  ;event all days.
  nda := ndtooned(nda)
  nda := nda(ind(.not.ismissing(nda)))  ;index for z.
  pt = pht(nda)
  at = ampt(nda)
  ;select amp >1. and phase.
  ip = ind(pt.eq.nphase.and.at.gt.1.)
  ip@_FillValue = -999
  if (any(ip.eq.ip@_FillValue).ne."Missing")
    nda := nda(ip)
  end if

  return timt(nda)

end
;----------------------------------------------------------------------
  lats = 20
  latn = 70
  lonw = 150
  lone = 240

  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"
  zdir = ddir+"/"+"ERA5_z500_daily_anom_1979-2022_r2.5.nc"
;;jet events
  fdir = ddir+"/"+"IO_date_omi.csv"

  omidir = "/home/sunming/data5/cuixy/data/clim_indices/OMIp.csv"
;-----------------------------------------------------------
; main code.
;-----------------------------------------------------------
;;times
  f = addfile(zdir,"r")
  z = f->z(:,{lats:latn},{lonw:lone})
  time = z&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))

  n1 = ind(year.eq.1979.and.mon.eq.10.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.5.and.day.eq.31)

  mon := mon(n1:n2)
  z := z(n1:n2,:,:)
  nt = ind(mon.eq.10.or.mon.eq.11.or.mon.eq.12.or.mon.eq.1\
           .or.mon.eq.2.or.mon.eq.3.or.mon.eq.4.or.mon.eq.5)
  z := z(nt,:,:)
  delete([/TIME,time,year,mon,day,nt/])
;; ONDJFMAM
  time = z&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  timi = cd_calendar(time,-2)
;;OMI
  delim = ","
  rmmf  = asciiread(omidir,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
  phaset = str_get_field(rmmf,4,delim)
  ampt = str_get_field(rmmf,7,delim)

  iyear := toint(iyear(1::))
  imon  := toint(imon(1::))
  iday  := toint(iday(1::))
  phaset := tointeger(phaset(1::))
  ampt := tofloat(ampt(1::))

  phase = new((/dimsizes(time)/),integer)
  amp   = new((/dimsizes(time)/),float)
  phase = 0
  amp   = 0.

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i)\
             .and.iday.eq.day(i))
    phase(i) = phaset(nt)
    amp(i) = ampt(nt)
    delete(nt)
  end do
  delete([/phaset,ampt,iyear,imon,iday/])
  print("complete OMI selling "+ systemfunc("date"))
;----------------------------------------------------------------------
;计算PCC
  print("IOfP begin. "+ systemfunc("date")) 
  n = 1
  IOfP_p1 = return_p3(timi,fdir,phase,n,amp) 
  n = 2
  IOfP_p2 = return_p3(timi,fdir,phase,n,amp) 
  n = 3
  IOfP_p3 = return_p3(timi,fdir,phase,n,amp) 
  n = 4
  IOfP_p4 = return_p3(timi,fdir,phase,n,amp) 
  n = 5
  IOfP_p5 = return_p3(timi,fdir,phase,n,amp) 
  n = 6
  IOfP_p6 = return_p3(timi,fdir,phase,n,amp) 
  n = 7
  IOfP_p7 = return_p3(timi,fdir,phase,n,amp) 
  n = 8
  IOfP_p8 = return_p3(timi,fdir,phase,n,amp)  
  print("IOfP done. "+ systemfunc("date"))

;;store.
  dq = str_get_dq()                
  fields = (/"date"/)
  fields = dq + fields + dq                      
  header = [/fields/]              
  format = "%d"
  head = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/pdate"

  m = 1
  ;IOfP
  csv_filename = head+"/IO_p"+m+"_a1.csv"    
  alist  = [/IOfP_p1/]       
  system("rm -rf " + csv_filename)               
  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format) 
  delete(alist) 

  m = 2
  ;IOfP
  csv_filename = head+"/IO_p"+m+"_a1.csv"    
  alist  = [/IOfP_p2/]       
  system("rm -rf " + csv_filename)               
  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format) 
  delete(alist) 

  m = 3
  ;IOfP
  csv_filename = head+"/IO_p"+m+"_a1.csv"    
  alist  = [/IOfP_p3/]       
  system("rm -rf " + csv_filename)               
  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format) 
  delete(alist) 

  m = 4
  ;IOfP
  csv_filename = head+"/IO_p"+m+"_a1.csv"    
  alist  = [/IOfP_p4/]       
  system("rm -rf " + csv_filename)               
  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format) 
  delete(alist) 

  m = 5
  ;IOfP
  csv_filename = head+"/IO_p"+m+"_a1.csv"    
  alist  = [/IOfP_p5/]       
  system("rm -rf " + csv_filename)               
  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format) 
  delete(alist) 

  m = 6
  ;IOfP
  csv_filename = head+"/IO_p"+m+"_a1.csv"    
  alist  = [/IOfP_p6/]       
  system("rm -rf " + csv_filename)               
  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format) 
  delete(alist) 

  m = 7
  ;IOfP
  csv_filename = head+"/IO_p"+m+"_a1.csv"    
  alist  = [/IOfP_p7/]       
  system("rm -rf " + csv_filename)               
  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format) 
  delete(alist) 

  m = 8
  ;IOfP
  csv_filename = head+"/IO_p"+m+"_a1.csv"    
  alist  = [/IOfP_p8/]       
  system("rm -rf " + csv_filename)               
  write_table(csv_filename, "w", header, "%s")   
  write_table(csv_filename, "a", alist, format) 
  delete(alist) 