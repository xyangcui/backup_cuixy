;-----------------------------------------------------
; calculate mutiple linear regression by RMM indices.
; Y: 250-hPa u anomalies. x1: RMM1 x2: RMM2
; time: 20240528
;-----------------------------------------------------
  lats = 5
  latn = 60
  lonw = 0
  lone = 360
  
  ddir = "/home/sunming/data5/cuixy/ENSO-MJO-Jet"
  udir = ddir+"/data/sf250_daily_anom_1979-2022_r2.5.nc"
  ninodir = "/home/sunming/data5/cuixy/data/clim_indices/Nino_3.4_index.csv"
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------
  f = addfile(udir,"r")

  ul = f->sf(:,{lats:latn},{lonw:lone})
  u = ul*10^(-6)
  copy_VarCoords(ul,u)
  lat = u&lat
  lon = u&lon 

  time = u&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  print("complete data selling "+ systemfunc("date"))
;-----------------------------------------------------
;read in RMM indices.
;-----------------------------------------------------
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/OMI.csv"

  delim = ","
  rmmf  = asciiread(rmmdir,-1,"string")
  iyear = str_get_field(rmmf,1,delim)
  imon  = str_get_field(rmmf,2,delim)
  iday  = str_get_field(rmmf,3,delim)
  rmm1t = str_get_field(rmmf,6,delim)
  rmm2t = str_get_field(rmmf,5,delim)

  iyear := toint(iyear(1::))
  imon  := toint(imon(1::))
  iday  := toint(iday(1::))
  rmm1t := tofloat(rmm1t(1::))
  rmm2t := tofloat(rmm2t(1::))

  rmm1 = new((/dimsizes(time)/),float)
  rmm2 = new((/dimsizes(time)/),float)
  rmm1 = 0.
  rmm2 = 0.

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i)\
             .and.iday.eq.day(i))
    rmm1(i) = rmm1t(nt)
    rmm2(i) = -1*rmm2t(nt)
    delete(nt)
  end do

  delete([/rmm1t,rmm2t,iyear,imon,iday/])
  print("complete RMM selling "+ systemfunc("date"))

  ;ut = u
  ;do i = 0, dimsizes(time) - 20
  ;  ut(i,:,{0:240}) = dim_avg_n(ut(i+3:i+7,:,{0:240}), 0)
  ;  ut(i,:,{240:360}) = dim_avg_n(ut(i+3:i+7,:,{240:360}), 0)
  ;end do 
  ;u = ut
;----------------------------------------------------------
; prepare data for regression. (select DJF).
  n1 = ind(year.eq.1979.and.mon.eq.12.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)
  year := year(n1:n2)
  mon  := mon(n1:n2)
  rmm1 := rmm1(n1:n2)
  rmm2 := rmm2(n1:n2)
  u    := u(n1:n2,:,:)

  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  year := year(nt)
  mon  := mon(nt)
  rmm1 := rmm1(nt)
  rmm2 := rmm2(nt)
  u    := u(nt,:,:)

  ni = ind(mon.eq.12)
  year(ni) = year(ni)+1
;; select for EI, LA or NE. Dec year + 1.

  ;delete(iyear)
  delim  = ","
  ninof  = asciiread(ninodir,-1,"string")
  iyear  = str_get_field(ninof,1,delim)
  ievent = str_get_field(ninof,6,delim)

  iyear  := toint(iyear(1::))
  ievent := tofloat(ievent(1::))

;;EI
  nei = ind(ievent.eq.1)
  yei = iyear(nei)
;;LA
  nla = ind(ievent.eq.-1)
  yla = iyear(nla)
;;NE
  nne = ind(ievent.eq.0)
  yne = iyear(nne)

;; here for EI.
  if (nty.eq."NE") then
    ny = yne
    else if (nty.eq."EI")
      ny = yei
    else 
      ny = yla
    end if
  end if
  nc = ind(year.eq.ny(0))
  do i = 1, dimsizes(ny)-1
    ntemp = ind(year.eq.ny(i))
    nc := array_append_record(nc, ntemp, 0)
    delete(ntemp)
  end do

  x = new((/dimsizes(nc),2/),float)
  x(:,0) = rmm1(nc)
  x(:,1) = rmm2(nc)

  ampt = sqrt(x(:,0)^2+x(:,1)^2)
  amp  = avg(ampt)
;; prepare an array for storing partial coefficents and R2.
  pc = new((/2,dimsizes(lat),dimsizes(lon)/),float)
  R2 = new((/dimsizes(lat),dimsizes(lon)/),float)
  pc = 0.
  R2 = 0.

  print("begin regression. "+ systemfunc("date"))
  do i = 0,dimsizes(lat)-1
    do j = 0,dimsizes(lon)-1
      b = reg_multlin_stats(u(nc,i,j), x, False)
      pc(0,i,j) = b(1)
      pc(1,i,j) = b(2)
      R2(i,j) = b@SSR/b@SST
    end do
    delete(b)
  end do
  print("complete regression. "+ systemfunc("date"))
;; calculate coamplitude and cophase.
  coAmp = R2
  coPha = R2
  coAmp = sqrt(pc(0,:,:)*pc(0,:,:)+pc(1,:,:)*pc(1,:,:))
  coPha = atan2(pc(1,:,:), pc(0,:,:))

  r2d = 45.0/atan(1.0)
  coPha := coPha* r2d 
;; convert to MJO phase.
  do i = 0,dimsizes(lat)-1
    do j = 0,dimsizes(lon)-1
      xt = coPha(i,j)
      if (xt.ge.0.and.xt.lt.45.) then
        coPha(i,j) = 5.
      else if (xt.ge.45..and.xt.lt.90.) then
        coPha(i,j) = 6.
      else if (xt.ge.90..and.xt.lt.135.) then
        coPha(i,j) = 7.
      else if (xt.ge.135..and.xt.lt.180.) then
        coPha(i,j) = 8.
      else if (xt.ge.-180..and.xt.lt.-135.) then
        coPha(i,j) = 1.
      else if (xt.ge.-135..and.xt.lt.-90.) then
        coPha(i,j) = 2.
      else if (xt.ge.-90..and.xt.lt.-45.) then
        coPha(i,j) = 3.
      else
        coPha(i,j) = 4.
      end if
      end if
      end if
      end if
      end if
      end if
      end if
    end do
  end do

  R2@long_name = "coefficient of determination"
  coAmp@long_name = "coAmplitude with unit"
  coPha@long_name = "coPhase in radians"

  copy_VarCoords(u(0,:,:), coAmp)
  copy_VarCoords(u(0,:,:), coPha)
  copy_VarCoords(u(0,:,:), R2)

  system("/bin/rm -f "+ddir+"/data/MLRegression_"+nty+"_sf.nc")   
  fout = addfile(ddir+"/data/MLRegression_"+nty+"_sf.nc","c")  

  fAtt               = True           
  fAtt@creation_date = systemfunc ("date")
  fileattdef( fout, fAtt )            

  fout->coAmp = coAmp                                
  fout->coPha = coPha                       
  fout->R2 = R2