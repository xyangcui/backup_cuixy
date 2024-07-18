;-----------------------------------------------------
; calculate mutiple linear regression by RMM indices.
; Y: 250-hPa u anomalies. x1: RMM1 x2: RMM2
; time: 20240528
;-----------------------------------------------------
  lats = 0
  latn = 90
  lonw = 0
  lone = 360
  
  ddir = "/home/sunming/data5/cuixy/global_jets/data"
  udir = ddir+"/u250_daily_anom_1979-2022_10-90_r2.5.nc"
;--------------------------------------------------------
; Main code.
;--------------------------------------------------------
  f = addfile(udir,"r")

  u = f->u(:,{lats:latn},{lonw:lone})
  lat = u&lat
  lon = u&lon 

  time = u&time
  TIME = cd_calendar(time, 0)
  year = toint(TIME(:,0))
  mon  = toint(TIME(:,1))
  day  = toint(TIME(:,2))
  print("complete data selling "+ systemfunc("date"))

  ;ut = u
  ;do i = 0, dimsizes(time) - 20
  ;  ut(i,:,:) = dim_avg_n(ut(i+3:i+7,:,:), 0)
  ;end do 
  ;u = ut
;-----------------------------------------------------
;read in RMM indices.
;-----------------------------------------------------
  rmmdir = "/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

  delim  = ","
  rmmf   = asciiread(rmmdir,-1,"string")
  iyear  = str_get_field(rmmf,1,delim)
  imon   = str_get_field(rmmf,2,delim)
  iday   = str_get_field(rmmf,3,delim)
  rmm1t  = str_get_field(rmmf,4,delim)
  rmm2t  = str_get_field(rmmf,5,delim)

  iyear  := toint(iyear(1::))
  imon   := toint(imon(1::))
  iday   := toint(iday(1::))
  rmm1t  := tofloat(rmm1t(1::))
  rmm2t  := tofloat(rmm2t(1::))

  rmm1 = new((/dimsizes(time)/),float)
  rmm2 = new((/dimsizes(time)/),float)
  rmm1 = 0.
  rmm2 = 0.

  do i = 0, dimsizes(time)-1
    nt = ind(iyear.eq.year(i).and.imon.eq.mon(i).and.iday.eq.day(i))
    rmm1(i) = rmm1t(nt)
    rmm2(i) = rmm2t(nt)
    delete(nt)
  end do

  nn = ispan(0, dimsizes(time)-1, 1)
  rmm1!0 = "timeind"
  rmm2!0 = "timeind"
  rmm1&timeind = nn
  rmm2&timeind = nn
  delete([/nn,rmm1t,rmm2t,iyear,imon,iday/])

  print("complete RMM selling "+ systemfunc("date"))
;----------------------------------------------------------
; prepare data for regression. (select DJF).
  n1 = ind(year.eq.1979.and.mon.eq.11.and.day.eq.1)
  n2 = ind(year.eq.2022.and.mon.eq.3.and.day.eq.1)
  mon  := mon(n1:n2)
  rmm1 := rmm1(n1:n2)
  rmm2 := rmm2(n1:n2)
  u    := u(n1:n2,:,:)

  nt = ind(mon.eq.12.or.mon.eq.1.or.mon.eq.2)
  rmm1 := rmm1(nt)
  rmm2 := rmm2(nt)
  u    := u(nt,:,:)

  ut    = u
  rmm1t = rmm1
  rmm2t = rmm2

  amp = sqrt(rmm1*rmm1+rmm2*rmm2)
  nt := ind(amp.gt.1.2)
  rmm1 := rmm1(nt)
  rmm2 := rmm2(nt)
  u    := u(nt,:,:)

  x = new((/dimsizes(nt),2/),float)
  x(:,0) = rmm1
  x(:,1) = rmm2
;; prepare an array for storing partial coefficents and R2.
  pc = new((/2,dimsizes(lat),dimsizes(lon)/),float)
  R2 = new((/dimsizes(lat),dimsizes(lon)/),float)
  pc = 0.
  R2 = 0.
  bb = new((/dimsizes(lat),dimsizes(lon)/),float)
  bb = 0.

  print("begin regression. "+ systemfunc("date"))
  do i = 0,dimsizes(lat)-1
    do j = 0,dimsizes(lon)-1
      b = reg_multlin_stats(u(:,i,j), x, False)
      pc(0,i,j) = b(1)
      pc(1,i,j) = b(2)
      bb(i,j) = b(0)
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
  copy_VarCoords(u(0,:,:), pc(0,:,:))

  system("/bin/rm -f "+ddir+"/MLRegression.nc")   
  fout = addfile(ddir+"/MLRegression.nc","c")  

  fAtt               = True           
  fAtt@creation_date = systemfunc ("date")
  fileattdef( fout, fAtt )            

  fout->coAmp = coAmp                                
  fout->coPha = coPha                       
  fout->R2 = R2
  fout->pc = pc

;; reconstruct u anomaly.
  print("begin reconstruct. "+ systemfunc("date"))

  Ru = ut
  Ru = 0.

  do i = 0, dimsizes(Ru&time)-1
    Ru(i,:,:) = bb + pc(0,:,:)*rmm1t(i) + pc(1,:,:)*rmm2t(i)
  end do

  copy_VarCoords(ut,Ru)

  Ru@long_name = "Reconstruct_u250 by MLRegression."

  system("/bin/rm -f "+ddir+"/Reconstruct_u250_anom_30-60.nc")   
  fout = addfile(ddir+"/Reconstruct_u250_anom_30-60.nc","c")  

  fAtt               = True           
  fAtt@creation_date = systemfunc ("date")
  fileattdef( fout, fAtt )            

  fout->u = Ru   

  print("finish reconstruct. "+ systemfunc("date"))                            