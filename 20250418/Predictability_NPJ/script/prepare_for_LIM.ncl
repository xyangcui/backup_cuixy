;---------------------------------------------------------------
;prepare for LIM state vector.
;input variable: sst heat sf200 sf850 
;OLR:   -15 - 15N   13EOFs 71%
;z200:   10 - 90N   17EOFs 71%
;u200:   NPJ domain 12EOFs 71%
;
;output
;file1: state vectors (ntime, PCs)
;file2: EOFs of u200
;file3: EOFs of z200
;file4: EOFs of OLR
;---------------------------------------------------------------

;---------------------------------------------------------------
;Function EOF
;---------------------------------------------------------------
function calc_eof(xAnom,neof,optEOF,optETS,dims)
local xAnom,neof,optEOF,optETS,rad,clat,xw,wx,eof,eof_ts,List
begin
  
  if (dims.eq.3) then
    rad    = 4.*atan(1.)/180.
    clat   = xAnom&lat            
    clat   = sqrt( cos(rad*clat) )            
    ;weights
    xw     = xAnom
    copy_VarMeta(xAnom, xw)
    ;reorder
    wx     = xw(lat|:,lon|:,time|:)                
    delete(xw)
  else
    wx     = xAnom(lon|:,time|:)
  end if
  ;standardize the variable.
  dimt = dimsizes(wx)
  ;wxt  = reshape(wx,(/dimt(0)*dimt(1),dimt(2)/))
  ;wxt  = dim_standardize_n(wxt,0,0)
  ;wxt := reshape(wxt,(/dimt(0),dimt(1),dimt(2)/))
  wxt  = wx*conform(wx, tofloat(clat), 0)
  copy_VarMeta(wx,wxt)
  ;calcualte eof.
  eof    = eofunc_Wrap(wxt, neof, optEOF)                                      
  ;calc pc.
  eof_ts = eofunc_ts_Wrap(wxt, eof, optETS)
  eof_ts = dim_standardize_n( eof_ts, 0, 1)  

  pcvar = eof@pcvar
  print("pcvar: "+sum(pcvar))

  List = [/eof_ts, eof/]

  return List
end


neofu = 100;25
neofz = 100;25
neoft = 100;25
neofo = 100;25

ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"
udir = ddir+"/ERA5_daily_anom_u200_1982-2022.nc"
s1dir = ddir+"/ERA5_daily_anom_sf200_1982-2022.nc"
s2dir = ddir+"/ERA5_daily_anom_sf850_1982-2022.nc"
odir = ddir+"/NOAA_daily_anom_olr_1982-2022.nc"
tdir = ddir+"/OISST_daily_anom_sst_1982-2022.nc"
;----------------------------------------------------------------
;main code
;----------------------------------------------------------------
  fu = addfile(udir,"r")
  fs1 = addfile(s1dir,"r")
  fs2 = addfile(s2dir,"r")
  fo = addfile(odir,"r")
  ft = addfile(tdir,"r")

  TIME=cd_calendar(fu->time,0)
  timi=cd_calendar(fu->time,-2)
  mon = TIME(:,1)
  nt = ind((timi.ge.19821201.and.timi.le.20220331)\
           .and.(mon.eq.12.or.mon.eq.1.or.mon.eq.2.or.mon.eq.3) )
  print(dimsizes(nt))

  u = fu->u(nt,{10:80},{100:240})
  s1 = fs1->sf(nt,{15:90},{0:360})
  s2 = fs2->sf(nt,{15:90},{0:360})
  t = ft->sst(nt,{15:60},{120:240})
  o = fo->olr(nt,{-20:15},{40:280})
  delete([/TIME,timi,mon,nt/])

;EOF
  print("u200 begin.")
  List_u = calc_eof(u,neofu,True,False,3)
  print("u200 end.")
  eof_u = List_u[1]
  ts_u  = List_u[0]
  delete([/List_u,u/])

  print("sf200 begin.")
  List_z = calc_eof(s1,neofz,True,False,3)
  print("sf200 end.")
  eof_z = List_z[1]
  ts_z  = List_z[0]
  delete([/List_z,s1/])

  print("sf850 begin.")
  List_s = calc_eof(s2,neofz,True,False,3)
  print("sf850 end.")
  eof_s2 = List_s[1]
  ts_s2  = List_s[0]
  delete([/List_s,s2/])

  print("sst begin.")
  List_t = calc_eof(t,neoft,True,False,3)
  print("sst end.")
  eof_t = List_t[1]
  ts_t  = List_t[0]
  delete([/List_t,t/])

  print("olr begin.")
  List_o = calc_eof(o,neofo,True,False,3)
  print("olr end.")
  eof_o = List_o[1]
  ts_o  = List_o[0]
  delete([/List_o,o/])

;store data. 0: pcs 1: times 
  DIR = ddir+"/pcs_u200.nc"
  system("/bin/rm -f "+DIR)
  fout = addfile(DIR,"c")
  fout->x  = (/ts_u/)  

  DIR = ddir+"/pcs_sf200.nc"
  system("/bin/rm -f "+DIR)
  fout = addfile(DIR,"c")
  fout->x  = (/ts_z/)  

  DIR = ddir+"/pcs_sf850.nc"
  system("/bin/rm -f "+DIR)
  fout = addfile(DIR,"c")
  fout->x  = (/ts_s2/) 

  DIR = ddir+"/pcs_sst.nc"
  system("/bin/rm -f "+DIR)
  fout = addfile(DIR,"c")
  fout->x  = (/ts_t/) 

  DIR = ddir+"/pcs_olr.nc"
  system("/bin/rm -f "+DIR)
  fout = addfile(DIR,"c")
  fout->x  = (/ts_o/)  
 

  ;eofs  
  DIR = ddir+"/eof_u200.nc"
  system("/bin/rm -f "+DIR)
  fout = addfile(DIR,"c")
  fout->x  = eof_u

  DIR = ddir+"/eof_sf200.nc"
  system("/bin/rm -f "+DIR)
  fout = addfile(DIR,"c")
  fout->x  = eof_z

  DIR = ddir+"/eof_sf850.nc"
  system("/bin/rm -f "+DIR)
  fout = addfile(DIR,"c")
  fout->x  = eof_s2

  DIR = ddir+"/eof_sst.nc"
  system("/bin/rm -f "+DIR)
  fout = addfile(DIR,"c")
  fout->x  = eof_t

  DIR = ddir+"/eof_olr.nc"
  system("/bin/rm -f "+DIR)
  fout = addfile(DIR,"c")
  fout->x  = eof_o
