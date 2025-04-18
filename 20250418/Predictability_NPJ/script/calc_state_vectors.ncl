;---------------------------------------------------------------
;prepare for LIM state vector.
;input variable: sst heat sf200 sf850 
;OLR:    -15 - 15N    20EOFs 71%
;sf200:   10 - 90N    15EOFs 71%
;u200:    NPJ domain  8EOFs 71%
;sst:                 16EOFs
;
;output
;file1: state vectors (ntime, PCs)
;file2: EOFs of u200
;file3: EOFs of z200
;file4: EOFs of OLR
;---------------------------------------------------------------
nu  = 20   ;u200
ns  = 16   ;sf200
ns2 = 12   ;sf850
nt  = 20   ;sst
no  = 25   ;olr

ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"
;ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data/NCEP"

  DIR  = ddir+"/pcs_u200.nc"
  fout = addfile(DIR,"r")
  ts_u = fout->x(0:nu-1,:)

  DIR  = ddir+"/pcs_sf200.nc"
  fout = addfile(DIR,"r")
  ts_s = fout->x(0:ns-1,:)

  DIR  = ddir+"/pcs_sf850.nc"
  fout = addfile(DIR,"r")
  ts_s2 = fout->x(0:ns2-1,:)

  DIR  = ddir+"/pcs_sst.nc"
  fout = addfile(DIR,"r")
  ts_t = fout->x(0:nt-1,:)

  DIR = ddir+"/pcs_olr.nc"
  fout = addfile(DIR,"r")
  ts_o  = fout->x(0:no-1,:)

;pcvar
  DIR  = ddir+"/eof_u200.nc"
  fout = addfile(DIR,"r")
  eof = fout->x
  var = eof@pcvar
  print("u200 First "+nu+" var: " +sum(var(0:nu-1)))
  delete([/DIR,fout,eof,var/])

  DIR  = ddir+"/eof_sf200.nc"
  fout = addfile(DIR,"r")
  eof = fout->x
  var = eof@pcvar
  print("sf200 First "+ns+" var: " +sum(var(0:ns-1)))
  delete([/DIR,fout,eof,var/])

  DIR  = ddir+"/eof_sf850.nc"
  fout = addfile(DIR,"r")
  eof = fout->x
  var = eof@pcvar
  print("sf850 First "+ns2+" var: " +sum(var(0:ns2-1)))
  delete([/DIR,fout,eof,var/])

  DIR  = ddir+"/eof_sst.nc"
  fout = addfile(DIR,"r")
  eof = fout->x
  var = eof@pcvar
  print("sst First "+nt+" var: " +sum(var(0:nt-1)))
  delete([/DIR,fout,eof,var/])


  DIR = ddir+"/eof_olr.nc"
  fout = addfile(DIR,"r")
  eof = fout->x
  var = eof@pcvar
  print("olr First "+no+" var: " +sum(var(0:no-1)))
  delete([/DIR,fout,eof,var/])

;store data. 0: pcs 1: times
  dim1 = dimsizes(ts_u)
  dim2 = dimsizes(ts_s)
  dim3 = dimsizes(ts_s2)
  dim4 = dimsizes(ts_t)
  dim5 = dimsizes(ts_o)

  ;u s t o
  ;ts_tmp = new((/dim1(0)+dim2(0)+dim4(0)+dim5(0),dim1(1)/),typeof(ts_u))
  ts_tmp = new((/dim1(0)+dim2(0)+dim3(0)+dim4(0)+dim5(0),dim1(1)/),typeof(ts_u))
  ;ts_tmp = new((/dim1(0)+dim2(0)+dim3(0)+dim5(0),dim1(1)/),typeof(ts_u))
  nt1 = 0
  nt2 = nt1 + dim1(0) - 1
  ts_tmp(nt1:nt2,:) = ts_u

  nt1 = nt2 + 1
  nt2 = nt1 + dim2(0) - 1
  ts_tmp(nt1:nt2,:) = ts_s

  nt1 = nt2 + 1
  nt2 = nt1 + dim3(0) - 1
  ts_tmp(nt1:nt2,:) = ts_s2

  nt1 = nt2 + 1
  nt2 = nt1 + dim4(0) - 1
  ts_tmp(nt1:nt2,:) = ts_t

  nt1 = nt2 + 1
  nt2 = nt1 + dim5(0) - 1
  ts_tmp(nt1:nt2,:) = ts_o

  ;pcs
  DIR = ddir+"/state_vectors.nc"
  system("/bin/rm -f "+DIR)
  fout = addfile(DIR,"c")
  fout->x  = (/ts_tmp/)   