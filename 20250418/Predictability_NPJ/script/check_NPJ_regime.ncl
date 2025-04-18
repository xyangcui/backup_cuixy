;---------------------------------------------------------------
;prepare for LIM state vector.
;input variable: sst heat sf200 sf850 
;SST:   -20 - 15N 15EOFs 72%
;Heat:  -20 - 15N 15EOFs 39%
;sf200:  15 - 90N 12EOFs 70%
;sf850:  15 - 90N 10EOFs 70%
;
;output
;file1: state vectors (ntime, PCs)
;file2: EOFs of SST
;file3: EOFs of Heat
;file4: EOFs of sf200
;file5: EOFs of sf850
;---------------------------------------------------------------

;---------------------------------------------------------------
;Function EOF
;---------------------------------------------------------------
function calc_eof(xAnom,neof,optEOF,optETS)
local xAnom,neof,optEOF,optETS,rad,clat,xw,wx,eof,eof_ts,List
begin
  rad    = 4.*atan(1.)/180.
  clat   = xAnom&lat            
  clat   = sqrt( cos(rad*clat) )            

  ;weights
  xw     = xAnom*conform(xAnom, clat, 1)
  copy_VarMeta(xAnom, xw)
  ;reorder
  wx     = xw(lat|:,lon|:,time|:)                
  delete(xw)

  ;calcualte eof.
  eof    = eofunc_Wrap(wx, neof, optEOF)      
                                
  ;calc pc.
  eof_ts = eofunc_ts_Wrap(wx, eof, optETS)
  eof_ei = dim_stddev_n(eof_ts, 1)
  eof_t  = eof_ts/conform_dims(dimsizes(eof_ts), eof_ei, 0)
  ;eof_ts = dim_standardize_n( eof_ts, 0, 1)  

  eoft = eof;eof*conform_dims(dimsizes(eof), eof_ei, 0)
  copy_VarMeta(eof, eoft)
  copy_VarMeta(eof_ts, eof_t)

  pcvar = eof@pcvar
  print("pcvar: "+sum(pcvar))

  List = [/eof_t, eoft/]

  return List
end

latn = 80
lats = 10
lonw = 100
lone = 260

ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"
;----------------------------------------------------------------
;main code
;----------------------------------------------------------------
  f = addfile(ddir+"/ERA5_daily_anom_u200_1979-2022.nc","r")

  TIME=cd_calendar(f->time,0)
  timi=cd_calendar(f->time,-2)
  mon = TIME(:,1)
  nt = ind((timi.ge.19791001.and.timi.le.20240331)\
           .and.(mon.eq.11.or.mon.eq.12.or.mon.eq.1.or.mon.eq.2.or.mon.eq.3) )
  print(dimsizes(nt))

  x = f->u(nt,{lats:latn},{lonw:lone})
 
  delete([/TIME,timi,mon,nt/])

;EOF
  print("u200 begin.")
  List_u = calc_eof(x,2,True,False)
  print("u200 end.")
  eof_u = List_u[1]
  ts_u  = List_u[0]
  delete([/List_u,x/])

  ;pcs
  DIR = ddir+"/NPJ_pcs.nc"
  system("/bin/rm -f "+DIR)
  fout = addfile(DIR,"c")
  fout->x  = (/ts_u/)   

  ;eofs  
  DIR = ddir+"/eof_u.nc"
  system("/bin/rm -f "+DIR)
  fout = addfile(DIR,"c")
  fout->x  = eof_u