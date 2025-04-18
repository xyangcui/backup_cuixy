

ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"

lats = 10
latn = 80
lonw = 100
lone = 240

  fu = addfile(ddir+"/ecmwf/ecmwf_pf_anom_u200_2003-2022_r1.5.nc","r")
  u = fu->u(:,:,{lats:latn},{lonw:lone})
  dims = dimsizes(u)
  ;eofs
  f = addfile(ddir+"/eof_u200.nc","r")
  eof = f->x(evn|0:1,lat|:,lon|:)
  eval= eof@eval

  ;calculate pcs
  rad    = 4.*atan(1.)/180.
  clat   = u&lat            
  clat   = sqrt( cos(rad*clat) )   
  ;weights
  xw     = u*conform(u, clat, 2)
  copy_VarMeta(u, xw)
  xw    := xw(lat|:,lon|:,lead_time|:,time|:)
  ts_u   = new((/2,dims(0),dims(1)/),double)
  do i =0,dims(0)-1
    wx = xw(:,:,:,i)
    eof_ts = eofunc_ts(wx, eof, False)
    ts_u(:,i,:) = eof_ts/conform(eof_ts, sqrt(eval(0:1)), 0)
    delete([/wx,eof_ts/])
  end do

  printVarSummary(ts_u)
  ;store pcs.
  valid_time = fu->valid_time(:,:)
  time = todouble(valid_time&time)
  valid_time&time := time

  system("/bin/rm -f "+ddir+"/ecmwf_pcs_u200.nc")

  fout = addfile(ddir+"/ecmwf_pcs_u200.nc","c")
  fout->x = ts_u
  fout->valid_time = valid_time