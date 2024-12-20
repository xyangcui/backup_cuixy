load "/home/sunming/data5/cuixy/script/ncl_func/calc_rws.ncl"
;----------------------------------------------------------
; calc Rossby wave source.
; time: 2024-09-04  authour: cui xiangyang
;----------------------------------------------------------

lats = 10
latn = 70
lonw = 10
lone = 250

ddir  = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"

umdir = ddir+"/ERA5_u200_daily_mean_1979-2022_r2.5.nc"
vmdir = ddir+"/ERA5_v200_daily_mean_1979-2022_r2.5.nc"

uidir = ddir+"/ERA5_u200_daily_anom_1979-2022_r2.5.nc"
vidir = ddir+"/v200_daily_anom_1979-2022_r2.5.nc"

rmdir = ddir+"/ERA5_rwsv_daily_mean_1979-2022_r2.5.nc"
ridir = ddir+"/ERA5_rwsv_daily_anom_1979-2022_r2.5.nc"
;------------------------------------------------------------
; Main code
;------------------------------------------------------------

;calc div u and v and relative vorticity.
  wts = filwgts_lanczos(201, 2, 1/100., 1/20., 1.0) 

;  fum = addfile(umdir,"r")
;  fvm = addfile(vmdir,"r")

;  ut = fum->u(:,{200},{-90:90},:)
;  vt = fvm->v(:,{200},{-90:90},:)

;  time = ut&time
;  lat  = ut&lat
;  lon  = ut&lon

;  TIME = cd_calendar(time, 0)
;  year = toint(TIME(:,0))
;  mon  = toint(TIME(:,1))
;  day  = toint(TIME(:,2))
;  ddd     = day_of_year(year, mon, day) 
;  yyyyddd = year*1000 + ddd     

;  n1 = ind(year.eq.1979.and.mon.eq.10.and.day.eq.1)
;  n2 = ind(year.eq.2022.and.mon.eq.5.and.day.eq.31)

;  mon := mon(n1:n2)
;  um := ut(n1:n2,:,:)
;  vm := vt(n1:n2,:,:)
;  nt = ind(mon.eq.11.or.mon.eq.12.or.mon.eq.1.or.mon.eq.2.or.mon.eq.3)
;  um := um(nt,:,:)
;  vm := vm(nt,:,:)

;  um := dim_avg_n_Wrap(um, 0)
;  vm := dim_avg_n_Wrap(vm, 0)

;  print("calc mean begin.")
;  vor = uv2vrF_Wrap(um,vm)

;  dv  = uv2dvF_Wrap(um,vm) 
;  uvd = dv2uvF_Wrap(dv)
;  delete(dv)

;  u = uvd(0,:,:)
;  v = uvd(1,:,:)
;  copy_VarCoords(um,u)
;  copy_VarCoords(vm,v)
;  delete([/uvd,um,vm/])

;  ui = wgt_runave_n_Wrap(ut, wts, 0, 0)
;  vi = wgt_runave_n_Wrap(vt, wts, 0, 0)
;  print("calc anom begin.")
;  xa = uv2vrF_Wrap(ui,vi)

;  dv  = uv2dvF_Wrap(ui,vi) 
;  uvd = dv2uvF_Wrap(dv)
;  delete(dv)

;  uc = uvd(0,:,:,:)
;  vc = uvd(1,:,:,:)
;  copy_VarCoords(ui,uc)
;  copy_VarCoords(vi,vc)
;  delete(uvd)

;  system("/bin/rm -f " + rmdir)    
 ; system("/bin/rm -f " + ridir)  
;
;  fout1 = addfile(rmdir, "c")
;  fout2 = addfile(ridir, "c")

;  setfileoption(fout1,"DefineMode",True)
;  setfileoption(fout2,"DefineMode",True)

;  dimNames1 = (/"lat", "lon"/)  
;  dimSizes1 = (/dimsizes(lat),  dimsizes(lon)/) 
;  dimUnlim1 = (/ False, False/)   
;  filedimdef(fout1,dimNames1,dimSizes1,dimUnlim1)
                        
;  filevardef(fout1, "lat"  ,typeof(lat),getvardims(lat))                          
;  filevardef(fout1, "lon"  ,typeof(lon),getvardims(lon))                          
;  filevardef(fout1, "vor"  ,typeof(vor),getvardims(vor))    
;  filevardef(fout1, "u"   ,typeof(u) ,getvardims(u))    
;  filevardef(fout1, "v"   ,typeof(v) ,getvardims(v))        

;  filevarattdef(fout1,"vor", vor)  
;  filevarattdef(fout1,"u", u)    
;  filevarattdef(fout1,"v", v)                                                               
;  filevarattdef(fout1,"lat"  ,lat)                     
;  filevarattdef(fout1,"lon"  ,lon)                                

;  setfileoption(fout1,"DefineMode",False)
   
;  fout1->lat    = (/lat/)
;  fout1->lon    = (/lon/) 
;  fout1->vor    = (/vor/)
;  fout1->u      = (/u/)
;  fout1->v      = (/v/)

;  print("store file 1 done.")

;;fout2

;  dimNames2 = (/"time", "lat", "lon"/)  
;  dimSizes2 = (/ -1   ,  dimsizes(lat),  dimsizes(lon)/) 
;  dimUnlim2 = (/ True , False, False/)
;  filedimdef(fout2,dimNames2,dimSizes2,dimUnlim2)

;  filevardef(fout2, "time" ,typeof(time),getvardims(time))                           
;  filevardef(fout2, "lat"  ,typeof(lat),getvardims(lat))                          
;  filevardef(fout2, "lon"  ,typeof(lon),getvardims(lon))                          
;  filevardef(fout2, "vor"  ,typeof(xa),getvardims(xa))    
;  filevardef(fout2, "u"   ,typeof(uc) ,getvardims(uc))    
;  filevardef(fout2, "v"   ,typeof(vc) ,getvardims(vc))        

;  filevarattdef(fout2,"vor", xa)  
;  filevarattdef(fout2,"u", uc)    
;  filevarattdef(fout2,"v", vc)                          
;  filevarattdef(fout2,"time" ,time)                                      
;  filevarattdef(fout2,"lat"  ,lat)                     
;  filevarattdef(fout2,"lon"  ,lon)                                

;  setfileoption(fout2,"DefineMode",False)

;  fout2->time   = (/time/)     
;  fout2->lat    = (/lat/)
;  fout2->lon    = (/lon/) 
;  fout2->vor    = (/xa/)
;  fout2->u      = (/uc/)
;  fout2->v      = (/vc/)

;  print("store file 2 done.")

;calc u and v mean.
  fm = addfile(rmdir,"r")
  fi = addfile(ridir,"r")

  vom = fm->vor({lats:latn},{lonw:lone})
  um  = fm->u({lats:latn},{lonw:lone})
  vm  = fm->v({lats:latn},{lonw:lone})

;read in ui and calc rws.
  voi = fi->vor(:,{lats:latn},{lonw:lone})
  ui = fi->u(:,{lats:latn},{lonw:lone})
  vi = fi->v(:,{lats:latn},{lonw:lone})

;calc RWS.
  rws = new((/dimsizes(ui&time),dimsizes(ui&lat),dimsizes(ui&lon),5/),float)

  do i = 0, dimsizes(ui&time) - 1

    if ( any(ismissing(ui(i,:,:))) ) then
      continue
    else
      rws(i,:,:,:) = calc_rws(ui(i,:,:),um,vi(i,:,:),vm,voi(i,:,:),vom)
    end if

  end do

  copy_VarCoords(ui, rws(:,:,:,0))
  rws@unit = "s^(-2)"
  rws@var = "0: general, 1: divi&vorm 2: i_mvor adv 3: divm&vori 4: m_ivor adv"

  rws!3 = "terms"
  terms = ispan(0,4,1)
  rws@terms = terms

  lat  = rws&lat
  lon  = rws&lon
  time = rws&time

  system("/bin/rm -f " + ddir+"/rws.nc") 

  setfileoption("nc","Format","NetCDF4Classic")
  setfileoption("nc","CompressionLevel",4)
  fout = addfile(ddir+"/rws.nc", "c")

  setfileoption(fout,"DefineMode",True)

  dimNames = (/"time", "lat", "lon", "terms"/)  
  dimSizes = (/ -1   ,  dimsizes(lat),  dimsizes(lon), dimsizes(terms)/) 
  dimUnlim = (/True, False, False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout, "terms" ,typeof(terms),getvardims(terms)) 
  filevardef(fout, "time" ,typeof(time),getvardims(time))                           
  filevardef(fout, "lat"  ,typeof(lat),getvardims(lat))                          
  filevardef(fout, "lon"  ,typeof(lon),getvardims(lon))                          
  filevardef(fout, "rws"  ,typeof(rws),getvardims(rws))           

  filevarattdef(fout,"rws", rws)                           
  filevarattdef(fout,"time" ,time)                                      
  filevarattdef(fout,"lat"  ,lat)                     
  filevarattdef(fout,"lon"  ,lon)                              

  setfileoption(fout,"DefineMode",False)

  fout->terms  = (/terms/) 
  fout->time   = (/time/)     
  fout->lat    = (/lat/)
  fout->lon    = (/lon/) 
  fout->rws    = (/rws/)

  print("store file done.")