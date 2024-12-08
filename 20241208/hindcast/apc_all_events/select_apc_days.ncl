;-----------------------------------------------
;read in csv file.
function return_pday(Dir)
begin
  delim = ","
  f = asciiread(Dir,-1,"string")
  p = str_get_field(f,1,delim)
  p := toint(p(1::))

  return p
end

procedure store_model_values_cf(outdir,time,lead,lat,lon,var,varname)
local outdir,time,lead,lat,lon,var,varname
begin
  system("/bin/rm -f "+outdir)

  setfileoption("nc", "Format", "NetCDF4Classic")
  setfileoption("nc","CompressionLevel",5)

  fout = addfile(outdir,"c")
  setfileoption(fout,"DefineMode",True)
  
  fAtt               = True           
  fAtt@creation_date = systemfunc ("date")
  fileattdef( fout, fAtt )

  dimNames = (/"initial_time", "lead_time", "lat", "lon"/)  
  dimSizes = (/ -1   ,  dimsizes(lead), dimsizes(lat),  dimsizes(lon)/) 
  dimUnlim = (/ True , False, False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout, "initial_time" ,typeof(time),getvardims(time)) 
  filevardef(fout, "lead_time" ,typeof(lead),getvardims(lead))                          
  filevardef(fout, "lat"  ,typeof(lat),getvardims(lat))                          
  filevardef(fout, "lon"  ,typeof(lon),getvardims(lon))                          
  filevardef(fout, varname  ,typeof(var)  ,getvardims(var))    
                                                                                             
  filevarattdef(fout,"initial_time" ,time)  
  filevarattdef(fout,"lead_time" ,lead)                   
  filevarattdef(fout,"lat"  ,lat)                   
  filevarattdef(fout,"lon"  ,lon)                     
  filevarattdef(fout,varname ,var)   

  setfileoption(fout,"DefineMode",False)

  fout->initial_time  = (/time/) 
  fout->lead_time     = (/lead/)     
  fout->lat           = (/lat/)
  fout->lon           = (/lon/) 
  fout->$varname$     = (/var/)
end

procedure store_model_values_pf(outdir,ens,time,lead,lat,lon,var,varname)
local outdir,time,lead,lat,lon,var,varname
begin
  system("/bin/rm -f "+outdir)

  setfileoption("nc", "Format", "NetCDF4Classic")
  setfileoption("nc","CompressionLevel",5)

  fout = addfile(outdir,"c")
  setfileoption(fout,"DefineMode",True)
  
  fAtt               = True           
  fAtt@creation_date = systemfunc ("date")
  fileattdef( fout, fAtt )

  dimNames = (/"number", "initial_time", "lead_time", "lat", "lon"/)  
  dimSizes = (/dimsizes(ens), -1   ,  dimsizes(lead), dimsizes(lat),  dimsizes(lon)/) 
  dimUnlim = (/False, True , False, False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout, "number" ,typeof(ens),getvardims(ens))
  filevardef(fout, "initial_time" ,typeof(time),getvardims(time)) 
  filevardef(fout, "lead_time" ,typeof(lead),getvardims(lead))
  filevardef(fout, "lat"  ,typeof(lat),getvardims(lat))                          
  filevardef(fout, "lon"  ,typeof(lon),getvardims(lon))                          
  filevardef(fout, varname  ,typeof(var)  ,getvardims(var))    

  filevarattdef(fout,"number" ,ens)                                                                                           
  filevarattdef(fout,"initial_time" ,time)  
  filevarattdef(fout,"lead_time" ,lead)                   
  filevarattdef(fout,"lat"  ,lat)                   
  filevarattdef(fout,"lon"  ,lon)                     
  filevarattdef(fout,varname ,var)   

  setfileoption(fout,"DefineMode",False)

  fout->number        = (/ens/)
  fout->initial_time  = (/time/) 
  fout->lead_time     = (/lead/)     
  fout->lat           = (/lat/)
  fout->lon           = (/lon/) 
  fout->$varname$     = (/var/)
end

procedure store_model_values_ob(outdir,time,lead,lat,lon,var,varname)
local outdir,time,lead,lat,lon,var,varname
begin
  system("/bin/rm -f "+outdir)

  setfileoption("nc", "Format", "NetCDF4Classic")
  setfileoption("nc","CompressionLevel",5)

  fout = addfile(outdir,"c")
  setfileoption(fout,"DefineMode",True)
  
  fAtt               = True           
  fAtt@creation_date = systemfunc ("date")
  fileattdef( fout, fAtt )

  dimNames = (/"time", "lead_time", "lat", "lon"/)  
  dimSizes = (/ -1   , dimsizes(lead), dimsizes(lat),  dimsizes(lon)/) 
  dimUnlim = (/ True , False, False, False, False/)   
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout, "time" ,typeof(time),getvardims(time)) 
  filevardef(fout, "lead_time" ,typeof(lead),getvardims(lead))                           
  filevardef(fout, "lat"  ,typeof(lat),getvardims(lat))                          
  filevardef(fout, "lon"  ,typeof(lon),getvardims(lon))                          
  filevardef(fout, varname  ,typeof(var)  ,getvardims(var))    
                                                                                             
  filevarattdef(fout,"time" ,time)  
  filevarattdef(fout,"lead_time" ,lead)                   
  filevarattdef(fout,"lat"  ,lat)                   
  filevarattdef(fout,"lon"  ,lon)                     
  filevarattdef(fout,varname ,var)   

  setfileoption(fout,"DefineMode",False)

  fout->time       = (/time/)  
  fout->lead_time  = (/lead/)    
  fout->lat        = (/lat/)
  fout->lon        = (/lon/) 
  fout->$varname$  = (/var/)
end

model_name = "ecmwf"
ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data/"

cfdir = ddir+"hindcast/ecmwf/"+Year+"/data/"+model_name+"_cf_z500_anom_r2.5.nc"
pfdir = ddir+"hindcast/ecmwf/"+Year+"/data/"+model_name+"_pf_z500_anom_r2.5.nc"
obdir = ddir+"hindcast/ecmwf/"+model_name+"_z500_daily_anom_1979-2022.nc"
;obdir = ddir+"hindcast/ecmwf/"+Year+"/data/"+model_name+"_z500_daily_anom_1979-2022.nc"

sdir = new((/8/),string)
do i = 1,8
  sdir(i-1) = ddir+"pdate/"+"IO_p"+i+"_a1.csv"
end do
;--------------------------------------------------------------
;main code.
  f1 = addfile(cfdir,"r") ;control forecast
  f2 = addfile(pfdir,"r") ;perturbed forecast
  f3 = addfile(obdir,"r") ;observation

  xc = f1->gh     ;unit:gpm
  xp = f2->gh     ;unit:gpm
  xot = f3->z     ;unit:m2/s2
  xo = xot/9.8    ;unit:gpm
  copy_VarCoords(xot,xo)
  copy_VarAtts(xot,xo)

  number = xp&number
  number2= toint(number)
  copy_VarAtts(number, number2)  

  model_time = xc&initial_time
  model_time2= todouble(model_time)
  copy_VarAtts(model_time, model_time2)
  printVarSummary(model_time2)
  model_TIME = cd_calendar(model_time2, -2)
  obser_time = xo&time
  obser_TIME = cd_calendar(obser_time, -2)
  delete([/model_time,xc&initial_time,xp&initial_time,xp&number/])
  xc&initial_time = model_time2
  xp&initial_time = model_time2
  xp&number = number2

;read in csv file.
  do i = 0,7
    print("phase"+i+"begin. "+systemfunc("date +%m-%d"))
    sp = return_pday(sdir(i))
    mt_sp = new(dimsizes(sp),integer)
    ot_sp = mt_sp
    ;return related values.
    do t = 0,dimsizes(sp)-1
      nt = ind(model_TIME.eq.sp(t))
      mt_sp(t) = nt
    end do
    igood = ind(.not.ismissing(mt_sp))
    mt_sp := mt_sp(igood)
    delete(igood)

    model_sp = model_TIME(mt_sp)
    do t= 0,dimsizes(model_sp)-1
      nt = ind(obser_TIME.eq.model_sp(t))
      ot_sp(t) = nt
    end do
    igood = ind(.not.ismissing(ot_sp))
    ot_sp := ot_sp(igood)
    delete(igood)

    ph = i+1
    ;IOsP model
    xc_sp = xc(mt_sp,:,:,:)
    xp_sp = xp(:,mt_sp,:,:,:)

    x_sp = new((/11,dimsizes(xc_sp&initial_time),dimsizes(xc_sp&lead_time),dimsizes(xc_sp&lat),dimsizes(xc_sp&lon)/),float)
    x_sp = 0.
    copy_VarCoords(xc_sp,x_sp(0,:,:,:,:))
    x_sp!0 = "number"
    x_sp&number = ispan(0, 10, 1)
    x_sp(0,:,:,:,:) = xc_sp
    x_sp(1:,:,:,:,:) = xp_sp
    store_model_values_pf(ddir+"hindcast/"+model_name+"/"+Year+"/"+"IO_p"+ph+"_a1.nc",x_sp&number,x_sp&initial_time,\
                          x_sp&lead_time,x_sp&lat,x_sp&lon,x_sp,"gh")

    ;IOsP obser
    xo_sp = xc_sp
    do j = 1,dimsizes(xc&lead_time)
      xo_sp(:,j-1,:,:) = (/xo(ot_sp+j,:,:)/)
    end do
    copy_VarCoords(xc_sp, xo_sp)   
    store_model_values_cf(ddir+"hindcast/"+model_name+"/"+Year+"/"+"IO_ob_p"+ph+"_a1.nc",\
                          xo_sp&initial_time,xo_sp&lead_time,xo_sp&lat,xo_sp&lon,xo_sp,"z")

    delete([/xc_sp,xp_sp,xo_sp,mt_sp,ot_sp,sp,x_sp,model_sp/])
    print("phase"+i+"done. "+systemfunc("date +%m-%d"))

  end do