;------------------------------------------------
;prepare LIM hindcast u200.
;prepare obs for LIM and ecmwf.
;------------------------------------------------

procedure obs_for_model(obs,TIME2,indir,lats,latn,lonw,lone,mname,outdir)
local obs,TIME2,indir,lats,latn,lonw,lone,mname,outdir
begin
  print(mname +" obs begin")
  f = addfile(indir,"r")
  u = f->u(:,:,{lats:latn},{lonw:lone})
  dims = dimsizes(u)
  step = u&lead_time
  if step(0).eq.24 then
    um = new((/dims(0),dims(1)+1,dims(2),dims(3)/),float)
  else
    um = new(dimsizes(u),float)
  end if
  dims := dimsizes(um)
  timi = u&time
  time = todouble(timi)
  copy_VarAtts(timi, time)
  TIMI = cd_calendar(time,-2)

  do i = 0, dimsizes(TIMI)-1
    nt = ind(TIME2.eq.TIMI(i))
    um(i,:,:,:) = obs(nt:nt+dims(1)-1,:,:)
    delete(nt)
  end do
  u&time := time
  if step(0).eq.24 then
    copy_VarMeta(u(:,0,:,:), um(:,0,:,:))
    um!1 = "lead_time"
    um&lead_time = ispan(0,1104,24)
  else
    copy_VarMeta(u, um)
  end if
  printVarSummary(um)

  ;store data.
  system("/bin/rm -f " + outdir)
  setfileoption("nc", "Format", "NetCDF4Classic")
  setfileoption("nc","CompressionLevel",4)
  fout  = addfile(outdir, "c")
  fout->u = um

  print(mname +" obs end")
end 


lats = 10
latn = 80
lonw = 100
lone = 240


ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"
;mdir
mdir = ddir+"/"+model+"/"+model+"_pf_mean_u200.nc"

;LIM
ldir = ddir+"/LIM_hindcast.nc"

odir = ddir+"/ERA5_daily_anom_u200_1979-2024.nc"
o2dir= ddir+"/ecmwf/ERA5_daily_mean_u200_1979-2024.nc"

;-----------------------------------------------
;Main code.
;-----------------------------------------------
  fo = addfile(odir,"r")
  uo = fo->u(:,{lats:latn},{lonw:lone})
  TIME = cd_calendar(uo&time, -2)
  timi=cd_calendar(fo->time,0)
  mon = timi(:,1)
  nT = ind((TIME.ge.19821201.and.TIME.le.20220331)\
           .and.(mon.eq.12.or.mon.eq.1.or.mon.eq.2.or.mon.eq.3) )
;  uo := uo(nt,:,:)
  delete([/TIME,timi,mon/])

  TIME = cd_calendar(uo&time, -2)
  dims = dimsizes(uo)
  print(dimsizes(TIME))

  fo2 = addfile(o2dir,"r")
  uo2 = fo2->u(:,{lats:latn},{lonw:lone})
  TIME2 = cd_calendar(uo2&time,-2)

;;obs for kma
  ;obs_for_model(uo2,TIME2,mdir,lats,latn,lonw,lone,\
  ;              model,ddir+"/"+model+"/"+model+"_obs_mean_u200.nc")

;LIM hindcast u.
  print("LIM obs begin")
  n = 20-1
  ;get pcs.
  f = addfile(ldir,"r")
  x = f->x(0:n,:,:)
  ;get eofs.
  fl = addfile(ddir+"/eof_u200.nc","r")
  eof= fl->x(0:n,{lats:latn},{lonw:lone})
  eval= eof@eval(0:n)
  rad    = 4.*atan(1.)/180.
  clat   = uo&lat            
  clat   = sqrt( cos(rad*clat) ) 
  eofs= eof*sqrt(conform(eof, tofloat(eval), 0))
  copy_VarMeta(eof,eofs)
  u = new((/dimsizes(x(0,:,0)),47,dims(1),dims(2)/),float)
  do j = 0,dimsizes(x(0,:,0))-1
    u(j,:,:,:) = doubletofloat(eof2data_n(eofs,x(:,j,:),0))/conform(u(j,:,:,:), tofloat(clat), 1)
  end do

  copy_VarMeta(uo(nT,:,:),u(:,0,:,:))
  ;copy_VarCoords(uc, u(0,:,0,0))
  u!1 = "lead_time"
  u&lead_time = ispan(0,46*24,24)

  um = new(dimsizes(u),typeof(u))
  TIMI = cd_calendar(u&time,-2)
  do i = 0, dimsizes(TIMI)-1
    nt = ind(TIME.eq.TIMI(i))
    um(i,:,:,:) = uo(nt:nt+46,:,:)
    delete(nt)
  end do
  copy_VarMeta(u, um)
  um!1 = "lead_time"
  um&lead_time = ispan(0,46*24,24)
  printVarSummary(um)

  ;store data.
  system("/bin/rm -f " + ddir+"/LIM_hindcast_u200.nc")
  setfileoption("nc", "Format", "NetCDF4Classic")
  setfileoption("nc","CompressionLevel",4)
  fout  = addfile(ddir+"/LIM_hindcast_u200.nc", "c")
  fout->u = u

  system("/bin/rm -f " + ddir+"/LIM_obs_u200.nc")
  setfileoption("nc", "Format", "NetCDF4Classic")
  setfileoption("nc","CompressionLevel",4)
  fout  = addfile(ddir+"/LIM_obs_u200.nc", "c")
  fout->u = um

  delete([/f,u,TIMI,um/])
  print("LIM obs end")