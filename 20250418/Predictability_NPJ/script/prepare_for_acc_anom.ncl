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
  dimstep = dimsizes(step)
  fstep = step(dimstep-1)
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
    ;um&lead_time = ispan(0,fstep+24,24)
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
odir= ddir+"/ecmwf/ERA5_daily_anom_u200_1979-2024.nc"

;-----------------------------------------------
;Main code.
;-----------------------------------------------
  fo = addfile(odir,"r")
  uo = fo->u(:,{lats:latn},{lonw:lone})
  TIME2 = cd_calendar(uo&time,-2)

;;obs for kma
  obs_for_model(uo,TIME2,mdir,lats,latn,lonw,lone,\
                model,ddir+"/"+model+"/"+model+"_obs_anom_u200_v2.nc")