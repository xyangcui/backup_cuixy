;-------------------------------------------------------
;acc skill.
;-------------------------------------------------------\
procedure calc_acc(event,mdir,odir,lats,latn,lonw,lone,mname)
local mdir,odir,lats,latn,lonw,lone,mname,event
begin
  ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"
  print("data preparation begin.")
  ;read in model value.
  f = addfile(mdir,"r")
  u = f->u(:,:,{lats:latn},{lonw:lone})
  valid_time = f->valid_time 
  TIMI = cd_calendar(valid_time,-2)
  ;read in obs value.
  fo = addfile(odir,"r")
  steps = u&lead_time
  if steps(0).eq.24 then
    ut = fo->u(:,1:,{lats:latn},{lonw:lone})
  else
    ut = fo->u(:,:,{lats:latn},{lonw:lone})
  end if
  dims = dimsizes(u)
  print("data preparation done.")

  uECW = new(dimsizes(u(0,:,:,:)),float)
  uECW = 1.
  do step = 0,dims(1)-1
    tmp = TIMI(:,step) ;step's dates
    ;get indices.
    Da  = new(dimsizes(tmp),integer)
    do i = 0, dimsizes(tmp)-1
      Dt = ind(event.eq.tmp(i))
      if (.not.ismissing(Dt)) then
        Da(i) = i
      end if
    end do
    igood = ind(.not.ismissing(Da))
    Da := Da(igood) ;here is the steps indice.
    print("step "+step+": "+dimsizes(Da))
    um = u(Da,step,:,:)
    uo = ut(Da,step,:,:)

    ;ACC skill.
    um_m = dim_avg_n_Wrap(um, 0)
    um := um - conform_dims(dimsizes(um), um_m, (/1,2/))
    uo_m = dim_avg_n_Wrap(uo, 0)
    uo := uo - conform_dims(dimsizes(uo), uo_m, (/1,2/))
    acc1 = dim_sum_n(um*uo,0)
    acc2 = sqrt(dim_sum_n(uo*uo,0))
    acc3 = sqrt(dim_sum_n(um*um,0))
    uECW(step,:,:) = doubletofloat(acc1/(acc2*acc3))

    delete([/um,uo,igood,Da,tmp/])
  end do

  copy_VarCoords(u(0,:,:,:), uECW)
  system("/bin/rm -f " + ddir+"/"+"acc_"+model+".nc")
  fout  = addfile(ddir+"/"+"acc_"+model+".nc", "c")
  fout->u = uECW
  print("acc "+mname+" end.")

end


lats = 10
latn = 80
lonw = 100
lone = 240
ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"
;take jet events. as day0.
rdir = ddir+"/NPJ_index.csv"
delim  = ","
rf   = asciiread(rdir,-1,"string")
NPJ  = str_get_field(rf,2,delim)
dat  = str_get_field(rf,1,delim)
NPJ := tofloat(NPJ(1:))
dat := todouble(dat(1:))

event = dat(ind(NPJ.gt.1.5.or.NPJ.lt.-1.5)) ;NPJ events' days.


if model.ne."LIM" then

  mdir = ddir+"/"+model+"/"+model+"_pf_anom_u200.nc"
  odir = ddir+"/"+model+"/"+model+"_obs_anom_u200_v2.nc"
  calc_acc(event,mdir,odir,lats,latn,lonw,lone,model)

else

  mdir = ddir+"/LIM_hindcast_u200_v2.nc"
  odir = ddir+"/LIM_obs_u200_v2.nc"
  ;read in model value.
  f = addfile(mdir,"r")
  u = f->u(:,:,{lats:latn},{lonw:lone})
  TIME = cd_calendar(u&time,-2)
  ;read in obs value.
  fo = addfile(odir,"r")
  ut = fo->u(:,:,{lats:latn},{lonw:lone})

  ;give ut a valid time.
  odir = ddir+"/ERA5_daily_anom_u200_1979-2024.nc"
  ft = addfile(odir,"r")
  timet= cd_calendar(ft->time,-2)
  TIMI = new((/dimsizes(TIME),47/),double)
  do i = 0, dimsizes(TIME)-1
    nt = ind(timet.eq.TIME(i))
    TIMI(i,:) = timet(nt:nt+46)
    delete(nt)
  end do
  uECW = new(dimsizes(u(0,:,:,:)),float)
  uECW = 1.

  do step = 0,46
    tmp = TIMI(:,step) ;step's dates
    ;get indices.
    Da  = new(dimsizes(tmp),integer)
    do i = 0, dimsizes(tmp)-1
      Dt = ind(event.eq.tmp(i))
      if (.not.ismissing(Dt)) then
        Da(i) = i
      end if
    end do
    igood = ind(.not.ismissing(Da))
    Da := Da(igood) ;here is the steps indice.
    print("step "+step+": "+dimsizes(Da))
    um = u(Da,step,:,:)
    uo = ut(Da,step,:,:)

    ;ACC skill.
    if dimsizes(Da).eq.1 then
      uECW(step,:,:) = 0.
    else
      um_m = dim_avg_n_Wrap(um, 0)
      um := um - conform_dims(dimsizes(um), um_m, (/1,2/))
      uo_m = dim_avg_n_Wrap(uo, 0)
      uo := uo - conform_dims(dimsizes(uo), uo_m, (/1,2/))
      acc1 = dim_sum_n(um*uo,0)
      acc2 = sqrt(dim_sum_n(uo*uo,0))
      acc3 = sqrt(dim_sum_n(um*um,0))
      uECW(step,:,:) = doubletofloat(acc1/(acc2*acc3))
    end if

    delete([/um,uo,igood,Da,tmp/])
  end do

  copy_VarCoords(u(0,:,:,:), uECW)
  system("/bin/rm -f " + ddir+"/acc_LIM.nc")
  fout  = addfile(ddir+"/acc_LIM.nc", "c")
  fout->u = uECW
  print("acc LIM end.")

end if