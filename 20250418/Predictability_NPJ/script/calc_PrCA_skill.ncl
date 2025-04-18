
  ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data/"

  ;observation
  f = addfile(ddir+"/PrCA_obs_variates.nc","r")
  vo = f->v  ;(mode,l_time,time)

  ;models
  f = addfile(ddir+"/PrCA_variates.nc","r")
  v = dim_avg_n(f->v,1) ;(l_time,time,mode)
  dims = dimsizes(v)

  ac = new((/dims(2),dims(0)/),float)
  printVarSummary(ac)

  do i = 0, dims(2)-1  ;mode
    do j = 0, dims(0)-1  ;l_time
      vm_m = dim_avg(vo(i,j,:))
      vom  = vo(i,j,:) - vm_m

      v_m = dim_avg(v(j,:,i))
      vm  = v(j,:,i) - v_m

      acc1 = dim_sum_n(vom*vm,0)
      acc2 = sqrt(dim_sum_n(vom*vom,0))
      acc3 = sqrt(dim_sum_n(vm*vm,0))
      ac(i,j)= doubletofloat(acc1/(acc2*acc3))
    end do
  end do

  ac!0 = "mode"
  ac&mode = ispan(1,20,1)
  ac!1 = "step"
  ac&step = ispan(0,46,1)

  system("/bin/rm -f " + ddir+"/"+"PrCA_ac.nc")
  fout  = addfile(ddir+"/"+"PrCA_ac.nc", "c")
  fout->ac = ac