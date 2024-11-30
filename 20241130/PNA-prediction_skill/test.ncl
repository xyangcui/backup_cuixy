  ddir = "/home/sunming/data5/cuixy/MJO_PC_diversiy/data"

  odir = ddir+"/"+"boostrap_pc_omi_alldays.nc"

  fout = addfile(odir,"r")           
  pc   = fout->pc
  print(max(pc))
  m_value2 = dim_pqsort_n(pc, 2, 1) 

  std1 = dim_stddev_n(pc, 1)
  print(std1)
  mean1 = dim_avg_n(pc, 1)
  print(mean1)
  sigh = 0.99

  nl2  = toint(dimsizes(pc(1,:,7))*sigh) - 1
  nh2  = toint(dimsizes(pc(0,:,7))*1) - 1
  
  pcf  = pc(0,nl2,7)
  pcs  = pc(1,nl2,7)

  print(pcs)
  print(pcf)