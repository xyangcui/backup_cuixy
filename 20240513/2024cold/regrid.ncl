load "/home/sunming/data5/cuixy/script/ncl_func/create_nc_file_3D.ncl"
;--------------------------------------------------
; 1. regrid from gauss grid to fixed grid(73x144).
; 2. final dimensions (19800101-20240228,73,144).
; time: 2024-03-29  author: cui xiangyang 
;--------------------------------------------------

ddir = "/home/sunming/data5/cuixy/2024cold/data/air.2m.gauss.1980-2024.nc"
tdir = "/home/sunming/data5/cuixy/2024cold/data/ncep_t2m_mean_daily_1980-2024.nc"
;------------------------------------------------------------
; Main code.
;------------------------------------------------------------

  f = addfile(ddir,"r")

  t = f->air

  nt = g2fsh_Wrap(t, (/73,144/))

  printVarSummary(nt)

  create_nc_file(nt,nt&time,nt&lat,nt&lon,tdir,"t")