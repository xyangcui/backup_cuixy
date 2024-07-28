#! /bin/bash

export ddir="/home/sunming/data5/cuixy/DATA/ERA-5/daily/uprs/raw_data"
export odir="/home/sunming/data5/cuixy/global_jets/data"

     cdo -O sellevel,250,300,350,400,450,\
500 -sellonlatbox,-180,180,-30,90 -mergetime\
 ${ddir}/day_u_flt_??.nc ${odir}/ERA5_daily_mean_ulevel_1979-2022.nc
