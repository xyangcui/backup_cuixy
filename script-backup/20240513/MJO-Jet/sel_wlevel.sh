#! /bin/bash

export ddir="/home/sunming/data5/cuixy/DATA/ERA-5/daily/wprs/post_data"
export odir="/home/sunming/data5/cuixy/global_jets/data/NAJS_temp_budget"

cdo -O sellevel,250,300,350,400,450,500 -sellonlatbox,180,360,0,90 \
 ${ddir}/ERA5_w_1979-2022_regrid.nc ${odir}/ERA5_daily_mean_wlevel_1979-2022_r2.5.nc

export ddir="/home/sunming/data5/cuixy/DATA/ERA-5/daily/tprs/post_data"
cdo -O sellevel,200,250,300,350,400,450,500,550 -sellonlatbox,180,360,0,90 \
 ${ddir}/ERA5_t_1979-2022_regrid.nc ${odir}/ERA5_daily_mean_tlevel_1979-2022_r2.5.nc

export ddir="/home/sunming/data5/cuixy/DATA/ERA-5/daily/uprs/post_data"
cdo -O sellevel,250,300,350,400,450,500 -sellonlatbox,180,360,0,90 \
 ${ddir}/ERA5_u_1979-2022_regrid.nc ${odir}/ERA5_daily_mean_ulevel_1979-2022_r2.5.nc

export ddir="/home/sunming/data5/cuixy/DATA/ERA-5/daily/vprs/post_data"
cdo -O sellevel,250,300,350,400,450,500 -sellonlatbox,180,360,0,90 \
 ${ddir}/ERA5_v_1979-2022_regrid.nc ${odir}/ERA5_daily_mean_vlevel_1979-2022_r2.5.nc