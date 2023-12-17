#! /bin/bash

export ddir="/home/sunming/data5/cuixy/DATA/ERA-5/daily/tprs/post_data"
export odir="/home/sunming/data5/cuixy/global_jets/data"

#for date in $( seq -w 21 );do
#     cdo -O sellevel,200,250,300,350,400,450,\
#500,600,700,850,925,1000 -sellonlatbox,-180,180,-30,90 \
#${ddir}/day_t_flt_${date}.nc ${odir}/day_t_flt_${date}.nc

#done

     cdo -O sellevel,200,250,300,350,400,450,\
500,600,700,850,925,1000 -sellonlatbox,-180,180,-30,90 \
${ddir}/ERA5_t_1979-2022_regrid.nc ${odir}/ERA5_daily_mean_tlevel_1979-2022.nc


