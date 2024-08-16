#! /bin/bash

export ddir="/home/sunming/data5/cuixy/DATA/ERA-5/daily/tprs/post_data/"
export odir="/home/sunming/data5/cuixy/global_jets/data/NAJ"

  cdo -O sellevel,850 \
    -sellonlatbox,-180,180,0,90 \
  ${ddir}/ERA5_t_1979-2022_regrid.nc ${odir}/ERA5_daily_mean_t850_1979-2022_r2.5.nc
