#!/bin/bash
# InFo
# select u250 from 1979 to 2022.
# author: cui xiangyang. time: 2023-12-12
export Ddir="/home/sunming/data5/cuixy/DATA/ERA-5/daily/zprs"
export odir="/home/sunming/data5/cuixy/global_jets/data"

for date in $(seq -w 01 20);do
    cdo -sellevel,250 ${Ddir}/raw_data/day_z_flt_${date}.nc ${odir}/tmp_day_z_flt_${date}.nc
done
cdo -b 32 copy  ${Ddir}/raw_data_2020/z_202001-202006.nc ${odir}/tmp1_z_202001-202006.nc
cdo -b 32 copy  ${Ddir}/raw_data_2020/z_202007-202012.nc ${odir}/tmp1_z_202007-202012.nc

cdo dayavg  ${odir}/tmp1_z_202001-202006.nc ${odir}/tmp_z_202001-202006.nc
cdo dayavg  ${odir}/tmp1_z_202007-202012.nc ${odir}/tmp_z_202007-202012.nc

cdo -sellevel,250 -mergetime ${odir}/tmp_z_*.nc ${odir}/tmp_day_z_flt_21.nc

cdo -O -mergetime ${odir}/tmp_day_z_flt_??.nc ${odir}/ERA5_daily_z250_1979-2020.nc

rm -f ${odir}/tmp*.nc
