#!/bin/bash
# InFo
# select v250 from 1979 to 2022.
# author: cui xiangyang. time: 2023-11-21

  echo "prepare v700 begin. $(date)"

  export Ddir="/home/sunming/data5/cuixy/DATA/ERA-5/daily/vprs/raw_data"
  export odir="/home/sunming/data5/cuixy/global_jets/data"

  for date in $(seq -w 01 23);do
    cdo sellevel,700  ${Ddir}/day_v_flt_${date}.nc ${odir}/tmp_day_v_flt_${date}.nc
  done

  cdo mergetime ${odir}/tmp_day_v_flt_??.nc ${odir}/tmp_ERA5_daily_v700.nc

  cdo seldate,1979-01-01,2022-12-31  ${odir}/tmp_ERA5_daily_v700.nc ${odir}/ERA5_daily_v700_1979-2022_all.nc

  echo "store date done.  $(date)"

  cdo remapbil,r144x73 ${odir}/ERA5_daily_v700_1979-2022_all.nc ${odir}/ERA5_daily_v700_1979-2022_r2.5.nc

  rm -f ${odir}/tmp*.nc
  rm -f ${odir}/ERA5_daily_v700_1979-2022_all.nc

  ncl -Q 'var="v"' ./calc_filtered_values.ncl

  echo "calc filtered values done. $(date)"

### for u700.

  echo "prepare u700 begin. $(date)"

  export Ddir="/home/sunming/data5/cuixy/DATA/ERA-5/daily/uprs/raw_data"
  export odir="/home/sunming/data5/cuixy/global_jets/data"

  for date in $(seq -w 01 23);do
    cdo sellevel,700  ${Ddir}/day_u_flt_${date}.nc ${odir}/tmp_day_u_flt_${date}.nc
  done

  cdo mergetime ${odir}/tmp_day_u_flt_??.nc ${odir}/tmp_ERA5_daily_u700.nc

  cdo seldate,1979-01-01,2022-12-31  ${odir}/tmp_ERA5_daily_u700.nc ${odir}/ERA5_daily_u700_1979-2022_all.nc

  echo "store date done. $(date)"

  cdo remapbil,r144x73 ${odir}/ERA5_daily_u700_1979-2022_all.nc ${odir}/ERA5_daily_u700_1979-2022_r2.5.nc

  rm -f ${odir}/tmp*.nc
  rm -f ${odir}/ERA5_daily_u700_1979-2022_all.nc

  ncl -Q 'var="u"' ./calc_filtered_values.ncl

  echo "calc filtered values done. $(date)"
