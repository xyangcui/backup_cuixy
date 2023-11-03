#! /bin/bash

export datadir="/home/sunming/data5/cuixy/data/MJO_pan-Atlantic"
export fnamebe="ERA5"
export var="${1}"
export var2="10m_v_component_of_wind" 
declare -a flist
declare -a flist2
#short2flt and mergetime.

export yr=1979
while [ $yr -le 2022 ] 
do
  num=$[ $yr - 1978 ]
  flist[num]="${datadir}/ERA5_${var}/${fnamebe}_${var}_hourly_${yr}.nc"

  let yr=yr+1
done

cdo -O -b 32 -mergetime ${flist[@]}  ${datadir}/${fnamebe}_${var}/tmp1_${fnamebe}_${var}_hourly_1979-2022.nc

#daily mean values and retrive 11.2019-----3.2022.
if [ "${var}" = "2m_temperature" ];then

  cdo -O -daymean -seldate,1979-11-01,2022-03-31  ${datadir}/${fnamebe}_${var}/tmp1_${fnamebe}_${var}_hourly_1979-2022.nc \
${datadir}/${fnamebe}_${var}/tmp_${fnamebe}_${var}_1979-2022.nc

  ncrename -O -d longitude,lon -v longitude,lon -d latitude,lat -v latitude,lat ${datadir}/${fnamebe}_${var}/tmp_${fnamebe}_${var}_1979-2022.nc \
${datadir}/ERA5_${var}/${fnamebe}_${var}_1979-2022.nc

elif [ "${var}" = "total_precipitation" ];then

  cdo -O -daysum -shifttime,-1hour -seldate,1979-11-01,2022-03-31 -selmon,11,12,1,2,3 ${datadir}/${fnamebe}_${var}/tmp1_${fnamebe}_${var}_hourly_1979-2022.nc \
${datadir}/${fnamebe}_${var}/tmp_${fnamebe}_${var}_1979-2022.nc

  ncrename -O -d longitude,lon -v longitude,lon -d latitude,lat -v latitude,lat ${datadir}/${fnamebe}_${var}/tmp_${fnamebe}_${var}_1979-2022.nc \
${datadir}/ERA5_${var}/${fnamebe}_${var}_1979-2022.nc

else
# calculate wind speed and retrive maximum wind speed.
  
  for yr in {1979..2022}; do
    num=$[ $yr - 1978 ]
    flist2[num]="${datadir}/ERA5_${var2}/${fnamebe}_${var2}_hourly_${yr}.nc"
  done
  
  cdo -O -b 32 -mergetime ${flist2[@]} ${datadir}/${fnamebe}_${var2}/tmp1_${fnamebe}_${var2}_hourly_1979-2022.nc

  cdo merge  ${datadir}/${fnamebe}_${var}/tmp1_${fnamebe}_${var}_hourly_1979-2022.nc ${datadir}/${fnamebe}_${var2}/tmp1_${fnamebe}_${var2}_hourly_1979-2022.nc \
  ${datadir}/${fnamebe}_${var}/tmp1_${fnamebe}_s10_hourly_1979-2022.nc

  cdo expr,'s10=sqrt(u10*u10+v10*v10)' ${datadir}/${fnamebe}_${var}/tmp1_${fnamebe}_s10_hourly_1979-2022.nc \
  ${datadir}/${fnamebe}_${var}/tmp2_${fnamebe}_s10_hourly_1979-2022.nc

  cdo daymean -seldate,1979-11-01,2022-03-31 ${datadir}/${fnamebe}_${var}/tmp2_${fnamebe}_s10_hourly_1979-2022.nc ${datadir}/${fnamebe}_${var}/tmp_${fnamebe}_s10_1979-2022.nc 

  ncrename -O -d longitude,lon -v longitude,lon -d latitude,lat -v latitude,lat ${datadir}/${fnamebe}_${var}/tmp_${fnamebe}_s10_1979-2022.nc \
${datadir}/ERA5_${var}/${fnamebe}_s10_1979-2022.nc
  
  rm -f ${datadir}/${fnamebe}_${var2}/tmp1*.nc

fi

rm -f ${datadir}/${fnamebe}_${var}/tmp*.nc
unset datadir;unset fnamebe;unset var;unset var2

