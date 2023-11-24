#! /bin/bash
#Information
#Usage: This script is used to ks 250hPa uwind from uwind file and merge to a single file.
#times:1979-2022
#author: cui xiangyang 
#time: 2023-10-21
#release: v1
#-----------------------------------------------------------------------
DATADIR="/home/sunming/data5/cuixy/DATA/ERA-5"
VARDIR="daily/uprs/raw_data"
VARPRE="day_u_flt_"
VAREND="nc"

OUTDIR="/home/sunming/data5/cuixy/data/MJO_pan-Atlantic"
OUTNAM="u250.daily.1x1.nc"
#-----------------------------------------------------------------------
#Judge for the existence of directory.

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi
#-----------------------------------------------------------------------
#Main code.

num=1

while [ $num -le 23 ]

do

    inum=$( printf "%02d" ${num} )
    cdo -sellonlatbox,-80,55,5,85 -sellevel,250 ${DATADIR}/${VARDIR}/${VARPRE}${inum}.nc ${OUTDIR}/tmp.${inum}.nc

    let num=num+1
done

#-----------------------------------------------------------------------
#concactenente files into a single file.

ncrcat -O ${OUTDIR}/tmp.*.nc ${OUTDIR}/${OUTNAM}
ncrename -O -d latitude,lat -d longitude, lon -v latitude,lat -v longitude,lon ${OUTDIR}/${OUTNAM} ${OUTDIR}/${OUTNAM}


rm -f ${OUTDIR}/tmp.*.nc
 

