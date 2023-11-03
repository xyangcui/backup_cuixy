#! /bin/bash
#---------------------------------------------------------
# This script is execute shell script. calculating u250 anomalies
# and plotting.
# Information
# author: cui xiangyang 
# time: 2023-10-3
#---------------------------------------------------------
# calculate
export DATADIR="/home/sunming/data5/cuixy/data/MJO_pan-Atlantic/u250.daily.1x1.nc"
export OUTDIR="/home/sunming/data5/cuixy/data/MJO_pan-Atlantic/u250.anom.1x1.nc"
export var="u"
export SCRIPTDIR="/home/sunming/data5/cuixy/script/MJO_pan-Atlantic/ncl_script"
export outdir="/home/sunming/data5/cuixy/data/MJO_pan-Atlantic/u250.anom.NDJFMA.1x1.nc"
export latmin=15
export latmax=75
export lonmin=-80
export lonmax=0


if [ ! -f $OUTDIR ]; then	
    ncl ${SCRIPTDIR}/u250anom.ncl
fi

if [ ! -f $outdir ]; then
    cdo -O -selmon,1,2,3,4,11,12 -seldate,1979-11-01,2021-04-30 ${OUTDIR} ${outdir} 
fi

#---------------------------------------------------------
# plot
export rmmdir="/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
export plotdir="/home/sunming/data5/cuixy/plot/"
export i=1



    export plotdir="${plotdir}/lagu250_${i}"
    ncl ${SCRIPTDIR}/lagplot.ncl
    export plotdir="/home/sunming/data5/cuixy/plot/"


unset DATADIR;unset outdir;unset var;unset SCRIPTDIR;unset OUTDIR
unset plotdir;unset rmmdir
#--------------------------------------------------------

