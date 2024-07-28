#! /bin/bash
#---------------------------------------------------------
# This script is execute shell script. calculating jet indices 
# and plotting.
# Information
# author: cui xiangyang 
# time: 2023-10-22
#---------------------------------------------------------
# calculate
DATADIR1="/home/sunming/data5/cuixy/data/MJO_pan-Atlantic/u250.daily.1x1.nc"
outdir="/home/sunming/data5/cuixy/data/MJO_pan-Atlantic/jetindices.nc"
var="u"
SCRIPTDIR="/home/sunming/data5/cuixy/script/MJO_pan-Atlantic/ncl_script"
DATADIR="/home/sunming/data5/cuixy/data/MJO_pan-Atlantic/u250.NDJFMA.1x1.nc"

export DATADIR;export outdir;export var;export SCRIPTDIR
if [ ! -f $outdir ]; then
    cdo -selmon,1,2,3,4,11,12 -seldate,1979-11-01,2021-04-30 ${DATADIR1} ${DATADIR} 
    ncl  ${SCRIPTDIR}/Jet_indices.ncl
fi

#---------------------------------------------------------
# plot
var1="jets"
var2="jetl"
var3="jetz"
rmmdir="/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
plotdir="/home/sunming/data5/cuixy/plot/jetindices"

export var1;export var2;export var3; export rmmdir;export plotdir

ncl  ${SCRIPTDIR}/plot.ncl

unset DATADIR;unset outdir;unset var;unset SCRIPTDIR
unset plotdir;unset var1;unset var2;unset var3;unset rmmdir
#--------------------------------------------------------
