#! /bin/bash
###############################################################################
# This script calls for two ncl scripts.
# 1. jet_indices.ncl to calculate jet indices.
# 2. plot.ncl to plot
# $1: jet location; $2: latmin $3: latmax $4: lonmin $5: lonmax $6: box.
# Information
# author: cui xiangyang time: 2023-11-07
###############################################################################

export wkdir="/home/sunming/data5/cuixy/global_jets"
export indir="${wkdir}/data/ERA5_daily_u250_1979-2022_all.nc"
export udir="${wkdir}/data/ERA5_daily_u250_1979-2022.nc"
export outdir="${wkdir}/data/jet_indices_${1}_${6}.nc"
export plotdir="${wkdir}/plot/${1}_${6}"
export rmmdir="/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"

export var="u"
export var1="jets"
export var2="jetl"

export clat="-46"
export clon="6"

export latmin="$[ ${clat} - ${2} ]"
export latmax="$[ ${clat} + ${3} ]"
export lonmin="$[ ${clon} - ${4} ]"
export lonmax="$[ ${clon} + ${5} ]"

export n1="${7}"
export n2="${8}"

##############################################################################

if [ ! -f $outdir ]; then

    ncl ${wkdir}/script/Jet_indices.ncl
fi
echo "jet indices finish."

ncl  ${wkdir}/script/plot_backup.ncl

echo "plot finish."
