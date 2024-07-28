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
export indir="${wkdir}/data/jet_indices_${1}.nc"
export outdir="${wkdir}/data/monte_carlo/jet_indices_${1}.nc"
export rmmdir="/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
export udir="${wkdir}/data/ERA5_daily_u250_1979-2022.nc"
export plotdir="${wkdir}/plot/MEJS"

export var1="jets"
export var2="jetl"

export n1="${6}"
export n2="${7}"

export latmin="$[  ${2} ]"
export latmax="$[  ${3} ]"
export lonmin="$[  ${4} ]"
export lonmax="$[  ${5} ]"

##############################################################################

ncl -Q ${wkdir}/script/monte_carlo/calc_stat.ncl

ncl -Q ${wkdir}/script/monte_carlo/plot_backup.ncl
