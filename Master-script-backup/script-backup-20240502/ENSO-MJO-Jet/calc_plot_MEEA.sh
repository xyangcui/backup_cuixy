#! /bin/bash
###############################################################################
# This script calls for two ncl scripts.
# 1. jet_indices.ncl to calculate jet indices.
# 2. plot.ncl to plot
# $1: jet location; $2: latmin $3: latmax $4: lonmin $5: lonmax $6: box.
# Information
# author: cui xiangyang time: 2023-11-07
###############################################################################

export wkdir="/home/sunming/data5/cuixy/ENSO-MJO-Jet"
export indir="${wkdir}/data/jet_indices_${1}.nc"
export outdir="${wkdir}/data/monte_carlo/jet_indices_${1}.nc"
export rmmdir="/home/sunming/data5/cuixy/data/clim_indices/rmm.74toRealtime.csv"
export ninodir="/home/sunming/data5/cuixy/data/clim_indices/Nino_3.4_index.csv"

export var1="jets"
export var2="jetl"

export n1="${2}"
export n2="${3}"

##############################################################################

ncl -Q ${wkdir}/script/calc_stat.ncl
