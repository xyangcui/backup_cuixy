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
  export indir="${wkdir}/data/u250_daily_mean_1979-2022.nc"
  export outdir="${wkdir}/data/jet_indices_${1}.nc"

  export lats="$[  ${2} ]"
  export latn="$[  ${3} ]"
  export lonw="$[  ${4} ]"
  export lone="$[  ${5} ]"

  echo "jet indices begin."
  
  if [ ! -f $outdir ]; then
    ncl ${wkdir}/script/Jet_indices.ncl
  fi

  echo "jet indices done."

