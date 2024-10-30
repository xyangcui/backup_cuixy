#! /bin/bash

  export sdir="/home/sunming/data5/cuixy/MJO_PC_diversiy/script"

  echo "w begin."
  ncl -Q 'var="w"' ${sdir}/calc_filtered_vvalues.ncl
  echo "w done."

  echo "u begin."
  ncl -Q 'var="u"' ${sdir}/calc_filtered_vvalues.ncl
  echo "u done."

  echo "q begin."
  ncl -Q 'var="q"' ${sdir}/calc_filtered_vvalues.ncl
  echo "q done."

  exit 0