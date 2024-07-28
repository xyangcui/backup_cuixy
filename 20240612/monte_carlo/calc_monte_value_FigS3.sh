#! /bin/bash

  export wkdir="/home/sunming/data5/cuixy/global_jets/script/monte_carlo"

  echo "calc meri monte value begin. $(date)"
  ncl -Q ${wkdir}/calc_monte_value_meri_FigS3.ncl
  echo "calc meri monte value done. $(date)"


  echo "calc olr monte value begin. $(date)"
  ncl -Q ${wkdir}/calc_monte_value_olr_FigS3.ncl
  echo "calc olr monte value done. $(date)"