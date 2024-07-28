#! /bin/bash

  export wkdir="/home/sunming/data5/cuixy/ENSO-MJO-Jet/script"


  echo "Neutral year begin."
  ncl -Q 'nty="NE"' ${wkdir}/calc_MLRegression_MJOJet.ncl

  echo "EI year begin."
  ncl -Q 'nty="EI"' ${wkdir}/calc_MLRegression_MJOJet.ncl

  echo "LA year begin."
  ncl -Q 'nty="LA"' ${wkdir}/calc_MLRegression_MJOJet.ncl