#! /bin/bash

  export sdir="/home/sunming/data5/cuixy/global_jets/script"

# ulow tlow
  echo "ulow,tlow begin."
  ncl -Q 'tin="low"' 'uin="low"' 'varname="ultl"' $sdir/calc_budget_terms_seperate.ncl
  echo "ulow,tlow end."

# ulow tiso
  echo "ulow,tiso begin."
  ncl -Q 'tin="iso"' 'uin="low"' 'varname="ulti"' $sdir/calc_budget_terms_seperate.ncl
  echo "ulow,tiso end."

# ulow thigh
  echo "ulow,thigh begin."
  ncl -Q 'tin="high"' 'uin="low"' 'varname="ulth"' $sdir/calc_budget_terms_seperate.ncl
  echo "ulow,thigh end."

# uiso tlow
  echo "uiso,tlow begin."
  ncl -Q 'tin="low"' 'uin="iso"' 'varname="uitl"' $sdir/calc_budget_terms_seperate.ncl
  echo "uiso,tlow end."

# uiso tiso
  echo "uiso,tiso begin."
  ncl -Q 'tin="iso"' 'uin="iso"' 'varname="uiti"' $sdir/calc_budget_terms_seperate.ncl
  echo "uiso,tiso end."

# uiso thigh
  echo "uiso,thigh begin."
  ncl -Q 'tin="high"' 'uin="iso"' 'varname="uith"' $sdir/calc_budget_terms_seperate.ncl
  echo "uiso,thigh end."

# uhigh tlow
  echo "uhigh,tlow begin."
  ncl -Q 'tin="low"' 'uin="high"' 'varname="uhtl"' $sdir/calc_budget_terms_seperate.ncl
  echo "uhigh,tlow end."

# uhigh tiso
  echo "uhigh,tiso begin."
  ncl -Q 'tin="iso"' 'uin="high"' 'varname="uhti"' $sdir/calc_budget_terms_seperate.ncl
  echo "uhigh,tiso end."

# uhigh thigh
  echo "uhigh,thigh begin."
  ncl -Q 'tin="high"' 'uin="high"' 'varname="uhth"' $sdir/calc_budget_terms_seperate.ncl
  echo "uhigh,thigh end."
