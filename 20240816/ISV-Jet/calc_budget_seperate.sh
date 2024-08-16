#! /bin/bash

  export sdir="/home/sunming/data5/cuixy/global_jets/script"

# ulow tlow
  echo "ulow,tlow begin."
  ncl -Q 'varname="ultl"' $sdir/calc_budget_terms_seperate.ncl
  echo "ulow,tlow end."

# ulow tiso
  echo "ulow,tiso begin."
  ncl -Q 'varname="ulti"' $sdir/calc_budget_terms_seperate.ncl
  echo "ulow,tiso end."

# ulow thigh
  echo "ulow,thigh begin."
  ncl -Q 'varname="ulth"' $sdir/calc_budget_terms_seperate.ncl
  echo "ulow,thigh end."

# uiso tlow
  echo "uiso,tlow begin."
  ncl -Q 'varname="uitl"' $sdir/calc_budget_terms_seperate.ncl
  echo "uiso,tlow end."

# uiso tiso
  echo "uiso,tiso begin."
  ncl -Q 'varname="uiti"' $sdir/calc_budget_terms_seperate.ncl
  echo "uiso,tiso end."

# uiso thigh
  echo "uiso,thigh begin."
  ncl -Q 'varname="uith"' $sdir/calc_budget_terms_seperate.ncl
  echo "uiso,thigh end."

# uhigh tlow
  echo "uhigh,tlow begin."
  ncl -Q 'varname="uhtl"' $sdir/calc_budget_terms_seperate.ncl
  echo "uhigh,tlow end."

# uhigh tiso
  echo "uhigh,tiso begin."
  ncl -Q 'varname="uhti"' $sdir/calc_budget_terms_seperate.ncl
  echo "uhigh,tiso end."

# uhigh thigh
  echo "uhigh,thigh begin."
  ncl -Q 'varname="uhth"' $sdir/calc_budget_terms_seperate.ncl
  echo "uhigh,thigh end."


### vdtdy
# ulow tlow
  echo "ulow,tlow begin."
  ncl -Q 'varname="vltl"' $sdir/calc_budget_terms_seperate-v.ncl
  echo "ulow,tlow end."

# ulow tiso
  echo "ulow,tiso begin."
  ncl -Q 'varname="vlti"' $sdir/calc_budget_terms_seperate-v.ncl
  echo "ulow,tiso end."

# ulow thigh
  echo "ulow,thigh begin."
  ncl -Q 'varname="vlth"' $sdir/calc_budget_terms_seperate-v.ncl
  echo "ulow,thigh end."

# uiso tlow
  echo "uiso,tlow begin."
  ncl -Q 'varname="vitl"' $sdir/calc_budget_terms_seperate-v.ncl
  echo "uiso,tlow end."

# uiso tiso
  echo "uiso,tiso begin."
  ncl -Q 'varname="viti"' $sdir/calc_budget_terms_seperate-v.ncl
  echo "uiso,tiso end."

# uiso thigh
  echo "uiso,thigh begin."
  ncl -Q 'varname="vith"' $sdir/calc_budget_terms_seperate-v.ncl
  echo "uiso,thigh end."

# uhigh tlow
  echo "uhigh,tlow begin."
  ncl -Q 'varname="vhtl"' $sdir/calc_budget_terms_seperate-v.ncl
  echo "uhigh,tlow end."

# uhigh tiso
  echo "uhigh,tiso begin."
  ncl -Q 'varname="vhti"' $sdir/calc_budget_terms_seperate-v.ncl
  echo "uhigh,tiso end."

# uhigh thigh
  echo "uhigh,thigh begin."
  ncl -Q 'varname="vhth"' $sdir/calc_budget_terms_seperate-v.ncl
  echo "uhigh,thigh end."

