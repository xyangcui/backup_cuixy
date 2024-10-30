#! /bin/bash
########################################################### 
  export dir="/home/sunming/data5/cuixy/MJO_PC_diversiy/" 
  models=("ncep" "cma" "isac" "ecmwf")

  /home/app/anaconda3/envs/ncl/bin/ncl -Q 'model_name="isac"' ${dir}/script/hindcast/select_apc_days.ncl
  for i in {1..8}
  do 
    /home/app/anaconda3/envs/ncl/bin/ncl -Q i=$i 'model_name="isac"' ${dir}/script/hindcast/calc_apc.ncl
  done

  /home/app/anaconda3/envs/ncl/bin/ncl -Q 'model_name="ncep"' ${dir}/script/hindcast/select_apc_days.ncl
  for i in {1..8}
  do 
    /home/app/anaconda3/envs/ncl/bin/ncl -Q i=$i 'model_name="ncep"' ${dir}/script/hindcast/calc_apc.ncl
  done

  /home/app/anaconda3/envs/ncl/bin/ncl -Q 'model_name="cma"' ${dir}/script/hindcast/select_apc_days.ncl
  for i in {1..8}
  do 
    /home/app/anaconda3/envs/ncl/bin/ncl -Q i=$i 'model_name="cma"' ${dir}/script/hindcast/calc_apc.ncl
  done

  /home/app/anaconda3/envs/ncl/bin/ncl -Q 'model_name="ecmwf"' ${dir}/script/hindcast/select_apc_days.ncl
  for i in {1..8}
  do 
    /home/app/anaconda3/envs/ncl/bin/ncl -Q i=$i 'model_name="ecmwf"' ${dir}/script/hindcast/calc_apc.ncl
  done
