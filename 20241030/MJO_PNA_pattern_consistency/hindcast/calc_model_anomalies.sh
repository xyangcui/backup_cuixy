#! /bin/bash
########################################################### 
  export dir="/home/sunming/data5/cuixy/MJO_PC_diversiy/" 
  models=("ncep" "cma" "isac" "ecmwf")
  for model in ${models[*]}
  do 
    echo "$model begin. $(date "+%Y-%m-%d %H:%M:%S")"
    /home/sunming/miniconda3/envs/py39/bin/python ${dir}/script/hindcast/calc_model_anomalies.py ${model}
    echo "$model done. $(date "+%Y-%m-%d %H:%M:%S")"
  done
