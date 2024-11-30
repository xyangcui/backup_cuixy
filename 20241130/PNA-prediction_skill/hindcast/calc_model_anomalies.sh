#! /bin/bash
########################################################### 
  export dir="/home/sunming/data5/cuixy/MJO_PC_diversiy" 
  year=("2021" "2022" "2023" "2019" "2018" "2017" "2016")

 
  #echo "${year[0]}, ${year[1]} and ${year[2]} begin. $(date "+%Y-%m-%d %H:%M:%S")"
  #/home/sunming/miniconda3/envs/py39/bin/python ${dir}/script/hindcast/calc_model_anomalies.py ${year[0]} & \
  #  /home/sunming/miniconda3/envs/py39/bin/python ${dir}/script/hindcast/calc_model_anomalies.py ${year[1]} & \
  #    /home/sunming/miniconda3/envs/py39/bin/python ${dir}/script/hindcast/calc_model_anomalies.py ${year[2]}
  #echo "${year[0]}, ${year[1]} and ${year[2]} done. $(date "+%Y-%m-%d %H:%M:%S")"

  echo "${year[3]} and ${year[4]} begin. $(date "+%Y-%m-%d %H:%M:%S")"
  /home/sunming/miniconda3/envs/py39/bin/python ${dir}/script/hindcast/calc_model_anomalies.py ${year[3]} & \
    /home/sunming/miniconda3/envs/py39/bin/python ${dir}/script/hindcast/calc_model_anomalies.py ${year[4]} 
  echo "${year[3]} and ${year[4]} done. $(date "+%Y-%m-%d %H:%M:%S")"


  echo "${year[5]} and ${year[6]} begin. $(date "+%Y-%m-%d %H:%M:%S")"
  /home/sunming/miniconda3/envs/py39/bin/python ${dir}/script/hindcast/calc_model_anomalies.py ${year[5]} & \
    /home/sunming/miniconda3/envs/py39/bin/python ${dir}/script/hindcast/calc_model_anomalies.py ${year[6]}
  echo "${year[5]} and ${year[6]} done. $(date "+%Y-%m-%d %H:%M:%S")"

  exit 0