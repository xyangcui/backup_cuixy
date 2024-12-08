#! /bin/bash
########################################################### 
  export dir="/home/sunming/data5/cuixy/MJO_PC_diversiy/" 
  years=("2021" "2022" "2023" "2019" "2018" "2017" "2016")

  /home/app/anaconda3/envs/ncl/bin/ncl -Q 'Year="2020"' ${dir}/script/hindcast/apc_all_events/select_apc_days.ncl
  for i in {1..8}
  do 
    /home/app/anaconda3/envs/ncl/bin/ncl -Q i=$i 'Year="2020"' ${dir}/script/hindcast/apc_all_events/calc_apc_ens.ncl
  done

#Year 2021
  /home/app/anaconda3/envs/ncl/bin/ncl -Q 'Year="2021"' ${dir}/script/hindcast/apc_all_events/select_apc_days.ncl
  for i in {1..8}
  do 
    /home/app/anaconda3/envs/ncl/bin/ncl -Q i=$i 'Year="2021"' ${dir}/script/hindcast/apc_all_events/calc_apc_ens.ncl
  done

#Year 2022
  /home/app/anaconda3/envs/ncl/bin/ncl -Q 'Year="2022"' ${dir}/script/hindcast/apc_all_events/select_apc_days.ncl
  for i in {1..8}
  do 
    /home/app/anaconda3/envs/ncl/bin/ncl -Q i=$i 'Year="2022"' ${dir}/script/hindcast/apc_all_events/calc_apc_ens.ncl
  done

#Year 2023
  /home/app/anaconda3/envs/ncl/bin/ncl -Q 'Year="2023"' ${dir}/script/hindcast/apc_all_events/select_apc_days.ncl
  for i in {1..8}
  do 
    /home/app/anaconda3/envs/ncl/bin/ncl -Q i=$i 'Year="2023"' ${dir}/script/hindcast/apc_all_events/calc_apc_ens.ncl
  done

#Year 2019
  /home/app/anaconda3/envs/ncl/bin/ncl -Q 'Year="2019"' ${dir}/script/hindcast/apc_all_events/select_apc_days.ncl
  for i in {1..8}
  do 
    /home/app/anaconda3/envs/ncl/bin/ncl -Q i=$i 'Year="2019"' ${dir}/script/hindcast/apc_all_events/calc_apc_ens.ncl
  done

#Year 2018
  /home/app/anaconda3/envs/ncl/bin/ncl -Q 'Year="2018"' ${dir}/script/hindcast/apc_all_events/select_apc_days.ncl
  for i in {1..8}
  do 
    /home/app/anaconda3/envs/ncl/bin/ncl -Q i=$i 'Year="2018"' ${dir}/script/hindcast/apc_all_events/calc_apc_ens.ncl
  done

#Year 2017
  /home/app/anaconda3/envs/ncl/bin/ncl -Q 'Year="2017"' ${dir}/script/hindcast/apc_all_events/select_apc_days.ncl
  for i in {1..8}
  do 
    /home/app/anaconda3/envs/ncl/bin/ncl -Q i=$i 'Year="2017"' ${dir}/script/hindcast/apc_all_events/calc_apc_ens.ncl
  done

#Year 2016
  /home/app/anaconda3/envs/ncl/bin/ncl -Q 'Year="2016"' ${dir}/script/hindcast/apc_all_events/select_apc_days.ncl
  for i in {1..8}
  do 
    /home/app/anaconda3/envs/ncl/bin/ncl -Q i=$i 'Year="2016"' ${dir}/script/hindcast/apc_all_events/calc_apc_ens.ncl
  done
