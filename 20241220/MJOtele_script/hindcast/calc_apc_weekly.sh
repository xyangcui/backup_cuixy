#! /bin/bash
########################################################### 
  export dir="/home/sunming/data5/cuixy/MJO_PC_diversiy/" 
  years=("2021" "2022" "2023" "2019" "2018" "2017" "2016")

  #/home/app/anaconda3/envs/ncl/bin/ncl -Q 'Year="2020"' ${dir}/script/hindcast/select_apc_days.ncl
  for i in {1..8}
  do 
    /home/app/anaconda3/envs/ncl/bin/ncl -Q i=$i ${dir}/script/hindcast/calc_apc_ens_weekly.ncl
  done