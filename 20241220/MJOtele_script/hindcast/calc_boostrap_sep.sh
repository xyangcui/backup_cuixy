#! /bin/bash
########################################################### 
  export dir="/home/sunming/data5/cuixy/MJO_PC_diversiy" 

  echo week1 begin. `date`
  /home/app/anaconda3/envs/ncl/bin/ncl -Q ${dir}/script/hindcast/calc_boostrap_apc_sep_models_sp_w1.ncl
  echo week1 done. `date`

  echo week2 begin. `date`
  /home/app/anaconda3/envs/ncl/bin/ncl -Q ${dir}/script/hindcast/calc_boostrap_apc_sep_models_sp_w2.ncl
  echo week2 done. `date`

  echo week3 begin. `date`
  /home/app/anaconda3/envs/ncl/bin/ncl -Q ${dir}/script/hindcast/calc_boostrap_apc_sep_models_sp_w3.ncl
  echo week3 done. `date`

  echo week4 begin. `date`
  /home/app/anaconda3/envs/ncl/bin/ncl -Q ${dir}/script/hindcast/calc_boostrap_apc_sep_models_sp_w4.ncl
  echo week4 done. `date`