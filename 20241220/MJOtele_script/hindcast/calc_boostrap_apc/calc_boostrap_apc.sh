#! /bin/bash
########################################################### 
  export dir="/home/sunming/data5/cuixy/MJO_PC_diversiy" 

  /home/app/anaconda3/envs/ncl/bin/ncl -Q ${dir}/script/hindcast/calc_boostrap_apc/calc_boostrap_apc_w1.ncl

  /home/app/anaconda3/envs/ncl/bin/ncl -Q ${dir}/script/hindcast/calc_boostrap_apc/calc_boostrap_apc_w2.ncl

  /home/app/anaconda3/envs/ncl/bin/ncl -Q ${dir}/script/hindcast/calc_boostrap_apc/calc_boostrap_apc_w3.ncl

  /home/app/anaconda3/envs/ncl/bin/ncl -Q ${dir}/script/hindcast/calc_boostrap_apc/calc_boostrap_apc_w4.ncl
