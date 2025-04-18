#! /bin/bash

#ncl -Q 'model="cma"'   prepare_for_acc_anom.ncl 
#ncl -Q 'model="cnrm"'  prepare_for_acc_anom.ncl 
#ncl -Q 'model="eccc"'  prepare_for_acc_anom.ncl 
#ncl -Q 'model="ecmwf"' prepare_for_acc_anom.ncl 
#ncl -Q 'model="hmcr"'  prepare_for_acc_anom.ncl 
#ncl -Q 'model="isac"'  prepare_for_acc_anom.ncl 
#ncl -Q 'model="kma"'   prepare_for_acc_anom.ncl 
#ncl -Q 'model="ncep"'  prepare_for_acc_anom.ncl 
#ncl -Q 'model="ukmo"'  prepare_for_acc_anom.ncl 

ncl -Q 'model="cma"'   calc_acc_alldays.ncl 
ncl -Q 'model="cnrm"'  calc_acc_alldays.ncl 
ncl -Q 'model="eccc"'  calc_acc_alldays.ncl 
ncl -Q 'model="ecmwf"' calc_acc_alldays.ncl 
ncl -Q 'model="hmcr"'  calc_acc_alldays.ncl 
ncl -Q 'model="isac"'  calc_acc_alldays.ncl
ncl -Q 'model="kma"'   calc_acc_alldays.ncl 
ncl -Q 'model="ncep"'  calc_acc_alldays.ncl 
ncl -Q 'model="ukmo"'  calc_acc_alldays.ncl 
ncl -Q 'model="LIM"'   calc_acc_alldays.ncl 

#ncl -Q 'model="cma"'   calc_acc_index.ncl 
#ncl -Q 'model="cnrm"'  calc_acc_index.ncl 
#ncl -Q 'model="eccc"'  calc_acc_index.ncl 
#cl -Q 'model="ecmwf"' calc_acc_index.ncl 
#cl -Q 'model="hmcr"'  calc_acc_index.ncl 
#ncl -Q 'model="isac"'  calc_acc_index.ncl
#ncl -Q 'model="kma"'   calc_acc_index.ncl 
#ncl -Q 'model="ncep"'  calc_acc_index.ncl 
#ncl -Q 'model="ukmo"'  calc_acc_index.ncl 
#ncl -Q 'model="LIM"'   calc_acc_index.ncl 