#! /bin/bash

export wkdir1="/home/sunming/data5/cuixy/global_jets"

#bash ${wkdir1}/script/calc_plot.sh NAJ 25 60 -60 0 box2 100 20
#echo "NAJ_box2_20-100"
#bash ${wkdir1}/script/calc_plot.sh NAJ 25 60 -60 0 box2 60  30
#echo "NAJ_box2_30-60"

#bash ${wkdir1}/script/calc_plot.sh NAJ 15 75 -90 0 box3 100 20
#echo "NAJ_box3_20-100"

#bash ${wkdir1}/script/calc_plot.sh NAJ 15 75 -90 0 box3 60  30
#echo "NAJ_box3_30-60"

bash ${wkdir1}/script/calc_plot.sh NAJ 25 60 -60 0 box2 60  30
echo "NAJ_box2"

bash ${wkdir1}/script/calc_plot.sh NAJ 15 75 -90 0 box3 100 20
echo "NAJ_box3"