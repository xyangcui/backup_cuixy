#! /bin/bash

export sdir="/home/sunming/data5/cuixy/global_jets/script"

for j in $( seq -w 1979 2022 );do

    echo "calc year ${j} begin." $(date)

	ncl -Q j=${j} nwa=100. nwb=20. 'var="t"'\
	 ${sdir}/calc_filtered_levelvalues.ncl

	echo "calc year ${j} done." $(date)

done