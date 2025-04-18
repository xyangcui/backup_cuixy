#! /bin/bash

/home/sunming/miniconda3/envs/py39/bin/python PrcAnalysis.py

if [ -e "/home/sunming/data5/cuixy/Subpre_NPJ/data/PrCA_patterns.nc" ]; then

	echo "PrC_analysis ran successfully, next to get obs invariates."

	/home/sunming/miniconda3/envs/py39/bin/python PrCA_obs_variates_all.py
	/home/sunming/miniconda3/envs/py39/bin/python PrCA_obs_variates.py

else
	echo "PrC_analysis ran unsuccessfully, next to retry."

	/home/sunming/miniconda3/envs/py39/bin/python PrcAnalysis.py

fi
