#! /bin/bash

ncl -Q 'var="dtdt"' ./calc_filtered_temp_budget.ncl

ncl -Q 'var="udtdx"' ./calc_filtered_temp_budget.ncl

ncl -Q 'var="vdtdy"' ./calc_filtered_temp_budget.ncl

ncl -Q 'var="wdtdp"' ./calc_filtered_temp_budget.ncl
