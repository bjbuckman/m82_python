#!/bin/bash
# Reruns least squares

cd ../gp_out

for dir in m82_05*w*
do
	gp_dir=$dir
	cd ${gp_dir}
	python ../../plots/gp2d_leastsquares_v4.py ${gp_dir}
	echo ${gp_dir}
	cd ../
done
