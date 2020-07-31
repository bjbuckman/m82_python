#!/bin/bash
# Reruns least squares

cd ../gp_out

for dir in m82_02*  #{10..18}??
do
	gp_dir=$dir
	cd ${gp_dir}
	python ../../plots/gp2d_leastsquares_ff_v2.py ${gp_dir}
	echo ${gp_dir}
	cd ..
done
