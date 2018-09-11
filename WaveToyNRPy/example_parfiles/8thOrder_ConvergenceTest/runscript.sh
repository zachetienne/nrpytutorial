#!/bin/bash
for i in planewave_along_3D_diagonal-dx_0.4__FD8-RK8.par planewave_along_3D_diagonal-dx_0.2__FD8-RK8.par; do
	taskset -c 0,1,2,3 ./cactus_etilgrmhd-FD8 $i
	DIRNAME=`echo $i|sed "s/.par//g"`
	./convert_IOASCII_1D_to_gnuplot.sh $DIRNAME
done
