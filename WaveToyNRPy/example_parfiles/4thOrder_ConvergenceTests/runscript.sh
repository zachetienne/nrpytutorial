#!/bin/bash
for i in planewave_along_3D_diagonal-dx_0.4__FD4-RK4.par planewave_along_3D_diagonal-dx_0.2__FD4-RK4.par 1D-planewave-dx_0.4__FD4-RK4.par 1D-planewave-dx_0.2__FD4-RK4.par; do
	taskset -c 0,1,2,3 ./cactus_etilgrmhd-FD4 $i
	DIRNAME=`echo $i|sed "s/.par//g"`
	./convert_IOASCII_1D_to_gnuplot.sh $DIRNAME
done
