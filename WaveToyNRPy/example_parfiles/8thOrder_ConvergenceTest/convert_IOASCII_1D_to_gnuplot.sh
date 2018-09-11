#!/bin/bash

cat $1/uuGF.xl | awk '{if($1=="#Time") printf("\n"); print $0}' | awk -v PREC="quad" '{if(NF==3) { printf("%s %.15e\n",$0,sin($2/sqrt(3.)-$1)); } else { print $0 }}' > $1/uuGF.gnuplot
