set term post enh color "Helvetica" fontscale 1.5
#set encoding utf8

set out "num-ex.ps"
set xlabel "x location"
set ylabel "(Numerical - Exact) solution"
set title "Data at time 40.0"
p [-10:10] "simplewave_sine-0.1-1D/uuGF.xl-parsed"  u 1:($2-sin($1-$3)) i 401 w lp ti "Delta x = 0.1", \
  	   "simplewave_sine-0.05-1D/uuGF.xl-parsed" u 1:($2-sin($1-$3)) i 801 w lp ti "Delta x = 0.05"

#p [-10:10] "simplewave_sine-0.1-1D/uuGF.xl-parsed" i 801 u 1:($2-sin($1-$3)) w l, \
#  "simplewave_sine-0.05-1D/uuGF.xl-parsed" i 401 u 1:($2-sin($1-$3)) w l title "dx = 0.05"
!ps2pdf13 num-ex.ps

set out "num-ex_scaled.ps"
p [-10:10] "simplewave_sine-0.1-1D/uuGF.xl-parsed"  u 1:($2-sin($1-$3))       i 401 w lp ti "Delta x = 0.1", \
  	   "simplewave_sine-0.05-1D/uuGF.xl-parsed" u 1:(($2-sin($1-$3))*16.) i 801 w lp ti "Delta x = 0.05"
!ps2pdf13 num-ex_scaled.ps
