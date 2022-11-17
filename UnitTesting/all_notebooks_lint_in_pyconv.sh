#!/bin/bash

rm -rf pyconv
mkdir pyconv

rm -f /tmp/joblist.txt /tmp/joblist2.txt
for i in *.ipynb; do
    # FIXME: Don't know why Tutorial-TOV-Piecewise_Polytrope_EOSs.ipynb.py causes pylint to choke...
    # FIXME: Don't know why Tutorial-WeylScalarsInvariants-Cartesian.ipynb causes pylint to choke...
    if [ $i != "Tutorial-BaikalETK.ipynb" ] && [ $i != "NRPyPlus_Tutorial.ipynb" ] && [ $i != "Tutorial-How_NRPy_Computes_Finite_Difference_Coeffs.ipynb" ] && [ $i != "Tutorial-TOV-Piecewise_Polytrope_EOSs.ipynb" ] && [ $i != "Tutorial-WeylScalarsInvariants-Cartesian.ipynb" ] ; then
       echo "jupyter nbconvert --to python $i --stdout |grep -v \"^\# \" > pyconv/$i.py" >> /tmp/joblist.txt  # ignore lines that start with #.
       echo "echo \"Linting complete: $i\" ; pylint --disable=trailing-newlines,reimported,ungrouped-imports --output=pyconv/$i.py.txt pyconv/$i.py" >> /tmp/joblist2.txt
    fi
done

parallel --jobs 16 < /tmp/joblist.txt
parallel --jobs 16 < /tmp/joblist2.txt

echo "TOP OFFENDERS:"
for i in pyconv/*.txt; do echo `grep rated $i` $i;done |sort -k7 -g|head -n20
