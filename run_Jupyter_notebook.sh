#!/bin/bash

if [ -z "$1" ]; then
    echo "Correct usage of this script:"
    echo "./run_Jupyter_notebook.sh [Jupyter notebook (a file with .ipynb extension)]"
    exit
fi

if [ ! -f $1 ]; then
    echo "You input: ./run_Jupyter_notebook.sh" $1
    echo "Jupyter notebook" \"$1\" "not found!"
    exit
fi

# NRPy+ Jupyter notebooks are completely Python 2/3 cross-compatible.
#   However `jupyter nbconvert` will refuse to run if the notebook
#   was generated using a different kernel. Here we fool Jupyter
#   to think the notebook was written using the native python kernel.
PYTHONMAJORVERSION=`python -c "import sys;print(sys.version_info[0])"`
if (( $PYTHONMAJORVERSION == 3 )); then
    cat $1 | sed "s/   \"name\": \"python2\"/   \"name\": \"python3\"/g" > $1-tmp ; mv $1-tmp $1
else
    cat $1 | sed "s/   \"name\": \"python3\"/   \"name\": \"python2\"/g" > $1-tmp ; mv $1-tmp $1
fi

time jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=-1 $1
