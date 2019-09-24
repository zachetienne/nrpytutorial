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

jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=-1 $1
