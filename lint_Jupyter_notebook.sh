#!/bin/bash

if [ -z "$1" ]; then
    echo "Correct usage of this script:"
    echo "./lint_Jupyter_notebook.sh [Jupyter notebook (a file with .ipynb extension)]"
    exit
fi

if [ ! -f $1 ]; then
    echo "You input: ./lint_Jupyter_notebook.sh" $1
    echo "Jupyter notebook" \"$1\" "not found!"
    exit
fi

jupyter nbconvert --to python $1 --stdout |grep -v "^\#" > $1.py # ignore lines that start with #.
pylint --disable=trailing-newlines,reimported,ungrouped-imports $1.py
