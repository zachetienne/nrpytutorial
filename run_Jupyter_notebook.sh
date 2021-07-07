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

if [ "$2" == "notimer" ]; then
    if jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=-1 $1;
    then
        echo
        # do nothing
    else
        echo BACKTRACE:
        echo git diff $1
        cat in_progress/Validation/out_GiRaFFEfood_NRPy_test.txt
        git diff $1 2>&1 | cat
        exit 1
    fi
else
    if time jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=-1 $1;
    then
        echo
        # do nothing
    else
        echo BACKTRACE:
        if jupyter nbconvert --to python $1;
        then
            FILENAME=`echo $1 | sed 's/ipynb/py/g'`
            echo $FILENAME
            if (( $PYTHONMAJORVERSION == 3 )); then
                ipython3 --log-level=DEBUG $FILENAME
            else
                ipython --log-level=DEBUG $FILENAME
            fi
            exit 1
        else
            echo ERROR: could not convert $1 to a Python script!
            exit 1
        fi
    fi
fi
