#!/bin/bash

rm -f /tmp/outNRPybench.txt
time for i in *.ipynb; do

    # First clean up any mess made by notebook.
    git clean -fdq
    echo Working on $i now ...

    # NRPy+ Jupyter notebooks are completely Python 2/3 cross-compatible.
    #   However `jupyter nbconvert` will refuse to run if the notebook
    #   was generated using a different kernel. Here we fool Jupyter
    #   to think the notebook was written using the native python kernel.
    PYTHONMAJORVERSION=`python -c "import sys;print(sys.version_info[0])"`
    if (( $PYTHONMAJORVERSION == 3 )); then
        cat $i | sed "s/   \"name\": \"python2\"/   \"name\": \"python3\"/g" > $i-tmp ; mv $i-tmp $i
    else
        cat $i | sed "s/   \"name\": \"python3\"/   \"name\": \"python2\"/g" > $i-tmp ; mv $i-tmp $i
    fi

    BENCH=$(/usr/bin/time -f %e jupyter nbconvert --log-level=0 --to notebook --inplace --execute --ExecutePreprocessor.timeout=-1 $i 2>&1)
    echo $BENCH $i | tee -a /tmp/outNRPybench.txt

done

sort -k1 -g /tmp/outNRPybench.txt
# Clean up any mess made by last notebook run
git clean -fdq
