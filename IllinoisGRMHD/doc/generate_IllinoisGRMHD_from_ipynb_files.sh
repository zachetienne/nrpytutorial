#!/bin/bash

NUM_JOBS=0

count=0
mkdir ../src
mkdir ../Convert_to_HydroBase
mkdir ../Convert_to_HydroBase/src
mkdir ../ID_converter_ILGRMHD
mkdir ../ID_converter_ILGRMHD/src
for i in *.ipynb; do
    echo Executing $i ...
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

    jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=-1 $i &
    if ((count==$NUM_JOBS)); then
        wait
        count=0
    else
        let count+=1
    fi
done

wait
echo Finished!

# Alternative approach if gnu parallel is installed:
# rm -f /tmp/joblist.txt
# for i in *.ipynb; do
#     echo jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=-1 $i >> /tmp/joblist.txt
# done
# parallel --jobs 8 < /tmp/joblist.txt
