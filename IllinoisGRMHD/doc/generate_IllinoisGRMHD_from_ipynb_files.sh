#!/bin/bash

NUM_JOBS=4

count=0
for i in *.ipynb; do
    echo Executing $i ...
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
