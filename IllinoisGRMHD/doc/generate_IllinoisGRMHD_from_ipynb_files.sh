#!/bin/bash

for i in *.ipynb; do
    echo Executing $i ...
    jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=-1 $i
done

echo Finished!
