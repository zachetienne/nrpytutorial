#!/bin/bash

set -e # Error out if any commands complete with an error.

# Skip Baikal and all Start-to-Finish notebooks, except ScalarWave for now
# Tutorial-Start_to_Finish-ScalarWave*.ipynb
# Let's try all but Psi4 Start-to-Finish and Baikal notebooks. Also cmdlinehelper yields whitespace differences in Python 2.7
for i in Tutorial-[A]*.ipynb Tutorial-[C-RT-Z]*.ipynb Tutorial-B[B-Z]*.ipynb Tutorial-S[A-SU-Z]*.ipynb NRPyPN/PN*.ipynb; do
    # For some reason (as of ~July 20, 2020) the hydro-without-hydro notebook takes too long in Travis, causing a timeout:
    #   Also as of Aug 3, 2020 the new WaveToyNRPy notebook is broken, as it seems Travis doesn't support parallel codegens
    if [ $i != "Tutorial-Start_to_Finish-BSSNCurvilinear-Neutron_Star-Hydro_without_Hydro.ipynb" ] && [ $i != "Tutorial-ETK_thorn-WaveToyNRPy.ipynb" ]; then
        ./run_Jupyter_notebook.sh $i notimer
        cat $i | sed "s/\\\r\\\n/\\\n/g" > $i-new && mv $i-new $i
        git diff $i |grep -v "image/png"|grep -E "^\-|^\+"|grep -v  '^\-\-\-'| \
            grep -v "metadata\":"|grep -v "\"execution\":"|grep -v "\"iopub."| \
            grep -v "\"shell.execute"|grep -v "\"version\":"|grep -v "   }"$|grep -v "   },"$| cdiff |cat
        #    git diff $i | grep -v "image/png" | cdiff | cat
        #    echo Number of lines different in the git diff: `git diff|grep -v image/png|wc -l`
    fi
done

# GiRaFFE unit tests:
# for i in in_progress/Tutorial-Start_to_Finish-GiRaFFE_NRPy-1D_tests-staggered.ipynb in_progress/Tutorial-Start_to_Finish_UnitTest*; do
#     ./run_Jupyter_notebook.sh $i notimer
#     cat $i | sed "s/\\\r\\\n/\\\n/g" > $i-new && mv $i-new $i
#     git diff $i |grep -v "image/png"|grep -E "^\-|^\+"|grep -v  '^\-\-\-'|cdiff |cat
# done
# ./run_Jupyter_notebook.sh Tutorial-Finite_Difference_Derivatives.ipynb && git diff Tutorial-Finite_Difference_Derivatives.ipynb && \
# ./run_Jupyter_notebook.sh Tutorial-Numerical_Grids.ipynb && git diff Tutorial-Numerical_Grids.ipynb && \
# ./run_Jupyter_notebook.sh Tutorial-Coutput__Parameter_Interface.ipynb && git diff Tutorial-Coutput__Parameter_Interface.ipynb && \
# ./run_Jupyter_notebook.sh Tutorial-cmdline_helper.ipynb && git diff Tutorial-cmdline_helper.ipynb && \
# ./run_Jupyter_notebook.sh Tutorial-Symbolic_Tensor_Rotation.ipynb && git diff Tutorial-Symbolic_Tensor_Rotation.ipynb
