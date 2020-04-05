#!/bin/bash
./run_Jupyter_notebook.sh Tutorial-Finite_Difference_Derivatives.ipynb && git diff Tutorial-Finite_Difference_Derivatives.ipynb && \
./run_Jupyter_notebook.sh Tutorial-Numerical_Grids.ipynb && git diff Tutorial-Numerical_Grids.ipynb && \
./run_Jupyter_notebook.sh Tutorial-Coutput__Parameter_Interface.ipynb && git diff Tutorial-Coutput__Parameter_Interface.ipynb && \
./run_Jupyter_notebook.sh Tutorial-cmdline_helper.ipynb && git diff Tutorial-cmdline_helper.ipynb && \
./run_Jupyter_notebook.sh Tutorial-Symbolic_Tensor_Rotation.ipynb && git diff Tutorial-Symbolic_Tensor_Rotation.ipynb
