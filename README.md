# NRPy+, SENRv2, and the NRPy+ Jupyter Tutorial
[![Build Status](https://travis-ci.org/zachetienne/nrpytutorial.svg?branch=master)](https://travis-ci.org/zachetienne/nrpytutorial)
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/zachetienne/nrpytutorial/master?filepath=NRPyPlus_Tutorial.ipynb)

This repository houses
* [The newest version of NRPy+: Python-Based Code Generation for Numerical Relativity... and Beyond](https://arxiv.org/abs/1712.07658),
* The second version of SENR, the Simple, Efficient Numerical Relativity code (see the "Colliding Black Holes" Start-to-Finish tutorial module), and 
* The NRPy+ Jupyter Tutorial: An Introduction to Python-Based Code Generation for Numerical Relativity... and Beyond!

To explore the NRPy+ tutorial without downloading, check out the [NRPy+ tutorial mybinder](https://mybinder.org/v2/gh/zachetienne/nrpytutorial/master?filepath=NRPyPlus_Tutorial.ipynb).

If you would like to explore the NRPy+ tutorial on your local computer, you'll need to install Python, Jupyter, Sympy, and Matplotlib. Once they are installed, [you may find this useful](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/execute.html)

In certain circumstances, developers may wish to execute one of these Jupyter notebooks from the command line. For example, when the notebook constructs an [Einstein Toolkit](https://einsteintoolkit.org) thorn. In such a case, the following command should be useful:

jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=-1 [Jupyter notebook file]
