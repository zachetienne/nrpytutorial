# NRPyPN: Validated Post-Newtonian Expressions for Input into TwoPunctures, Wolfram Mathematica, SymPy, or Highly-Optimized C Codes

## Author
Zachariah B. Etienne (http://astro.phys.wvu.edu/zetienne/)

## Special Thanks
Special thanks to Peter Diener and Roland Haas for reviewing
NRPyPN in preparation for its inclusion into the Einstein
Toolkit. Also special thanks to Antoni Ramos-Buades for sharing
the Mathematica notebooks that the paper
"Simple procedures to reduce eccentricity of binary black hole simulations"
Phys. Rev. D 99, 023003 (2019)
used, so that expressions in NRPyPN could be validated.

## Purpose

NRPyPN primarily focuses on implementation and validation of
post-Newtonian expressions, with the immediate goal of generating
high-PN-order tangential and radial momenta for binary black hole
initial data with minimal eccentricity. These momenta can be
directly injected into e.g., TwoPunctures to set up quasicircular
binary black hole initial data.

NRPyPN bases its approach on
"Simple procedures to reduce eccentricity of binary black hole simulations",
Ramos-Buades, Husa, and Pratten,
https://arxiv.org/abs/1810.00036,
Phys. Rev. D 99, 023003 (2019)

and

"Post-Newtonian Quasicircular Initial Orbits for Numerical Relativity",
Healy, Lousto, Nakano, and Zlochower,
https://arxiv.org/abs/1702.00872,
Class. Quant. Grav. 34 (2017) 14, 145011 


## Installation instructions

Prerequisites:

* Python 3.6+ preferred, though earlier versions are supported
* pip, the Python package manager, which should come with Python.
* (*Optional*) Pandoc (https://pandoc.org/), to enable PDF conversion of NRPyPN notebooks

Python packages:

* SymPy 1.2+
* Jupyter

### Quick install from the command line (bash shell)

* First set up a virtual environment:

python3 -m venv nrpyvirtualenv
source nrpyvirtualenv/bin/activate
pip install -U sympy jupyter

* Next navigate to Cactus/arrangements/EinsteinInitialData/NRPyPN, and run:

jupyter notebook

The NRPyPN.ipynb notebook both contains the Table of Contents and provides
a simple interface for generating quasicircular 


## Using NRPyPN in the Einstein Toolkit, with TwoPunctures

1. Follow the above installation instructions, launch Jupyter, then open NRPyPN.ipynb
2. Scroll down to the bottom of NRPyPN.ipynb and insert the desired black hole binary parameters
3. Click the "fast-forward" button at the top of the Jupyter notebook, then click "Restart and run all cells"
4. The PN tangential and radial momenta, for insertion into TwoPunctures, will be output at the bottom.

## License:
BSD 2-Clause

Copyright (c) 2020, Zachariah Etienne
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## Required Citation

1) Bibtex entry:

@misc{NRPyPN,
 author       = {Etienne, Zachariah B.},
 title        = {NRPyPN: Validated Post-Newtonian Expressions for Input into Wolfram Mathematica, SymPy, or Highly Optimized C Codes},
 month        = nov,
 year         = 2020,
 url          = {https://github.com/zachetienne/nrpytutorial/blob/master/NRPyPN/}
}

## Suggested Citation

1) Bibtex entry:

@article{Habib:2020dba,
    author = "Habib, Sarah and Ramos-Buades, Antoni and Huerta, E.A. and Husa, Sascha and Haas, Roland and Etienne, Zachariah",
    title = "{Initial Data and Eccentricity Reduction Toolkit for Binary Black Hole Numerical Relativity Waveforms}",
    eprint = "2011.08878",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "11",
    year = "2020"
}
