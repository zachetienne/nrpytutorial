# Common parameters for scalar wave evolutions
#
# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
#
# License: BSD 2-Clause

# COMPLETE DOCUMENTATION (JUPYTER NOTEBOOKS):
# START PAGE (start here!):  ../NRPy+_Tutorial.ipynb
# THIS MODULE: ../Tutorial-Scalarwave.ipynb

# Step P1: Import needed NRPy+ core modules:
import NRPy_param_funcs as par    # NRPy+: Parameter interface

thismodule = __name__

# Parameters common to/needed by all ScalarWave Python modules

# Step P2: Define the C parameter wavespeed. The `wavespeed`
#          variable is a proper SymPy variable, so it can be
#          used in below expressions. In the C code, it acts
#          just like a usual parameter, whose value is
#          specified in the parameter file.
wavespeed = par.Cparameters("REAL", thismodule, "wavespeed",1.0)
