# Common parameters for VacuumMaxwell evolutions

# Author: Terrence Pierre Jacques
#         terrencepierrej **at** gmail **dot* com

# License: BSD 2-Clause

# COMPLETE DOCUMENTATION (JUPYTER NOTEBOOKS):
# START PAGE (start here!):  ../NRPy+_Tutorial.ipynb
# THIS MODULE: ../Tutorial-VacuumMaxwell_Flat_Cartesian_ID.ipynb

# Step P1: Import needed NRPy+ core modules:
import NRPy_param_funcs as par    # NRPy+: Parameter interface

# The name of this module ("CommonParams") is given by __name__:
thismodule = __name__

# Parameters common to/needed by all VacuumMaxwell Python modules

# Step P2: Define the C parameters amp, lam, time, and wavespeed.
#          These variables proper SymPy variables, so they can be
#          used in SymPy expressions. In the C code, it acts
#          just like a usual parameter, whose value is
#          specified in the parameter file.

# amplitude
amp = par.Cparameters("REAL",thismodule,"amp", default_vals=1.0)

# lambda
lam = par.Cparameters("REAL",thismodule,"lam", default_vals=1.0)

time = par.Cparameters("REAL",thismodule,"time", default_vals=0.0)

wavespeed = par.Cparameters("REAL",thismodule,"wavespeed", default_vals=1.0)
