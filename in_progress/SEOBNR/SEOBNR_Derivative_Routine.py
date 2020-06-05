# As documented in the NRPy+ tutorial module
#   SEOBNR_Derivative_Routine.ipynb,
#   this module computes partial derivatives of
#   an input list of expressions with respect
#   to an input list of free variables.

# Authors: Zachariah B. Etienne & Tyler Knowles
#         zachetie **at** gmail **dot* com

# Step 1.a: import all needed modules from Python/NRPy+:
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import os, sys                    # Standard Python modules for multiplatform OS-level functions

# Step 1.?: check system path so can use outputC; #TylerK: remove and put outputC back with other imports
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

from outputC import *             # TylerK: check what is imported and remove *; also find appropriate description
