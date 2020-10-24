# indexedexpNRPyPN.py: functions related to indexed expressions,
# including e.g., tensors and pseudotensors.
# *** This is a stripped-down version of the indexedexp.py module
#     in the NRPy+ root directory.

# Step 1: Load needed modules
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import sys                       # Standard Python module for multiplatform OS-level functions

thismodule = __name__

def zerorank1(DIM=-1):
    if DIM == -1:
        DIM = 3  # default to 3D
    return [sp.sympify(0) for i in range(DIM)]

def zerorank3(DIM=-1):
    if DIM == -1:
        DIM = 3  # default to 3D
    return [[[sp.sympify(0) for i in range(DIM)] for j in range(DIM)] for k in range(DIM)]

def declarerank1(objname, DIM=-1):
    if DIM==-1:
        DIM = 3  # default to 3D
    return [sp.sympify(objname + str(i)) for i in range(DIM)]

class NonInvertibleMatrixError(ZeroDivisionError):
    """ Matrix Not Invertible; Division By Zero """

# Define the rank-3 version of the Levi-Civita symbol.
def LeviCivitaSymbol_dim3_rank3():
    LeviCivitaSymbol = zerorank3(DIM=3)

    for i in range(3):
        for j in range(3):
            for k in range(3):
                # From https://codegolf.stackexchange.com/questions/160359/levi-civita-symbol :
                LeviCivitaSymbol[i][j][k] = (i - j) * (j - k) * (k - i) * sp.Rational(1,2)
    return LeviCivitaSymbol

if __name__ == "__main__":
    import doctest
    sys.exit(doctest.testmod()[0])
