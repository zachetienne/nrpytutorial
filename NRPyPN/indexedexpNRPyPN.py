# indexedexpNRPyPN.py: functions related to indexed expressions,
# including e.g., tensors and pseudotensors.
# *** This is a stripped-down version of the indexedexp.py module
#     in the NRPy+ root directory.

# Step 1: Load needed modules
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import sys                       # Standard Python module for multiplatform OS-level functions
import re                        # Standard Python module for regular expressions

thismodule = __name__

def declare_indexedexp(rank, symbol=None, dimension=None):
    """ Generate an indexed expression of specified rank and dimension
    """
    if not dimension or dimension == -1:
        dimension = 3  # default to 3 dimension
    if symbol is not None:
        if not isinstance(symbol, str) or not re.match(r'[\w_]', symbol):
            raise ValueError('symbol must be an alphabetic string')
    if dimension is not None:
        if not isinstance(dimension, int) or not dimension > 0:
            raise ValueError('dimension must be a positive integer')
    loop_index = ['str(%s)' % chr(97 + n) for n in range(rank)]
    indexing = ' + '.join(loop_index)
    interior = 'sympify(0)' if not symbol \
          else 'sympify(\'%s\' + %s)' % (symbol, indexing)
    indexedexp = '[' * rank + interior
    for i in range(1, rank + 1):
        indexedexp += ' for %s in range(%s)]' % (loop_index[rank - i][4:-1], dimension)
    return eval(indexedexp, {'sympify': sp.sympify}, {})

def zerorank1(DIM=-1):
    return declare_indexedexp(rank=1, dimension=DIM)

def declarerank1(symbol, DIM=-1):
    return declare_indexedexp(rank=1, symbol=symbol, dimension=DIM)

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
