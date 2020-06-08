from outputC import nrpyAbs      # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: parameter interface
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends

thismodule = __name__

TINYDOUBLE = par.Cparameters("REAL", thismodule, "TINYDOUBLE", 1e-100)

def min_noif(a,b):
    # Returns the minimum of a and b
    if a==sp.sympify(0):
        return sp.Rational(1,2) * (b-nrpyAbs(b))
    if b==sp.sympify(0):
        return sp.Rational(1,2) * (a-nrpyAbs(a))
    return sp.Rational(1,2) * (a+b-nrpyAbs(a-b))

def max_noif(a,b):
    # Returns the maximum of a and b
    if a==sp.sympify(0):
        return sp.Rational(1,2) * (b+nrpyAbs(b))
    if b==sp.sympify(0):
        return sp.Rational(1,2) * (a+nrpyAbs(a))
    return sp.Rational(1,2) * (a+b+nrpyAbs(a-b))

def coord_leq_bound(x,xstar):
    # Returns 1.0 if x <= xstar, 0.0 otherwise.
    # Requires appropriately defined TINYDOUBLE
    return min_noif(x-xstar-TINYDOUBLE,0.0)/(x-xstar-TINYDOUBLE)

def coord_geq_bound(x,xstar):
    # Returns 1.0 if x >= xstar, 0.0 otherwise.
    # Requires appropriately defined TINYDOUBLE
    return max_noif(x-xstar+TINYDOUBLE,0.0)/(x-xstar+TINYDOUBLE)

def coord_less_bound(x,xstar):
    # Returns 1.0 if x < xstar, 0.0 otherwise.
    # Requires appropriately defined TINYDOUBLE
    return min_noif(x-xstar,0.0)/(x-xstar-TINYDOUBLE)

def coord_greater_bound(x,xstar):
    # Returns 1.0 if x > xstar, 0.0 otherwise.
    # Requires appropriately defined TINYDOUBLE
    return max_noif(x-xstar,0.0)/(x-xstar+TINYDOUBLE)
