from outputC import *            # NRPy+: Core C code output module

thismodule = __name__

SMALLDOUBLE = par.Cparameters("#define", thismodule, "SMALLDOUBLE", 1e-10)
# This value must be set such that SMALLDOUBLE is less than the smallest grid resolution 
# and greater than abs((machine epsilon) * (the largest grid extent))

def min_noif(a,b):
    # Returns the minimum of a and b
    return sp.Rational(1,2) * (a+b-nrpyAbs(a-b))

def max_noif(a,b):
    # Returns the maximum of a and b
    return sp.Rational(1,2) * (a+b+nrpyAbs(a-b))

def coord_leq_bound(x,xstar): 
    # Returns 1.0 if x <= xstar, 0.0 otherwise. 
    # Requires appropriately defined SMALLDOUBLE
    return min_noif(x-xstar-SMALLDOUBLE,0.0)/min_noif(x-xstar-SMALLDOUBLE,-SMALLDOUBLE)

def coord_geq_bound(x,xstar): 
    # Returns 1.0 if x >= xstar, 0.0 otherwise. 
    # Requires appropriately defined SMALLDOUBLE
    return max_noif(x-xstar+SMALLDOUBLE,0.0)/max_noif(x-xstar+SMALLDOUBLE,SMALLDOUBLE)

def coord_less_bound(x,xstar): 
    # Returns 1.0 if x < xstar, 0.0 otherwise. 
    # Requires appropriately defined SMALLDOUBLE
    return min_noif(x-xstar,0.0)/min_noif(x-xstar,-SMALLDOUBLE)

def coord_greater_bound(x,xstar): 
    # Returns 1.0 if x > xstar, 0.0 otherwise. 
    # Requires appropriately defined SMALLDOUBLE
    return max_noif(x-xstar,0.0)/max_noif(x-xstar,SMALLDOUBLE)
