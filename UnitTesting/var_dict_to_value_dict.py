# WARNING: Importing more than the bare minimum with mpmath will result in errors on eval() below.
# This is because we need SymPy to evaluate that expression, not mpmath.
from mpmath import mp, mpf, sqrt, pi, mpc
from random import seed, random
from sympy import cse, simplify
import UnitTesting.standard_constants as standard_constants


# Takes in a variable dictionary [var_dict] a precision value [precision] and a seed value [seed_value], and returns
# a dictionary with each expression in [var_dict] evaluated according to parameters (seed, precision).
# Throws an [AttributeError] if the variable list being passed in has no sympy symbols

# Called by run_test

def var_dict_to_value_dict(var_dict):

    if var_dict == dict():
        return dict()

    # Setting precision
    mp.dps = standard_constants.precision

    # List all the free symbols in the expressions in [var_dict].
    free_symbols_list = list(sum(var_dict.values()).free_symbols)

    # Sort free symbols list based off the alphanumeric strings for each variable
    free_symbols_list.sort(key=lambda v: str(v))

    # Set the random seed according to seed in trustedValuesDict:
    seed(standard_constants.seed)

    # Creating dictionary entry for each variable and its pseudorandom value in [0,1) as determined by seed
    variable_dictionary = dict()

    for var in free_symbols_list:
        if str(var) == "M_PI":
            variable_dictionary[var] = mpf(pi)
        # Then make sure M_SQRT1_2 is set to its correct value, 1/sqrt(2), to the desired number of significant digits:
        elif str(var) == "M_SQRT1_2":
            variable_dictionary[var] = mpf(1/sqrt(2))
        # All other free variables are set to random numbers
        else:
            variable_dictionary[var] = sqrt(mpf(random()))

    value_dict = dict()
    # Evaluating each expression using the values in variable_dictionary
    for var, expression in var_dict.items():
        # Copying variable_dictionary into a new variable dictionary
        new_var_dict = dict(variable_dictionary)

        # Using sympy's cse algorithm to optimize our value substitution
        replace, reduce = cse(expression, order='none')
        reduce = reduce[0]

        # Replacing old expressions with new expressions and putting result in new variable dictionary
        for new, old in replace:
            keys = old.free_symbols
            for key in keys:
                old = old.subs(key, new_var_dict[key])
            new_var_dict[new] = old

        # Evaluating expression after cse optimization
        keys = reduce.free_symbols
        for key in keys:
            reduce = reduce.subs(key, new_var_dict[key])

        # Adding our variable, value pair to our value_dict
        try:
            value_dict[var] = mpf(reduce)
        # If value is a complex number, store it as a mpc
        except TypeError:
            reduce = simplify(reduce)
            value_dict[var] = mpc(reduce)

    return value_dict
