# WARNING: Importing more than the bare minimum with mpmath will result in errors on eval() below.
# This is because we need SymPy to evaluate that expression, not mpmath.
from mpmath import mp, mpf, sqrt, pi, mpc, fabs
import random
from sympy import cse, N
import UnitTesting.standard_constants as standard_constants
import logging
import hashlib

# Called by run_test
# Uses self.variable_dict

# Setting precision
precision = standard_constants.precision


# cse_simplify_and_evaluate_sympy_expressions gets self.expanded_variable_dict by calling expand_variable_dict() on
# self.variable_dict. It then gets every free symbol in the resulting expanded_variable_dict and assigns them random,
# but consistent, mpf values. It then looks at each expression in self.expanded_variable_dict, uses SymPy's CSE
# algorithm to optimize our random value substitution, and gets a mpf (or mpc if the result is complex) value for each
# expression. It then checks whether any of these expressions are super close to zero -- if a value is, recalculate it
# at twice the precision. If the value drops off in magnitude, set it to exactly zero. Finally, return a dictionary
# where each variable is assigned its respective value.
def cse_simplify_and_evaluate_sympy_expressions(self):

    # If an empty variable dict is passed, return an empty dictionary
    if self.variable_dict == {}:
        return {}

    # Call expand_variable_dict
    expanded_variable_dict = expand_variable_dict(self.variable_dict)

    # Setting precision
    mp.dps = precision

    # Creating free_symbols_set, which stores all free symbols from all expressions.
    logging.debug(' Getting all free symbols...')
    free_symbols_set = set()
    for val in expanded_variable_dict.values():
        try:
            free_symbols_set = free_symbols_set | val.free_symbols
        except AttributeError:
            pass

    # Initializing free_symbols_dict
    free_symbols_dict = dict()

    logging.debug(' ...Setting each free symbol to a random value...')

    # Setting each variable in free_symbols_set to a random number in [0, 1) according to the hashed string
    # representation of each variable.
    for var in free_symbols_set:
        # Make sure M_PI is set to its correct value, pi
        if str(var) == "M_PI":
            free_symbols_dict[var] = mpf(pi)
        # Then make sure M_SQRT1_2 is set to its correct value, 1/sqrt(2)
        elif str(var) == "M_SQRT1_2":
            free_symbols_dict[var] = mpf(1/sqrt(2))
        # All other free variables are set to random numbers
        else:
            # Take the variable [var], turn it into a string, encode the string, hash the string using the md5
            # algorithm, turn the hash into a hex number, turn the hex number into an int, set the random seed to
            # that int. This ensures each variable gets a unique but consistent value.
            random.seed(int(hashlib.md5(str(var).encode()).hexdigest(), 16))
            # Store the random value in free_symbols_dict as a mpf
            free_symbols_dict[var] = mpf(random.random())

    # Initialize calculated_dict and simplified_expression_dict
    calculated_dict = dict()
    simplified_expression_dict = dict()

    logging.debug(' ...Calculating values for each variable based on free symbols...')

    # Evaluating each expression using the values in var_dict
    for var, expression in expanded_variable_dict.items():
        # Using SymPy's cse algorithm to optimize our value substitution
        replaced, reduced = cse(expression, order='none')

        # Calculate our result_value
        result_value = calculate_value(free_symbols_dict, replaced, reduced)

        # Check if result_value is near-zero, and double checking if it should be zero
        if fabs(result_value) != mpf('0.0') and fabs(result_value) < 10 ** ((-2.0/3)*precision):
            logging.info("Found |result| (" + str(fabs(result_value)) + ") close to zero. "
                         "Checking if indeed it should be zero.")
            new_result_value = calculate_value(free_symbols_dict, replaced, reduced, precision_factor=2)
            if fabs(new_result_value) < 10 ** (-(4.0/3) * precision):
                logging.info("After re-evaluating with twice the digits of precision, |result| dropped to " +
                             str(new_result_value) + ". Setting value to zero")
                result_value = mpf('0.0')

        # Store result_value in calculated_dict
        calculated_dict[var] = result_value

    return calculated_dict


# calculate_value takes in a dictionary of free symbols [free_symbols_dict], the outputs from SymPy's cse algorithm
# [replaced] and [reduced], and an optional argument precision factor that defaults to 1. It calculates a numerical
# value for reduced[0] using value substitution between [free_symbols_dict] and [replaced] at a precision of
# [precision] * [precision_factor]
def calculate_value(free_symbols_dict, replaced, reduced, precision_factor=1):

    # We only care about the first argument of [reduced]
    reduced = reduced[0]

    # Set precision to [precision] multiplied by [precision_factor]
    mp.dps = precision_factor * precision

    # Replacing old expressions with new expressions and storing result in free_symbols_dict
    for new, old in replaced:
        keys = old.free_symbols
        for key in keys:
            old = old.subs(key, free_symbols_dict[key])
        free_symbols_dict[new] = old

    # Evaluating expression after cse optimization substitution
    keys = reduced.free_symbols
    for key in keys:
        reduced = reduced.subs(key, free_symbols_dict[key])

    # Adding our variable, value pair to our calculated_dict
    try:
        res = mpf(reduced)
    # If value is a complex number, store it as a mpc
    except TypeError:
        res = mpc(N(reduced))

    # Set the precision back
    mp.dps = precision

    return res


# expand_variable_dict takes in a variable dictionary [variable_dict] and returns a dictionary that represents
# the expanded version of [variable_dict] according to the dimension of each variable in [variable_dict].
def expand_variable_dict(variable_dict):

    # Initialize the result dictionary
    result_dict = dict()

    # Iterate through all elements of variable_dict
    for var, expression_list in variable_dict.items():

        # Getting the dimension and length of expression list
        dim, length = get_variable_dimension(expression_list)

        # If list is a scalar, easy computation with no necessary indexing
        if dim == 0:
            result_dict[var] = expression_list
        # Otherwise, need to do more work
        else:
            # Initialize our counter of the correct dimension
            counter = '0' * dim
            # Call flatten on our expression list to get a flattened list
            flattened_list = flatten(expression_list)

            # Append next element to var list and increment counter
            for elt in flattened_list:
                result_dict[form_string(var, counter)] = elt
                counter = increment_counter(counter, length)

    return result_dict


# get_variable_dimension takes in a tensor [tensor] and returns a tuple containing the rank and dimension of [tensor]
# scalar -> rank 0 tensor with length 1 -> 0, 1
# vector with 5 elements -> rank 1 tensor with length 5 -> 1, 5
# tensor with NxN elements -> rank 2 tensor with length N -> 2, N
def get_variable_dimension(tensor):
    # Assume rank is 0, dim is 1
    rank = 0
    dim = 1

    # While an element of [tensor] is also a tensor, increment dim
    while isinstance(tensor, list):
        # Set length to the number of elements of [tensor]
        if rank == 0:
            dim = len(tensor)
        rank += 1
        tensor = tensor[0]

    return rank, dim


# flatten takes in a list [l] and an optional argument [fl] that represents the flattened list.
# It returns [fl], where [fl] is [l] but with all nested lists removed.
# Adapted from https://www.geeksforgeeks.org/python-convert-a-nested-list-into-a-flat-list/
def flatten(l, fl=None):

    if fl is None:
        fl = []

    for i in l:
        if type(i) == list:
            flatten(i, fl)
        else:
            fl.append(i)
    return fl


# form_string takes in a variable [var] and a counter [counter] and returns a string representing [var] at
# index [counter]
def form_string(var, counter):
    return_string = var
    for char in counter:
        return_string += '[' + char + ']'

    return return_string


# increment_counter takes in a counter [counter] and a length [length] and returns the next number after counter
# in base [length].
# Example: iter_counter('00', 2) -> '01'
#          iter_counter('02', 3) -> '10'
#          iter_counter('01111', 2) -> '10000'
def increment_counter(counter, length):

    # Set return_string to empty string, set num to 1
    return_string = ''
    num = 1

    # Loop backwards through each character [char] in [counter]
    for char in reversed(counter):

        # Add [num] to the integer representation of [char]
        digit = int(char) + num
        # If it's time to loop back around
        if digit == length:
            # Add a 0 to the return string, num = 1
            return_string += '0'
            num = 1
            # Add current digit to the return string, num = 0
        else:
            return_string += str(digit)
            num = 0

    # Return reversed return_string since we built it backwards
    return return_string[::-1]
