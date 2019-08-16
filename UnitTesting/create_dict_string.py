from mpmath import mp, mpf, mpc
from UnitTesting.standard_constants import precision

# [create_dict_string] takes in a value dictionary [value_dict] and returns a string representation of that dictionary
# that can be easily printed to the console or to a file.

# Called by calc_error, first_time_print


def create_dict_string(value_dict):

    # Setting proper precision
    mp.dps = precision

    # Initializing return_string
    return_string = ''

    # For each entry in the sorted dictionary, add properly formatted dictionary entry to return_string based on
    # type of variable (mpf, mpc, other)
    for var, num in sorted(value_dict.items(), key=lambda s: s[0].lower()):
        if type(num) == mpf:
            return_string += "'" + var + "': mpf('" + str(num) + "'), "
        elif type(num) == mpc:
            return_string += "'" + var + "': mpc(real='" + str(num.real) + "', imag='" + str(num.imag) + "'), "
        else:
            return_string += "'" + var + "': " + str(num) + ", "

    # Add dictionary brackets and remove extra ", " at the end
    return_string = '{' + return_string[0:-2] + '}'

    return return_string
