from get_variable_dimension import get_variable_dimension

# expand_variable_dict takes in a variable dictionary [variable_dict] and returns a dictionary that represents
# the expanded version of [variable_dict] according to the dimension of each tensor.

# Called by run_test


def expand_variable_dict(variable_dict):

    result_dict = dict()

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
            flattened_list = flatten(expression_list, [])

            for elt in flattened_list:
                # Append element to var list
                result_dict[form_string(var, counter)] = elt
                counter = increment_counter(counter, length)

    return result_dict

# Subfunctions:


# iter_counter takes in a counter [counter] and a length [length] and returns the next number after counter
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


# Used to form the proper string to be added to name list based on var and counter
def form_string(var, counter):
    return_string = var
    for char in counter:
        return_string += '[' + char + ']'

    return return_string


# Function used for removing nested lists in python.
# https://www.geeksforgeeks.org/python-convert-a-nested-list-into-a-flat-list/
def flatten(l, fl):
    for i in l:
        if type(i) == list:
            flatten(i, fl)
        else:
            fl.append(i)
    return fl



