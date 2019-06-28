import logging

from UnitTesting.calc_error import calc_error
from UnitTesting.first_time_print import first_time_print
from UnitTesting.evaluate_globals import evaluate_globals
from UnitTesting.expand_variable_dict import expand_variable_dict
from UnitTesting.var_dict_to_value_dict import var_dict_to_value_dict
from UnitTesting.is_first_time import is_first_time
from UnitTesting.create_trusted_globals_dict import create_trusted_globals_dict
from UnitTesting.standard_constants import precision
from time import time
from mpmath import mp


# run_test takes in :
# [self]- The unittest self object,
# [mod_dict]- The user-supplied dictionary of modules
# [trusted_values_dict]- The dictionary of trusted values
# [locs]- The current local variables in the workspace. Should ALWAYS be locals()
# It then runs a unittest, comparing calculated values with trusted values.
# Throws an [AssertionError] if [mod_dict] is empty
def run_test(self, mod_dict, trusted_values_dict, locs):

    # Can't use empty dictionaries
    assert mod_dict != dict()

    mp.dps = precision

    # Determining if this is the first time the code is run based of the existence of trusted values
    first_times = is_first_time(mod_dict, trusted_values_dict)

    # Creating trusted dictionary based off names of modules in ModDict
    trusted_dict = create_trusted_globals_dict(mod_dict, trusted_values_dict, first_times)

    # Timing how long evaluate_globals takes
    t = time()

    # Creating dictionary of expressions for all modules in ModDict
    result_dict = evaluate_globals(mod_dict, locs)

    # Printing the time it took to run evaluate_globals
    logging.debug(str(time()-t) + ' seconds to run evaluate_globals')

    del mod_dict

    # # If it is the first time for at least one module, sort the module dictionary based on first_times.
    # # This makes it so the new modules are done last. This makes it easy to copy the necessary modules' code.
    # if True in first_times:
    #     # https://stackoverflow.com/questions/13668393/python-sorting-two-lists
    #     first_times, result_mods = (list(x) for x in zip(*sorted(zip(first_times, result_dict))))
    #
    #     temp_dict = dict()
    #
    #     # Creates dictionary based on order of first_times
    #     for mod in result_mods:
    #         temp_dict[mod] = result_dict[mod]
    #
    #     # Updates resultDict to be in this new order
    #     result_dict = temp_dict
    #     del temp_dict, result_mods

    # Looping through each module in resultDict
    for mod in result_dict:

        var_dict = result_dict[mod]
        first_time = first_times[mod]

        if not first_time:
            logging.info('Currently working on module ' + mod + '...')

        # Generating variable list and name list for module
        new_dict = expand_variable_dict(var_dict)

        # Timing how long list_to_value_list takes
        t = time()

        # Calculating numerical list for module
        value_dict = var_dict_to_value_dict(new_dict)

        # Printing the time it took to run list_to_value_list
        logging.debug(str(time()-t) + ' seconds to run list_to_value_list')

        # If being run for the first time, print the code that must be copied into trusted_values_dict
        if first_time:
            first_time_print(mod, value_dict)

        # Otherwise, compare calculated values to trusted values
        else:

            # Calculates the error between mod_dict and trusted_dict[mod] for the current module
            values_identical = calc_error(mod, value_dict, trusted_dict[mod])

            # If at least one value differs, print exit message and fail the unittest
            if not values_identical:
                self.assertTrue(values_identical,
                                'Variable above has different calculated and trusted values. Follow '
                                'above instructions.')

            # If every value is the same, completed module.
            else:
                logging.info('Completed module ' + mod + ' with no errors.\n')
            self.assertTrue(values_identical)

    # If it's the first time for at least one module
    if True in first_times:
        self.assertTrue(False, 'Automatically failing due to first time for at least one module. Please see above'
                               'for the code to copy into your trusted_values_dict.')
