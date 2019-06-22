import logging
from UnitTesting.trusted_values_dict import trusted_values_dict
from mpmath import log10,fabs, mp
from datetime import date

# Takes in a module [mod], a calculated dictionary [calculated_dict], a trusted dictionary [trusted_dict], and
# a symbolic dictionary [symbolic_dict] and computes the error for each result-trusted pair for each respective index.
# Logs debug statements for each pair of values if the logging level is <= DEBUG
# and logs a failure message if logging level is <= ERROR.
# Returns a boolean [good] that represents if any two value pairs didn't differ
# by more than (precision/2) decimal places.

# Called by run_test


def calc_error(mod, calculated_dict, trusted_dict, output=True):

    # Precision for the module based off the set precision in trusted_values_dict
    precision = trusted_values_dict['precision']
    mp.dps = precision

    # Creating sets to easily compare the keys of resultDict and trustedDict
    calculated_set = set(calculated_dict)
    trusted_set = set(trusted_dict)

    # If resultDict and trustedDict have different variables, print the differing variables
    if calculated_set != trusted_set:
        if output:
            logging.error('\n\t' + mod + ': Calculated dictionary and trusted dictionary have different variables.')
            calculated_minus_trusted = calculated_set - trusted_set
            trusted_minus_calculated = trusted_set - calculated_set
            if calculated_minus_trusted != set([]):
                logging.error('\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' +
                              str(sorted(calculated_minus_trusted)))
            if trusted_minus_calculated != set([]):
                logging.error('\n\tTrusted Dictionary variables not in Calculated Dictionary: \n\t' +
                              str(sorted(trusted_minus_calculated)))
        return False

    del calculated_set, trusted_set

    # For each variable, print calculated and trusted values
    for var in sorted(calculated_dict):
        calculated_num = calculated_dict[var]
        trusted_num = trusted_dict[var]

        if output:
            logging.debug('\n' + mod + ': ' + var + ': Calculated: ' + str(calculated_num) + '\n' + mod + ': ' + var
                          + ': Trusted:    ' + str(trusted_num) + '\n')

        if trusted_num == 0:
            log10_relative_error = log10(fabs(calculated_num))
        elif calculated_num == 0:
            log10_relative_error = log10(fabs(trusted_num))
        else:
            log10_relative_error = log10(fabs((trusted_num - calculated_num) / trusted_num))

        good = (log10_relative_error < (precision / -2))
        if not good:
            if output:
                logging.info('\n\nVariable ' + "'" + var + "'" + ' in module ' + str(mod) + ' failed. Please check' +
                             ' values.\n\n' + 'If you are confident that the newly calculated values are correct, ' +
                             'comment out the old trusted values for ' + "'" + mod + "Globals'" +
                             ' in trusted_values_dict and copy the following code between the ##### into ' +
                             'trusted_values_dict. Make sure to fill out the TODO comment describing why the values' +
                             ' had to be changed. Then re-run test script.\n' + '#####\n\n# Generated on: ' +
                             str(date.today()) + '\n# Reason for changing values: TODO' + "\ntrusted_values_dict['" +
                             mod + "Globals'] = " + create_dict_string(calculated_dict) + '\n\n#####')
            return False

    return True


# Subfunction to properly format dict to print
def create_dict_string(calculated_dict):

    return_string = '{'

    for var, num in calculated_dict.items():
        return_string += "'" + var + "': mpf('" + str(num) + "'), "

    return_string = return_string[0:-2] + '}'

    return return_string
