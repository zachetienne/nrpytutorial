import logging
from mpmath import log10, fabs, mp
from datetime import date
from UnitTesting.standard_constants import precision
from UnitTesting.create_dict_string import create_dict_string


# calc_error loops through each item in self.calculated_dict and self.trusted_values_dict_entry,
# and makes sure that the difference between each respective value's number of decimal in the dictionaries places is
# less than 1/2 of the precision. It returns a boolean representing whether any variables differed.

# Called by run_test

# Uses self.calculated_dict, self.trusted_values_dict_entry, self.module_name


def calc_error(self):
    # Setting precision
    mp.dps = precision

    # Creating sets to easily compare the keys of calculated_dict and trusted_dict
    calculated_set = set(self.calculated_dict)
    trusted_set = set(self.trusted_values_dict_entry)

    logging.debug(' Checking that calculated and trusted dicts contain the same variables...')
    # If the sets differ, output the differing variables
    if calculated_set != trusted_set:
        logging.error(' {}: Calculated dictionary and trusted dictionary have different variables.'
                      .format(self.module_name))
        calculated_minus_trusted = calculated_set - trusted_set
        trusted_minus_calculated = trusted_set - calculated_set
        if calculated_minus_trusted != set([]):
            logging.error(' Calculated Dictionary variables not in Trusted Dictionary: ' +
                          str(sorted(calculated_minus_trusted)))
        if trusted_minus_calculated != set([]):
            logging.error(' Trusted Dictionary variables not in Calculated Dictionary: ' +
                          str(sorted(trusted_minus_calculated)))
        # Automatically fail and don't proceed
        return False

    logging.debug(' ...Success: same variables in both dicts.\n')

    # Initialize list of variables whose values differ
    bad_var_list = []

    logging.debug(' Comparing all calculated and trusted values...')
    # Loop through all variables in sorted order
    for var in sorted(self.calculated_dict):

        # Values to compare
        calculated_val = self.calculated_dict[var]
        trusted_val = self.trusted_values_dict_entry[var]

        # For each variable, print calculated and trusted values
        logging.debug('\n' + self.module_name + ': ' + var + ': Calculated: ' + str(
            calculated_val) + '\n' + self.module_name + ': ' + var
                      + ': Trusted:    ' + str(trusted_val) + '\n')

        # Calculate the error between both values
        if trusted_val == 0:
            log10_relative_error = log10(fabs(calculated_val))
        elif calculated_val == 0:
            log10_relative_error = log10(fabs(trusted_val))
        else:
            log10_relative_error = log10(fabs((trusted_val - calculated_val) / trusted_val))

        # Boolean determining if their difference is within the tolerance we accept
        good = (log10_relative_error < (precision / -2.0))

        # Store all variables who are not 'good'
        if not good:
            bad_var_list.append(var)

    # If we want to output and there exists at least one variable with error, print
    if bad_var_list:
        logging.error('''
\nVariable(s) {} in module {} failed. Please check values.
If you are confident that the newly calculated values are correct, comment out the old trusted values for
{} in your trusted_values_dict and copy the following code between the ##### into your trusted_values_dict.
Make sure to fill out the TODO comment describing why the values had to be changed. Then re-run test script.

#####

# Generated on: {}
# Reason for changing values: TODO
trusted_values_dict['{}'] = {}

#####
'''.format(bad_var_list, self.module_name, self.trusted_values_dict_name, date.today(),
           self.trusted_values_dict_name, create_dict_string(self.calculated_dict)))
    else:
        logging.debug(' ...Success: all variables identical.\n')

    # Return True if all variables are good, False otherwise
    return bad_var_list == []
