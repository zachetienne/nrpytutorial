from datetime import date
from UnitTesting.create_dict_string import create_dict_string
import logging
import os

# [first_time_print] takes in a module [mod], a value dictionary [value_dict], a path [path], and a boolean [write].
# It prints to the console the properly formatted trusted_values_dict entry based on [mod] and [value_dict].
# Additionally, if [write] is [True], it appends this output to the file [path]/trusted_values_dict.py

# Called by run_test

# Uses self.module_name, self.trusted_values_dict_name, self.calculated_dict, self.path


def first_time_print(self, write=True):
    logging.error('''
Module: {}
Please copy the following code between the ##### and paste it into your trusted_values_dict.py file for this module:

#####

# Generated on: {}
trusted_values_dict['{}'] = {}

#####
'''.format(self.module_name, date.today(), self.trusted_values_dict_name, create_dict_string(self.calculated_dict)))

    # If [write] is [True], write to [trusted_values_dict]
    if write:
        logging.debug(' Writing trusted_values_dict entry to trusted_values_dict.py...')
        fw = open(os.path.join(self.path, 'trusted_values_dict.py'), 'a')
        fw.write('''
# Generated on: {}
trusted_values_dict['{}'] = {}
'''.format(date.today(), self.trusted_values_dict_name, self.calculated_dict))
        fw.close()
        logging.debug(' ...Success: entry written to trusted_values_dict.py\n')
