# Step 1.a: Initialize core Python/UnitTesting modules
from UnitTesting.calc_error import calc_error
from UnitTesting.evaluate_globals import evaluate_globals
from UnitTesting.first_time_print import first_time_print
from UnitTesting.cse_simplify_and_evaluate_sympy_expressions import cse_simplify_and_evaluate_sympy_expressions
from UnitTesting.standard_constants import precision
from mpmath import mp
from importlib import import_module
import logging


def run_test(self):

    # Step 1: Setup

    logging.info(' Currently working on function ' + self.function + ' in module ' + self.module_name + '...\n')

    # Step 1.a: Set precision to the value defined in standard_constants
    mp.dps = precision

    # Step 1.b: Import trusted_values_dict from trusted_values_dict.py in self.path
    logging.info(' Importing trusted_values_dict...')
    self.trusted_values_dict = import_module('trusted_values_dict').trusted_values_dict
    logging.info(' ...Success: Imported trusted_values_dict.\n')

    # Step 1.c: Set boolean self.first_time based on existence of desired trusted_values_dict entry
    self.first_time = self.trusted_values_dict_name not in self.trusted_values_dict
    if self.first_time:
        logging.info(' Proper entry not in trusted_values_dict -- '
                     'this function in this module is being run for the first time.')

    # Step 1.e: Set trusted_values_dict_entry to its corresponding trusted_values_dict entry
    self.trusted_values_dict_entry = {} if self.first_time else self.trusted_values_dict[self.trusted_values_dict_name]

    # Step 2: Calculation

    # Step 2.a: Call evaluate_globals which calls self.function and gets expressions for all globals in self.global_list
    logging.info(' Calling evaluate_globals...')
    self.variable_dict = evaluate_globals(self)
    logging.info(' ...Success: evaluate_globals ran without errors.\n')

    # Step 2.b: Call cse_simplify_and_evaluate_sympy_expressions to assign each variable in each expression a random
    #           value and calculate the numerical result
    logging.info(' Calling cse_simplify_and_evaluate_sympy_expressions...')
    self.calculated_dict = cse_simplify_and_evaluate_sympy_expressions(self)
    logging.info(' ...Success: cse_simplify_and_evaluate_sympy_expressions ran without errors.\n')

    # Step 3: Comparison

    if self.first_time:
        # Step 3.a: Print self.calculated_dict in a nice format and append it to trusted_values_dict
        logging.info(' Calling first_time_print since it is being run for the first time...')
        first_time_print(self)
        logging.info(' ...Success: first_time_print ran without errors. Automatically failing due to first_time.\n')
        self.assertTrue(False)

    else:
        # Step 3.b: Call calc_error to calculate the error between the trusted values and the calculated values
        logging.info(' Calling calc_error...')
        values_identical = calc_error(self)

        # If there is an error large enough, fail
        if not values_identical:
            self.assertTrue(values_identical,
                            'Variable(s) above have different calculated and trusted values. Follow '
                            'instructions above.')
        # Otherwise, pass
        else:
            logging.info(' ...Success: calc_error ran without errors.\n')
