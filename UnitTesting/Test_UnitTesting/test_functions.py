import unittest
import logging
from sys import version_info

logging.basicConfig(level='INFO')

from mpmath import mpf, mpc, mp, pi, sqrt
from UnitTesting.standard_constants import precision
mp.dps = precision


class TestFunctions(unittest.TestCase):

    # The following extends the __init__() function in Python's
    #   unittest. For more, read this:
    # https://techoverflow.net/2020/04/21/how-to-fix-python-unittest-__init__-takes-1-positional-argument-but-2-were-given/
    def setUp(self):
        self.calculated_dict = {}
        self.function = ''
        self.global_list = []
        self.initialization_string = ''
        self.module = ''
        self.module_name = ''
        self.path = ''
        self.trusted_values_dict_entry = {}
        self.trusted_values_dict_name = ''
        self.variable_dict = {}
        print('')

    def test_calc_error(self):

        from UnitTesting.calc_error import calc_error
        from datetime import date

        mp.dps = precision

        self.module_name = 'TestModule'
        self.trusted_values_dict_name = 'TestModule__globals'

        self.calculated_dict = {}
        self.trusted_values_dict_entry = {}
        self.assertEqual(True, calc_error(self))

        self.calculated_dict = {'a': mpf('1.0')}
        self.trusted_values_dict_entry = {}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'ERROR:root: TestModule: Calculated dictionary and trusted dictionary have different variables.',
                   'ERROR:root: Calculated Dictionary variables not in Trusted Dictionary: ' + "['a']"]
        calc_error_helper(self, message, False)

        self.calculated_dict = {'a': mpf('1.0'), 'b': mpf('2.0')}
        self.trusted_values_dict_entry = {}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'ERROR:root: TestModule: Calculated dictionary and trusted dictionary have different variables.',
                   'ERROR:root: Calculated Dictionary variables not in Trusted Dictionary: ' + "['a', 'b']"]
        calc_error_helper(self, message, False)

        self.calculated_dict = {'a': mpf('1.0')}
        self.trusted_values_dict_entry = {'b': mpf('2.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'ERROR:root: TestModule: Calculated dictionary and trusted dictionary have different variables.',
                   'ERROR:root: Calculated Dictionary variables not in Trusted Dictionary: ' + "['a']",
                   'ERROR:root: Trusted Dictionary variables not in Calculated Dictionary: ' + "['b']"]
        calc_error_helper(self, message, False)

        self.calculated_dict = {'a': mpf('2.0'), 'b': mpf('3.0')}
        self.trusted_values_dict_entry = {'c': mpf('1.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'ERROR:root: TestModule: Calculated dictionary and trusted dictionary have different variables.',
                   'ERROR:root: Calculated Dictionary variables not in Trusted Dictionary: ' + "['a', 'b']",
                   'ERROR:root: Trusted Dictionary variables not in Calculated Dictionary: ' + "['c']"
                   ]
        calc_error_helper(self, message, False)

        self.calculated_dict = {'a': mpf('2.0'), 'b': mpf('3.0')}
        self.trusted_values_dict_entry = {'a': mpf('1.0'), 'c': mpf('4.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'ERROR:root: TestModule: Calculated dictionary and trusted dictionary have different variables.',
                   'ERROR:root: Calculated Dictionary variables not in Trusted Dictionary: ' + "['b']",
                   'ERROR:root: Trusted Dictionary variables not in Calculated Dictionary: ' + "['c']"
                   ]
        calc_error_helper(self, message, False)

        self.calculated_dict = {'a': mpf('1.0'), 'b': mpf('2.0')}
        self.trusted_values_dict_entry = {'c': mpf('3.0'), 'd': mpf('4.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'ERROR:root: TestModule: Calculated dictionary and trusted dictionary have different variables.',
                   'ERROR:root: Calculated Dictionary variables not in Trusted Dictionary: ' + "['a', 'b']",
                   'ERROR:root: Trusted Dictionary variables not in Calculated Dictionary: ' + "['c', 'd']"
                   ]
        calc_error_helper(self, message, False)

        self.calculated_dict = {'b': mpf('1.0'), 'a': mpf('2.0')}
        self.trusted_values_dict_entry = {'c': mpf('3.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'ERROR:root: TestModule: Calculated dictionary and trusted dictionary have different variables.',
                   'ERROR:root: Calculated Dictionary variables not in Trusted Dictionary: ' + "['a', 'b']",
                   'ERROR:root: Trusted Dictionary variables not in Calculated Dictionary: ' + "['c']"
                   ]
        calc_error_helper(self, message, False)

        self.calculated_dict = {'a': mpf('1.0')}
        self.trusted_values_dict_entry = {'a': mpf('1.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'DEBUG:root: ...Success: same variables in both dicts.\n',
                   'DEBUG:root: Comparing all calculated and trusted values...',
                   'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a' + ': Trusted:    ' + '1.0' + '\n',
                   'DEBUG:root: ...Success: all variables identical.\n'
                   ]
        calc_error_helper(self, message, True)

        self.calculated_dict = {'a': mpf('1.0'), 'b': mpf('2.0')}
        self.trusted_values_dict_entry = {'a': mpf('1.0'), 'b': mpf('2.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'DEBUG:root: ...Success: same variables in both dicts.\n',
                   'DEBUG:root: Comparing all calculated and trusted values...',
                   'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a' + ': Trusted:    ' + '1.0' + '\n',
                   'DEBUG:root:' + '\nTestModule: ' + 'b' + ': Calculated: ' + '2.0' + '\nTestModule: ' + 'b' + ': Trusted:    ' + '2.0' + '\n',
                   'DEBUG:root: ...Success: all variables identical.\n'
                   ]
        calc_error_helper(self, message, True)

        self.calculated_dict = {'a': mpf('1.0'), 'b': mpf('2.0')}
        self.trusted_values_dict_entry = {'b': mpf('2.0'), 'a': mpf('1.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'DEBUG:root: ...Success: same variables in both dicts.\n',
                   'DEBUG:root: Comparing all calculated and trusted values...',
                   'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a' + ': Trusted:    ' + '1.0' + '\n',
                   'DEBUG:root:' + '\nTestModule: ' + 'b' + ': Calculated: ' + '2.0' + '\nTestModule: ' + 'b' + ': Trusted:    ' + '2.0' + '\n',
                   'DEBUG:root: ...Success: all variables identical.\n'
                   ]
        calc_error_helper(self, message, True)

        self.calculated_dict = {'b': mpf('2.0'), 'a': mpf('1.0')}
        self.trusted_values_dict_entry = {'b': mpf('2.0'), 'a': mpf('1.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'DEBUG:root: ...Success: same variables in both dicts.\n',
                   'DEBUG:root: Comparing all calculated and trusted values...',
                   'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a' + ': Trusted:    ' + '1.0' + '\n',
                   'DEBUG:root:' + '\nTestModule: ' + 'b' + ': Calculated: ' + '2.0' + '\nTestModule: ' + 'b' + ': Trusted:    ' + '2.0' + '\n',
                   'DEBUG:root: ...Success: all variables identical.\n'
                   ]
        calc_error_helper(self, message, True)

        self.calculated_dict = {'a': mpf('1.0')}
        self.trusted_values_dict_entry = {'a': mpf('2.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'DEBUG:root: ...Success: same variables in both dicts.\n',
                   'DEBUG:root: Comparing all calculated and trusted values...',
                   'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a' + ': Trusted:    ' + '2.0' + '\n',
                   'ERROR:root:' +
                   '''
\nVariable(s) {} in module TestModule failed. Please check values.
If you are confident that the newly calculated values are correct, comment out the old trusted values for
TestModule__globals in your trusted_values_dict and copy the following code between the ##### into your trusted_values_dict.
Make sure to fill out the TODO comment describing why the values had to be changed. Then re-run test script.

#####

# Generated on: {}
# Reason for changing values: TODO
trusted_values_dict['TestModule__globals'] = {}

#####
'''.format("['a']", str(date.today()), "{'a': mpf('1.0')}")
                   ]
        calc_error_helper(self, message, False)

        ## Broken test: FAILS WITH precision > 30
        # self.calculated_dict = {'a': mpf('1.0')}
        # self.trusted_values_dict_entry = {'a': mpf('1.00000000000000010000000000000')}
        # message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
        #            'DEBUG:root: ...Success: same variables in both dicts.\n',
        #            'DEBUG:root: Comparing all calculated and trusted values...',
        #            'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a' + ': Trusted:    ' + '1.0000000000000001' + '\n',
        #            'DEBUG:root: ...Success: all variables identical.\n'
        #            ]
        # calc_error_helper(self, message, True)

        ## Broken test: FAILS WITH precision > 30
        # self.calculated_dict = {'a': mpf('0.0')}
        # self.trusted_values_dict_entry = {'a': mpf('0.0000000000000001')}
        # message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
        #            'DEBUG:root: ...Success: same variables in both dicts.\n',
        #            'DEBUG:root: Comparing all calculated and trusted values...',
        #            'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' + '\nTestModule: ' + 'a' + ': Trusted:    ' + '1.0e-16' + '\n',
        #            'DEBUG:root: ...Success: all variables identical.\n'
        #            ]
        # calc_error_helper(self, message, True)

        self.calculated_dict = {'b': mpf('0.0')}
        self.trusted_values_dict_entry = {'b': mpf('0.000000000000001')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'DEBUG:root: ...Success: same variables in both dicts.\n',
                   'DEBUG:root: Comparing all calculated and trusted values...',
                   'DEBUG:root:' + '\nTestModule: ' + 'b' + ': Calculated: ' + '0.0' + '\nTestModule: ' + 'b' + ': Trusted:    ' + '1.0e-15' + '\n',
                    'ERROR:root:' +
                   '''
\nVariable(s) {} in module TestModule failed. Please check values.
If you are confident that the newly calculated values are correct, comment out the old trusted values for
TestModule__globals in your trusted_values_dict and copy the following code between the ##### into your trusted_values_dict.
Make sure to fill out the TODO comment describing why the values had to be changed. Then re-run test script.

#####

# Generated on: {}
# Reason for changing values: TODO
trusted_values_dict['TestModule__globals'] = {}

#####
'''.format("['b']", str(date.today()), "{'b': mpf('0.0')}")
                   ]
        calc_error_helper(self, message, False)

        ## Broken test: FAILS WITH precision > 30
        # self.calculated_dict = {'a': mpf('0.0000000000000001')}
        # self.trusted_values_dict_entry = {'a': mpf('0.0')}
        # message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
        #            'DEBUG:root: ...Success: same variables in both dicts.\n',
        #            'DEBUG:root: Comparing all calculated and trusted values...',
        #            'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0e-16' + '\nTestModule: ' + 'a' + ': Trusted:    ' + '0.0' + '\n',
        #            'DEBUG:root: ...Success: all variables identical.\n'
        #            ]
        # calc_error_helper(self, message, True)

        self.calculated_dict = {'alpha': mpf('0.000000000000001')}
        self.trusted_values_dict_entry = {'alpha': mpf('0.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'DEBUG:root: ...Success: same variables in both dicts.\n',
                   'DEBUG:root: Comparing all calculated and trusted values...',
                   'DEBUG:root:' + '\nTestModule: ' + 'alpha' + ': Calculated: ' + '1.0e-15' + '\nTestModule: ' + 'alpha' + ': Trusted:    ' + '0.0' + '\n',
                   'ERROR:root:' +
                   '''
\nVariable(s) {} in module TestModule failed. Please check values.
If you are confident that the newly calculated values are correct, comment out the old trusted values for
TestModule__globals in your trusted_values_dict and copy the following code between the ##### into your trusted_values_dict.
Make sure to fill out the TODO comment describing why the values had to be changed. Then re-run test script.

#####

# Generated on: {}
# Reason for changing values: TODO
trusted_values_dict['TestModule__globals'] = {}

#####
'''.format("['alpha']", str(date.today()), "{'alpha': mpf('1.0e-15')}")
                   ]
        calc_error_helper(self, message, False)

        self.calculated_dict = {'f': mpf('123.012345678901234567890123456')}
        self.trusted_values_dict_entry = {'f': mpf('123.012345678901234567890123456')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'DEBUG:root: ...Success: same variables in both dicts.\n',
                   'DEBUG:root: Comparing all calculated and trusted values...',
                   'DEBUG:root:' + '\nTestModule: ' + 'f' + ': Calculated: ' + '123.012345678901234567890123456' + '\nTestModule: ' + 'f' + ': Trusted:    ' + '123.012345678901234567890123456' + '\n',
                   'DEBUG:root: ...Success: all variables identical.\n'
                   ]
        calc_error_helper(self, message, True)

        self.calculated_dict = {'f': mpf('123.0123456781012345678901')}
        self.trusted_values_dict_entry = {'f': mpf('123.0123456789012345678901')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'DEBUG:root: ...Success: same variables in both dicts.\n',
                   'DEBUG:root: Comparing all calculated and trusted values...',
                   'DEBUG:root:' + '\nTestModule: ' + 'f' + ': Calculated: ' + '123.0123456781012345678901' + '\nTestModule: ' + 'f' + ': Trusted:    ' + '123.0123456789012345678901' + '\n',
                   'ERROR:root:' +
                   '''
\nVariable(s) {} in module TestModule failed. Please check values.
If you are confident that the newly calculated values are correct, comment out the old trusted values for
TestModule__globals in your trusted_values_dict and copy the following code between the ##### into your trusted_values_dict.
Make sure to fill out the TODO comment describing why the values had to be changed. Then re-run test script.

#####

# Generated on: {}
# Reason for changing values: TODO
trusted_values_dict['TestModule__globals'] = {}

#####
'''.format("['f']", str(date.today()), "{'f': mpf('123.0123456781012345678901')}")
                   ]
        calc_error_helper(self, message, False)

        self.calculated_dict = {'a': mpf('0.0'), 'b': mpf('1.0')}
        self.trusted_values_dict_entry = {'a': mpf('0.1'), 'b': mpf('1.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'DEBUG:root: ...Success: same variables in both dicts.\n',
                   'DEBUG:root: Comparing all calculated and trusted values...',
                   'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' + '\nTestModule: ' + 'a' + ': Trusted:    ' + '0.1' + '\n',
                   'DEBUG:root:' + '\nTestModule: ' + 'b' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'b' + ': Trusted:    ' + '1.0' + '\n',
                   'ERROR:root:' +
                   '''
\nVariable(s) {} in module TestModule failed. Please check values.
If you are confident that the newly calculated values are correct, comment out the old trusted values for
TestModule__globals in your trusted_values_dict and copy the following code between the ##### into your trusted_values_dict.
Make sure to fill out the TODO comment describing why the values had to be changed. Then re-run test script.

#####

# Generated on: {}
# Reason for changing values: TODO
trusted_values_dict['TestModule__globals'] = {}

#####
'''.format("['a']", str(date.today()), "{'a': mpf('0.0'), 'b': mpf('1.0')}")
                   ]
        calc_error_helper(self, message, False)

        self.calculated_dict = {'a': mpf('0.0'), 'b': mpf('1.1')}
        self.trusted_values_dict_entry = {'a': mpf('0.1'), 'b': mpf('1.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'DEBUG:root: ...Success: same variables in both dicts.\n',
                   'DEBUG:root: Comparing all calculated and trusted values...',
                   'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' + '\nTestModule: ' + 'a' + ': Trusted:    ' + '0.1' + '\n',
                   'DEBUG:root:' + '\nTestModule: ' + 'b' + ': Calculated: ' + '1.1' + '\nTestModule: ' + 'b' + ': Trusted:    ' + '1.0' + '\n',
                   'ERROR:root:' +
                   '''
\nVariable(s) {} in module TestModule failed. Please check values.
If you are confident that the newly calculated values are correct, comment out the old trusted values for
TestModule__globals in your trusted_values_dict and copy the following code between the ##### into your trusted_values_dict.
Make sure to fill out the TODO comment describing why the values had to be changed. Then re-run test script.

#####

# Generated on: {}
# Reason for changing values: TODO
trusted_values_dict['TestModule__globals'] = {}

#####
'''.format("['a', 'b']", str(date.today()), "{'a': mpf('0.0'), 'b': mpf('1.1')}")
                   ]
        calc_error_helper(self, message, False)

        self.calculated_dict = {'a': mpf('0.0'), 'b': mpf('1.1')}
        self.trusted_values_dict_entry = {'a': mpf('0.0'), 'b': mpf('1.0')}
        message = ['DEBUG:root: Checking that calculated and trusted dicts contain the same variables...',
                   'DEBUG:root: ...Success: same variables in both dicts.\n',
                   'DEBUG:root: Comparing all calculated and trusted values...',
                   'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' + '\nTestModule: ' + 'a' + ': Trusted:    ' + '0.0' + '\n',
                   'DEBUG:root:' + '\nTestModule: ' + 'b' + ': Calculated: ' + '1.1' + '\nTestModule: ' + 'b' + ': Trusted:    ' + '1.0' + '\n',
                   'ERROR:root:' +
                   '''
\nVariable(s) {} in module TestModule failed. Please check values.
If you are confident that the newly calculated values are correct, comment out the old trusted values for
TestModule__globals in your trusted_values_dict and copy the following code between the ##### into your trusted_values_dict.
Make sure to fill out the TODO comment describing why the values had to be changed. Then re-run test script.

#####

# Generated on: {}
# Reason for changing values: TODO
trusted_values_dict['TestModule__globals'] = {}

#####
'''.format("['b']", str(date.today()), "{'a': mpf('0.0'), 'b': mpf('1.1')}")
                   ]
        calc_error_helper(self, message, False)

        logging.info(' All calc_error tests passed.')

    def test_create_dict_string(self):
        from UnitTesting.create_dict_string import create_dict_string

        calculated_dict = {}
        self.assertEqual("{}", create_dict_string(calculated_dict))

        calculated_dict = {'a': 0}
        self.assertEqual("{'a': 0}", create_dict_string(calculated_dict))

        calculated_dict = {'b': mpf('1.0')}
        self.assertEqual("{'b': mpf('1.0')}", create_dict_string(calculated_dict))

        calculated_dict = {'c': mpc(real='1.0', imag='2.0')}
        self.assertEqual("{'c': mpc(real='1.0', imag='2.0')}", create_dict_string(calculated_dict))

        calculated_dict = {'alpha': 4, 'beta': mpf('0.2'), 'gamma': mpc(real='3.14', imag='6.28')}
        self.assertEqual("{'alpha': 4, 'beta': mpf('0.2'), 'gamma': mpc(real='3.14', imag='6.28')}",
                         create_dict_string(calculated_dict))

        calculated_dict = {'beta': mpf('0.2'), 'gamma': mpc(real='3.14', imag='6.28'), 'alpha': 4}
        self.assertEqual("{'alpha': 4, 'beta': mpf('0.2'), 'gamma': mpc(real='3.14', imag='6.28')}",
                         create_dict_string(calculated_dict))

        calculated_dict = {'AZ': mpf('2.4287654'), 'ab': mpc(real='0.0', imag='123.1234123412341234')}
        self.assertEqual("{'ab': mpc(real='0.0', imag='123.1234123412341234'), 'AZ': mpf('2.4287654')}",
                         create_dict_string(calculated_dict))

        calculated_dict = {'betaU[0]': mpf('1.0'), 'betaU[1]': mpf('0.0'), 'betaU[2]': mpf('-0.733235010089164696012176136719')}
        self.assertEqual("{'betaU[0]': mpf('1.0'), 'betaU[1]': mpf('0.0'), 'betaU[2]': mpf('-0.733235010089164696012176136719')}",
                         create_dict_string(calculated_dict))

        logging.info(' All create_dict_string tests passed.')

    def ftest_create_test(self):
        pass

    def test_evaluate_globals(self):
        from UnitTesting.evaluate_globals import evaluate_globals
        import sympy as sp

        logging.getLogger().setLevel('CRITICAL')

        self.module = ''
        self.module_name = 'TestModule'
        self.initialization_string = ''
        self.function = ''
        self.global_list = []
        self.assertRaises(AssertionError, evaluate_globals, self)

        self.module = 'random_string'
        self.module_name = 'TestModule'
        self.initialization_string = ''
        self.function = ''
        self.global_list = []
        self.assertRaises(AssertionError, evaluate_globals, self)

        a, b, c = sp.symbols('a b c')

        self.module = 'module_for_testing'
        self.module_name = 'TestModule'
        self.initialization_string = ''
        self.function = 'function(create_gamma=False)'
        self.global_list = []
        self.assertEqual({}, evaluate_globals(self))

        self.module = 'module_for_testing'
        self.module_name = 'TestModule'
        self.initialization_string = ''
        self.function = 'function(create_gamma=False)'
        self.global_list = ['alpha']
        self.assertEqual({'alpha': a + b + c}, evaluate_globals(self))

        self.module = 'module_for_testing'
        self.module_name = 'TestModule'
        self.initialization_string = ''
        self.function = 'function(create_gamma=False)'
        self.global_list = ['alpha', 'betaU']
        self.assertEqual({'alpha': a + b + c, 'betaU': [0, a**2 + 2*b**2 + c**2, sp.sqrt(a + b)]},
                         evaluate_globals(self))

        self.module = 'module_for_testing'
        self.module_name = 'TestModule'
        self.initialization_string = ''
        self.function = 'function(create_gamma=False)'
        self.global_list = ['gamma']
        self.assertRaises(AttributeError, evaluate_globals, self)

        self.module = 'module_for_testing'
        self.module_name = 'TestModule'
        self.initialization_string = ''
        self.function = 'function(create_gamma=True)'
        self.global_list = ['alpha', 'betaU', 'gamma']
        self.assertEqual({'alpha': a + b + c, 'betaU': [0, a**2 + 2*b**2 + c**2, sp.sqrt(a + b)], 'gamma': sp.atan2(b, a)},
                         evaluate_globals(self))

        self.module = 'module_for_testing'
        self.module_name = 'TestModule'
        self.initialization_string = ''
        self.function = 'function2(create_gamma=False)'
        self.global_list = ['alpha2']
        self.assertRaises(NameError, evaluate_globals, self)

        self.module = 'module_for_testing'
        self.module_name = 'TestModule'
        self.initialization_string = ''
        self.function = 'function2(create_gamma=False)'
        self.global_list = ['betaU2']
        self.assertRaises(NameError, evaluate_globals, self)

        a2, b2, c2 = sp.symbols('a2 b2 c2')

        self.module = 'module_for_testing'
        self.module_name = 'TestModule'
        self.initialization_string = 'import module_for_testing as tm\ntm.init_function2()\n'
        self.function = 'function2(create_gamma=False)'
        self.global_list = ['alpha2', 'betaU2']
        self.assertEqual({'alpha2': a2 + b2 + c2, 'betaU2': [0, a2**2 + 2*b2**2 + c2**2, sp.sqrt(a2 + b2)]},
                         evaluate_globals(self))

        self.module = 'module_for_testing'
        self.module_name = 'TestModule'
        self.initialization_string = ''
        self.function = 'function2(create_gamma=False)'
        self.global_list = ['alpha2', 'betaU2']
        self.assertEqual({'alpha2': a2 + b2 + c2, 'betaU2': [0, a2**2 + 2*b2**2 + c2**2, sp.sqrt(a2 + b2)]},
                         evaluate_globals(self))

        self.module = 'module_for_testing'
        self.module_name = 'TestModule'
        self.initialization_string = ''
        self.function = 'function2(create_gamma=False)'
        self.global_list = ['alpha2', 'betaU2', 'gamma2']
        self.assertRaises(AttributeError, evaluate_globals, self)

        self.module = 'module_for_testing'
        self.module_name = 'TestModule'
        self.initialization_string = ''
        self.function = 'function2(create_gamma=True)'
        self.global_list = ['alpha2', 'betaU2', 'gamma2']
        self.assertEqual({'alpha2': a2 + b2 + c2, 'betaU2': [0, a2 ** 2 + 2 * b2 ** 2 + c2 ** 2, sp.sqrt(a2 + b2)],
                         'gamma2': sp.atan2(b2, a2)}, evaluate_globals(self))

        logging.getLogger().setLevel('INFO')
        logging.info(' All evaluate_globals tests passed.')

    def test_expand_variable_dict(self):
        from UnitTesting.cse_simplify_and_evaluate_sympy_expressions import expand_variable_dict

        variable_dict = dict()
        result_dict = dict()
        self.assertEqual(result_dict, expand_variable_dict(variable_dict))

        variable_dict = {'alpha': 1}
        result_tuple = {'alpha': 1}
        self.assertEqual(result_tuple, expand_variable_dict(variable_dict))

        variable_dict = {'alphaD': [1, 2]}
        result_tuple = {'alphaD[0]': 1, 'alphaD[1]': 2}
        self.assertEqual(result_tuple, expand_variable_dict(variable_dict))

        variable_dict = {'alphaDD': [[1, 2], [4, 3]]}
        result_tuple = {'alphaDD[0][0]': 1, 'alphaDD[0][1]': 2, 'alphaDD[1][0]':4 , 'alphaDD[1][1]':3}
        self.assertEqual(result_tuple, expand_variable_dict(variable_dict))

        variable_dict = {'aDD': [[1, 2, 3, 4, 5], [2, 3, 4, 5, 6], [2, 8, 9, 7, 6], [0, 0, 0, 0, 0], [3, 1, 4, 1, 5]]}
        result_tuple = {'aDD[0][0]': 1, 'aDD[0][1]': 2, 'aDD[0][2]': 3, 'aDD[0][3]': 4, 'aDD[0][4]': 5,
                        'aDD[1][0]': 2, 'aDD[1][1]': 3, 'aDD[1][2]': 4, 'aDD[1][3]': 5, 'aDD[1][4]': 6,
                        'aDD[2][0]': 2, 'aDD[2][1]': 8, 'aDD[2][2]': 9, 'aDD[2][3]': 7, 'aDD[2][4]': 6,
                        'aDD[3][0]': 0, 'aDD[3][1]': 0, 'aDD[3][2]': 0, 'aDD[3][3]': 0, 'aDD[3][4]': 0,
                        'aDD[4][0]': 3, 'aDD[4][1]': 1, 'aDD[4][2]': 4, 'aDD[4][3]': 1, 'aDD[4][4]': 5}
        self.assertEqual(result_tuple, expand_variable_dict(variable_dict))

        variable_dict = {'alpha': 4, 'beta': 5}
        result_tuple = {'alpha': 4, 'beta': 5}
        self.assertEqual(result_tuple, expand_variable_dict(variable_dict))

        variable_dict = {'alphaD': [1, 2], 'beta': 3}
        result_tuple = {'alphaD[0]': 1, 'alphaD[1]': 2, 'beta': 3}
        self.assertEqual(result_tuple, expand_variable_dict(variable_dict))

        logging.info(' All expand_variable_dict tests passed.')

    def test_first_time_print(self):
        from datetime import date

        self.module_name = 'TestModule'
        self.trusted_values_dict_name = 'TestModule__globals'
        self.calculated_dict = {}
        message = '''
Module: TestModule
Please copy the following code between the ##### and paste it into your trusted_values_dict.py file for this module:

#####

# Generated on: {}
trusted_values_dict['TestModule__globals'] = {}

#####
'''.format(date.today(), {})
        first_time_print_helper(self, message=message, write=False)

        self.calculated_dict = {'alpha': 0, 'beta': 1, 'gamma': 3}
        message = '''
Module: TestModule
Please copy the following code between the ##### and paste it into your trusted_values_dict.py file for this module:

#####

# Generated on: {}
trusted_values_dict['TestModule__globals'] = {}

#####
'''.format(date.today(), "{'alpha': 0, 'beta': 1, 'gamma': 3}")
        first_time_print_helper(self, message=message, write=False)

        self.calculated_dict = {'beta': 0, 'gamma': 1, 'alpha': 3}
        message = '''
Module: TestModule
Please copy the following code between the ##### and paste it into your trusted_values_dict.py file for this module:

#####

# Generated on: {}
trusted_values_dict['TestModule__globals'] = {}

#####
'''.format(date.today(), "{'alpha': 3, 'beta': 0, 'gamma': 1}")
        first_time_print_helper(self, message=message, write=False)

        self.calculated_dict = {'x': mpf('0.0'), 'y': mpf('1.23456789012345678912345678912')}
        message = '''
Module: TestModule
Please copy the following code between the ##### and paste it into your trusted_values_dict.py file for this module:

#####

# Generated on: {}
trusted_values_dict['TestModule__globals'] = {}

#####
'''.format(date.today(), "{'x': mpf('0.0'), 'y': mpf('1.23456789012345678912345678912')}")
        first_time_print_helper(self, message=message, write=False)

        self.calculated_dict = {'x': mpf('0.0'), 'y': mpf('1.23456789012345678912345')}
        message = '''
Module: TestModule
Please copy the following code between the ##### and paste it into your trusted_values_dict.py file for this module:

#####

# Generated on: {}
trusted_values_dict['TestModule__globals'] = {}

#####
'''.format(date.today(), "{'x': mpf('0.0'), 'y': mpf('1.23456789012345678912345')}")
        first_time_print_helper(self, message=message, write=False)

        self.calculated_dict = {'AZ': mpf('0.0'), 'ab': mpf('1.0')}
        message = '''
Module: TestModule
Please copy the following code between the ##### and paste it into your trusted_values_dict.py file for this module:

#####

# Generated on: {}
trusted_values_dict['TestModule__globals'] = {}

#####
'''.format(date.today(), "{'ab': mpf('1.0'), 'AZ': mpf('0.0')}")
        first_time_print_helper(self, message=message, write=False)

        self.module_name = 'different_name'
        self.trusted_values_dict_name = 'trusted_values_dict_name'
        self.calculated_dict = {}
        message = '''
Module: different_name
Please copy the following code between the ##### and paste it into your trusted_values_dict.py file for this module:

#####

# Generated on: {}
trusted_values_dict['trusted_values_dict_name'] = {}

#####
'''.format(date.today(), {})
        first_time_print_helper(self, message=message, write=False)

        logging.info(' All first_time_print tests passed.')

    def test_get_variable_dimension(self):
        from UnitTesting.cse_simplify_and_evaluate_sympy_expressions import get_variable_dimension

        rank0 = 4
        rank1 = [rank0, rank0+1, rank0]
        rank2 = [rank1, rank1, rank1]
        rank3 = [rank2]
        self.assertEqual(get_variable_dimension(rank0), (0, 1))
        self.assertEqual(get_variable_dimension(rank1), (1, 3))
        self.assertEqual(get_variable_dimension(rank2), (2, 3))
        self.assertEqual(get_variable_dimension(rank3), (3, 1))

        self.assertEqual(get_variable_dimension(rank3[0]), (2, 3))
        self.assertEqual(get_variable_dimension(rank2[0]), (1, 3))
        self.assertEqual(get_variable_dimension(rank2[1]), (1, 3))
        self.assertEqual(get_variable_dimension([rank2, rank2]), (3, 2))
        self.assertEqual(get_variable_dimension([[[[[rank0]]]]]), (5, 1))

        with self.assertRaises(IndexError):
            get_variable_dimension([])

        logging.info(' All get_variable_dimension tests passed.')

    def ftest_run_test(self):
        import sys
        from UnitTesting.run_test import run_test
        from UnitTesting.setup_trusted_values_dict import setup_trusted_values_dict

        logging.getLogger().setLevel('ERROR')
        setup_trusted_values_dict(sys.path[0])
        logging.getLogger().setLevel('INFO')

        self.module = 'module_for_testing'
        self.module_name = 'test_module'
        self.trusted_values_dict_name = 'test_module__function__create_gamma=False__globals'
        self.function = ''
        self.global_list = []
        self.initialization_string = ''
        self.path = sys.path[0]
        self.assertRaises(AssertionError, run_test, self)

        self.module = 'module_for_testing'
        self.module_name = 'test_module'
        self.trusted_values_dict_name = 'test_module__function__create_gamma=False__globals'
        self.function = 'function(create_gamma=False)'
        self.global_list = ['alpha', 'betaU']
        self.initialization_string = ''

        #self.assertIsNone(run_test(self))

        logging.info(' All run_test tests passed.')

    def ftest_setup_trusted_values_dict(self):

        # Tests the setup_trusted_values_dict and the file-writing portion of first_time_print
        from UnitTesting.setup_trusted_values_dict import setup_trusted_values_dict
        import sys
        import os
        import cmdline_helper as cmd

        logging.getLogger().setLevel('CRITICAL')

        path = sys.path[0]
        full_path = os.path.join(path, 'trusted_values_dict.py')

        cmd.delete_existing_files(full_path)

        self.assertFalse(os.path.exists(full_path))

        setup_trusted_values_dict(path)

        self.assertTrue(os.path.exists(full_path))

        with open(full_path, 'r') as file:
            expected_string = '''from mpmath import mpf, mp, mpc
from UnitTesting.standard_constants import precision

mp.dps = precision
trusted_values_dict = {}
'''
            self.assertEqual(file.read(), expected_string)

        setup_trusted_values_dict(path)

        self.assertTrue(os.path.exists(full_path))

        cmd.delete_existing_files(full_path)

        logging.getLogger().setLevel('INFO')

        logging.info(' All setup_trusted_values_dict tests passed.\n')

    def test_cse_simplify_and_evaluate_sympy_expressions(self):
        from UnitTesting.cse_simplify_and_evaluate_sympy_expressions import cse_simplify_and_evaluate_sympy_expressions
        import random
        from sympy import symbols
        import hashlib

        mp.dps = precision

        self.variable_dict = {}
        self.assertEqual({}, cse_simplify_and_evaluate_sympy_expressions(self))

        M_PI, M_SQRT1_2 = symbols('M_PI M_SQRT1_2')

        self.variable_dict = {'a': 1, 'b': M_PI}
        self.assertEqual({'a': mpf('1.0'), 'b': mpf(pi)}, cse_simplify_and_evaluate_sympy_expressions(self))

        self.variable_dict = {'a': M_PI}
        expected_result = {'a': mpf(pi)}
        actual_result = cse_simplify_and_evaluate_sympy_expressions(self)
        for key, val in expected_result.items():
            self.assertAlmostEqual(val, actual_result[key], 20)

        self.variable_dict = {'b': M_SQRT1_2}
        expected_result = {'b': mpf(1/sqrt(2))}
        actual_result = cse_simplify_and_evaluate_sympy_expressions(self)
        for key, val in expected_result.items():
            self.assertAlmostEqual(val, actual_result[key], 20)

        self.variable_dict = {'alpha': M_PI + M_SQRT1_2}
        expected_result = {'alpha': mpf(pi) + mpf(1/sqrt(2))}
        actual_result = cse_simplify_and_evaluate_sympy_expressions(self)
        for key, val in expected_result.items():
            self.assertAlmostEqual(val, actual_result[key], 20)

        x, y, z = symbols('x y z')
        symbs = {x: 0, y: 0, z: 0}

        for symb in symbs:
            random.seed(int(hashlib.md5(str(symb).encode()).hexdigest(), 16))
            symbs[symb] = mpf(random.random())

        self.variable_dict = {'a': x}
        expected_result = {'a': symbs[x]}
        actual_result = cse_simplify_and_evaluate_sympy_expressions(self)
        for key, val in expected_result.items():
            self.assertAlmostEqual(val, actual_result[key], 20)

        self.variable_dict = {'b': x + y}
        expected_result = {'b': symbs[x] + symbs[y]}
        actual_result = cse_simplify_and_evaluate_sympy_expressions(self)
        for key, val in expected_result.items():
            self.assertAlmostEqual(val, actual_result[key], 20)

        self.variable_dict = {'a': x, 'b': y, 'c': z}
        expected_result = {'a': symbs[x], 'b': symbs[y], 'c': symbs[z]}
        actual_result = cse_simplify_and_evaluate_sympy_expressions(self)
        for key, val in expected_result.items():
            self.assertAlmostEqual(val, actual_result[key], 20)

        self.variable_dict = {'a': x**2, 'b': (x + y)/z}
        expected_result = {'a': symbs[x]**2, 'b': (symbs[x] + symbs[y]) / symbs[z]}
        actual_result = cse_simplify_and_evaluate_sympy_expressions(self)
        for key, val in expected_result.items():
            self.assertAlmostEqual(val, actual_result[key])

        self.variable_dict = {'a': x**2 + 1, 'b': (x**2)**2, 'c': x**2 + y**2, 'd': 1-x**2, 'e': x**2 * z}
        expected_result = {'a': symbs[x]**2 + 1, 'b': symbs[x]**4, 'c': symbs[x]**2 + symbs[y]**2,
                           'd': 1-symbs[x]**2, 'e': symbs[x]**2*symbs[z]}
        actual_result = cse_simplify_and_evaluate_sympy_expressions(self)
        for key, val in expected_result.items():
            self.assertAlmostEqual(val, actual_result[key])

        logging.info(' All cse_simplify_and_evaluate_sympy_expressions tests passed')


def first_time_print_helper(self, message='', write=False):
    from UnitTesting.first_time_print import first_time_print

    if version_info[0] == 2 or version_info[1] < 4:

        from testfixtures import LogCapture
        with LogCapture() as logger:
            first_time_print(self, write)
        logger.check(('root', 'ERROR', message))

    else:

        with self.assertLogs(level='DEBUG') as logger:
            first_time_print(self, write)
        self.assertEqual(logger.output, ['ERROR:root:' + message])


def calc_error_helper(self, message, expected_result):
    from UnitTesting.calc_error import calc_error

    if version_info[0] == 2 or version_info[1] < 4:

        from testfixtures import LogCapture

        tuple_list = []

        for string in message:
            string_tuple = string.split(':', 2)
            string_tuple[0], string_tuple[1] = string_tuple[1], string_tuple[0]
            tuple_list.append(tuple(string_tuple))

        tuple_tuple = tuple(tuple_list)

        with LogCapture() as logger:
            self.assertTrue(expected_result == calc_error(self))
        logger.check(*tuple_tuple)

    else:

        with self.assertLogs(level='DEBUG') as logger:
            self.assertEqual(expected_result, calc_error(self))
        self.assertEqual(logger.output, message)


# Necessary for unittest class to work properly
if __name__ == '__main__':
    unittest.main()
