
import unittest
import logging

# TODO: Change level based on desired amount of output.
# ERROR -> Outputs minimal information -- only when there's an error
# INFO -> Outputs when starting and finishing a module, as well as everything in ERROR
# DEBUG -> Displays all pairs of values being compared, as well as everything in INFO
# NOTSET -> Displays symbolic dictionary for all modules, as well as everything in DEBUG
logging.basicConfig(level=logging.INFO)


class TestFunctions(unittest.TestCase):

    def test_calc_error(self):
        from calc_error import calc_error

        mod = 'TestModule'

        result_dict = dict()
        trusted_dict = dict()

        self.assertEqual(True, calc_error(mod, result_dict, trusted_dict))

        logging.info('All calc_error tests passed.')

    def test_create_trusted_globals_dict(self):
        from create_trusted_globals_dict import create_trusted_globals_dict
        from mpmath import mpf
        from trusted_values_dict import trusted_values_dict

        mod_dict = dict()
        first_times = dict()
        self.assertEqual(dict(), create_trusted_globals_dict(mod_dict, first_times))

        mod_dict = {'mod': []}
        first_times = {'mod': True}
        self.assertEqual({'mod': dict()}, create_trusted_globals_dict(mod_dict, first_times))

        mod_dict = {'mod': ['hello', 'world']}
        first_times = {'mod': True}
        self.assertEqual({'mod': dict()}, create_trusted_globals_dict(mod_dict, first_times))

        with self.assertRaises(KeyError):
            mod_dict = {'mod': ['hello', 'world']}
            first_times = {'mod': False}
            create_trusted_globals_dict(mod_dict, first_times)

        with self.assertRaises(AssertionError):
            mod_dict = {'mod': []}
            first_times = {'mod': True, 'random': False}
            create_trusted_globals_dict(mod_dict, first_times)

        with self.assertRaises(AssertionError):
            mod_dict = {'mod1': 'hello', 'mod2': 'world'}
            first_times = {'mod1': True}
            create_trusted_globals_dict(mod_dict, first_times)

        with self.assertRaises(AssertionError):
            mod_dict = {'a': 1, 'b': 2, 'c': 3}
            first_times = {'d': 4, 'e': 5, 'f': 6}
            create_trusted_globals_dict(mod_dict, first_times)

        with self.assertRaises(AssertionError):
            mod_dict = {'a': 1, 'b': 2, 'c': 3}
            first_times = {'d': 4, 'e': 5, 'f': 6}
            create_trusted_globals_dict(mod_dict, first_times)

        mod_dict = {'BrillLindquist': ['foo', 'bar']}
        first_times = {'BrillLindquist': True}
        self.assertEqual({'BrillLindquist': dict()}, create_trusted_globals_dict(mod_dict, first_times))

        mod_dict = {'BrillLindquist': ['foo', 'bar']}
        first_times = {'BrillLindquist': False}
        self.assertEqual({'BrillLindquist': {'alphaCart': mpf('0.12248333157451517615309'), 'betaCartU[0]': mpf('0.0'), 'betaCartU[1]': mpf('0.0'), 'betaCartU[2]': mpf('0.0'), 'BCartU[0]': mpf('0.0'), 'BCartU[1]': mpf('0.0'), 'BCartU[2]': mpf('0.0'), 'gammaCartDD[0][0]': mpf('66.657039107915231916559'), 'gammaCartDD[0][1]': mpf('0.0'), 'gammaCartDD[0][2]': mpf('0.0'), 'gammaCartDD[1][0]': mpf('0.0'), 'gammaCartDD[1][1]': mpf('66.657039107915231916559'), 'gammaCartDD[1][2]': mpf('0.0'), 'gammaCartDD[2][0]': mpf('0.0'), 'gammaCartDD[2][1]': mpf('0.0'), 'gammaCartDD[2][2]': mpf('66.657039107915231916559'), 'KCartDD[0][0]': mpf('0.0'), 'KCartDD[0][1]': mpf('0.0'), 'KCartDD[0][2]': mpf('0.0'), 'KCartDD[1][0]': mpf('0.0'), 'KCartDD[1][1]': mpf('0.0'), 'KCartDD[1][2]': mpf('0.0'), 'KCartDD[2][0]': mpf('0.0'), 'KCartDD[2][1]': mpf('0.0'), 'KCartDD[2][2]': mpf('0.0')}}
                         , create_trusted_globals_dict(mod_dict, first_times))

        mod_dict = {'BrillLindquist': ['foo', 'bar']}
        first_times = {'BrillLindquist': False}
        self.assertEqual({'BrillLindquist': trusted_values_dict['BrillLindquistGlobals']}
                         , create_trusted_globals_dict(mod_dict, first_times))

        mod_dict = {'BrillLindquist': ['foo', 'bar'], 'mod': ['hello world']}
        first_times = {'BrillLindquist': False, 'mod': True}
        self.assertEqual({'BrillLindquist': trusted_values_dict['BrillLindquistGlobals'], 'mod': dict()}, create_trusted_globals_dict(mod_dict, first_times))

        mod_dict = {'BrillLindquist': 1, 'ShiftedKerrSchild': 2}
        first_times = {'BrillLindquist': False, 'ShiftedKerrSchild': False}
        self.assertEqual({'BrillLindquist': trusted_values_dict['BrillLindquistGlobals']
                         , 'ShiftedKerrSchild': trusted_values_dict['ShiftedKerrSchildGlobals']}
                         , create_trusted_globals_dict(mod_dict, first_times))

        mod_dict = {'mod1': 1, 'mod2': 2, 'mod3': 3}
        first_times = {'mod1': True, 'mod2': True, 'mod3': True}
        self.assertEqual({'mod1': dict(), 'mod2': dict(), 'mod3': dict()},
                         create_trusted_globals_dict(mod_dict, first_times))

        logging.info('\nAll create_trusted_globals_dict tests passed.\n')

    def test_evaluate_globals(self):
        from evaluate_globals import evaluate_globals
        from functions_and_globals import functions_and_globals
        import NRPy_param_funcs as par
        from sympy import sqrt, symbols, sin
        import indexedexp as ixp
        import BSSN.BrillLindquist as BrillLindquist
        import BSSN.StaticTrumpet as StaticTrumpet

        self.assertEqual(dict(), evaluate_globals(dict(), dict()))

        mod_dict = {'BrillLindquist': functions_and_globals(['BrillLindquist(ComputeADMGlobalsOnly = True)'], ['alphaCart', 'betaCartU', 'BCartU', 'gammaCartDD', 'KCartDD'])}

        locs = dict(locals())

        thismodule = "Brill-Lindquist"
        BH1_posn_x, BH1_posn_y, BH1_posn_z = par.Cparameters("REAL", thismodule,
                                                             ["BH1_posn_x", "BH1_posn_y", "BH1_posn_z"])
        BH1_mass = par.Cparameters("REAL", thismodule, ["BH1_mass"])
        BH2_posn_x, BH2_posn_y, BH2_posn_z = par.Cparameters("REAL", thismodule,
                                                             ["BH2_posn_x", "BH2_posn_y", "BH2_posn_z"])
        BH2_mass = par.Cparameters("REAL", thismodule, ["BH2_mass"])
        Cartxyz = ixp.declarerank1("Cartxyz")
        Cartxyz0, Cartxyz1, Cartxyz2 = Cartxyz

        result_dict = {'BrillLindquist': {'alphaCart': (BH1_mass/(2*sqrt((-BH1_posn_x + Cartxyz0)**2 + (-BH1_posn_y + Cartxyz1)**2 + (-BH1_posn_z + Cartxyz2)**2)) + BH2_mass/(2*sqrt((-BH2_posn_x + Cartxyz0)**2 + (-BH2_posn_y + Cartxyz1)**2 + (-BH2_posn_z + Cartxyz2)**2)) + 1)**(-2), 'betaCartU': [0, 0, 0], 'BCartU': [0, 0, 0], 'gammaCartDD': [[(BH1_mass/(2*sqrt((-BH1_posn_x + Cartxyz0)**2 + (-BH1_posn_y + Cartxyz1)**2 + (-BH1_posn_z + Cartxyz2)**2)) + BH2_mass/(2*sqrt((-BH2_posn_x + Cartxyz0)**2 + (-BH2_posn_y + Cartxyz1)**2 + (-BH2_posn_z + Cartxyz2)**2)) + 1)**4, 0, 0], [0, (BH1_mass/(2*sqrt((-BH1_posn_x + Cartxyz0)**2 + (-BH1_posn_y + Cartxyz1)**2 + (-BH1_posn_z + Cartxyz2)**2)) + BH2_mass/(2*sqrt((-BH2_posn_x + Cartxyz0)**2 + (-BH2_posn_y + Cartxyz1)**2 + (-BH2_posn_z + Cartxyz2)**2)) + 1)**4, 0], [0, 0, (BH1_mass/(2*sqrt((-BH1_posn_x + Cartxyz0)**2 + (-BH1_posn_y + Cartxyz1)**2 + (-BH1_posn_z + Cartxyz2)**2)) + BH2_mass/(2*sqrt((-BH2_posn_x + Cartxyz0)**2 + (-BH2_posn_y + Cartxyz1)**2 + (-BH2_posn_z + Cartxyz2)**2)) + 1)**4]], 'KCartDD': [[0, 0, 0], [0, 0, 0], [0, 0, 0]]}}

        self.assertEqual(result_dict, evaluate_globals(mod_dict, locs))

        mod_dict = {'BrillLindquist': functions_and_globals(['BrillLindquist(ComputeADMGlobalsOnly = True)'], ['alphaCart', 'betaCartU', 'BCartU', 'gammaCartDD', 'KCartDD']),
                    'StaticTrumpet': functions_and_globals(['StaticTrumpet(ComputeADMGlobalsOnly = True)'], ['alphaSph', 'betaSphU', 'BSphU', 'gammaSphDD', 'KSphDD'])}

        locs = locals()

        r, th, ph = symbols('r th ph', real=True)
        M = par.Cparameters("REAL", thismodule, ["M"])

        result_dict = {'BrillLindquist': {'alphaCart': (BH1_mass/(2*sqrt((-BH1_posn_x + Cartxyz0)**2 + (-BH1_posn_y + Cartxyz1)**2 + (-BH1_posn_z + Cartxyz2)**2)) + BH2_mass/(2*sqrt((-BH2_posn_x + Cartxyz0)**2 + (-BH2_posn_y + Cartxyz1)**2 + (-BH2_posn_z + Cartxyz2)**2)) + 1)**(-2), 'betaCartU': [0, 0, 0], 'BCartU': [0, 0, 0], 'gammaCartDD': [[(BH1_mass/(2*sqrt((-BH1_posn_x + Cartxyz0)**2 + (-BH1_posn_y + Cartxyz1)**2 + (-BH1_posn_z + Cartxyz2)**2)) + BH2_mass/(2*sqrt((-BH2_posn_x + Cartxyz0)**2 + (-BH2_posn_y + Cartxyz1)**2 + (-BH2_posn_z + Cartxyz2)**2)) + 1)**4, 0, 0], [0, (BH1_mass/(2*sqrt((-BH1_posn_x + Cartxyz0)**2 + (-BH1_posn_y + Cartxyz1)**2 + (-BH1_posn_z + Cartxyz2)**2)) + BH2_mass/(2*sqrt((-BH2_posn_x + Cartxyz0)**2 + (-BH2_posn_y + Cartxyz1)**2 + (-BH2_posn_z + Cartxyz2)**2)) + 1)**4, 0], [0, 0, (BH1_mass/(2*sqrt((-BH1_posn_x + Cartxyz0)**2 + (-BH1_posn_y + Cartxyz1)**2 + (-BH1_posn_z + Cartxyz2)**2)) + BH2_mass/(2*sqrt((-BH2_posn_x + Cartxyz0)**2 + (-BH2_posn_y + Cartxyz1)**2 + (-BH2_posn_z + Cartxyz2)**2)) + 1)**4]], 'KCartDD': [[0, 0, 0], [0, 0, 0], [0, 0, 0]]},
                       'StaticTrumpet': {'alphaSph': r/(M + r), 'betaSphU': [M*r/(M + r)**2, 0, 0], 'BSphU': [0, 0, 0], 'gammaSphDD': [[(M/r + 1)**2, 0, 0], [0, r**2*(M/r + 1)**2, 0], [0, 0, r**2*(M/r + 1)**2*sin(th)**2]], 'KSphDD': [[-M/r**2, 0, 0], [0, M, 0], [0, 0, M*sin(th)**2]]}}

        self.assertEqual(result_dict, evaluate_globals(mod_dict, locs))

        logging.info('\nAll evaluate_globals tests passed.\n')

    def test_expand_variable_dict(self):
        from expand_variable_dict import expand_variable_dict

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

        logging.info('\nAll expand_variable_dict tests passed.\n')

    def test_functions_and_globals(self):
        from functions_and_globals import functions_and_globals

        basic_function_list = ['func1(), func2()']
        basic_global_list = ['x', 'y', 'z']

        self.assertEqual(functions_and_globals([], []), {'function_list': [], 'global_list': []})

        self.assertEqual(functions_and_globals([], basic_global_list),
                         {'function_list': [], 'global_list': basic_global_list})
        self.assertEqual(functions_and_globals(basic_function_list, []),
                         {'function_list': basic_function_list, 'global_list': []})
        self.assertEqual(functions_and_globals(basic_function_list, basic_global_list),
                         {'function_list': basic_function_list, 'global_list': basic_global_list})

        with self.assertRaises(AssertionError):
            functions_and_globals([1, 'hello', 'world'], [])

        with self.assertRaises(AssertionError):
            functions_and_globals(['hello', 'world'], [2])

        with self.assertRaises(AssertionError):
            functions_and_globals(['hello', 'world', 42], basic_global_list)

        with self.assertRaises(AssertionError):
            functions_and_globals('function()', [])

        with self.assertRaises(AssertionError):
            functions_and_globals([], 'glob')

        logging.info('\nAll functions_and_globals tests passed.\n')

    def test_get_variable_dimension(self):
        from get_variable_dimension import get_variable_dimension

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

        logging.info('\nAll get_variable_dimension tests passed.\n')

    def test_is_first_time(self):
        from is_first_time import is_first_time

        mod_dict = {'BrillLindquist': 'Hello World'}
        fake_mod_dict = {'fake_module': 'Goodbye World'}

        self.assertEqual(is_first_time(dict()), dict())

        self.assertEqual(is_first_time(mod_dict), {'BrillLindquist': False})
        self.assertEqual(is_first_time(fake_mod_dict), {'fake_module': True})

        large_mod_dict = {'BrillLindquist': 'Hello World', 'fake_module': 'Goodbye World'}

        self.assertEqual(is_first_time(large_mod_dict), {'BrillLindquist': False, 'fake_module': True})

        mod_dict_wrong_capitalization = {'brillLindquist': 2}

        self.assertEqual(is_first_time(mod_dict_wrong_capitalization), {'brillLindquist': True})

        logging.info('\nAll is_first_time tests passed.\n')

    def test_run_test(self):
        from run_test import run_test

        mod_dict = {}
        with self.assertRaises(AssertionError):
            run_test(self, mod_dict, locals())

    def test_var_dict_to_value_dict(self):
        from var_dict_to_value_dict import var_dict_to_value_dict
        from mpmath import mpf, sqrt, mp
        from random import random, seed
        from trusted_values_dict import trusted_values_dict
        from sympy.abc import x, y, z

        seed(trusted_values_dict['seed'])
        mp.dps = trusted_values_dict['precision']

        var_dict = dict()
        self.assertEqual(dict(), var_dict_to_value_dict(var_dict))

        first_val = mpf(sqrt(random()))
        second_val = mpf(sqrt(random()))
        third_val = mpf(sqrt(random()))

        var_dict = {'x': x}
        result_dict = {'x': first_val}
        self.assertEqual(result_dict, var_dict_to_value_dict(var_dict))

        var_dict = {'y': y}
        result_dict = {'y': first_val}
        self.assertEqual(result_dict, var_dict_to_value_dict(var_dict))

        var_dict = {'x': x, 'y': y}
        result_dict = {'x': first_val, 'y': second_val}
        self.assertEqual(result_dict, var_dict_to_value_dict(var_dict))

        var_dict = {'y': y, 'x': x}
        result_dict = {'y': second_val, 'x': first_val}
        self.assertEqual(result_dict, var_dict_to_value_dict(var_dict))

        var_dict = {'x': x, 'y': y, 'z': z}
        result_dict = {'x': first_val, 'y': second_val, 'z': third_val}
        self.assertEqual(result_dict, var_dict_to_value_dict(var_dict))

        var_dict = {'x+y': x+y}
        result_dict = {'x+y': first_val+second_val}
        self.assertEqual(result_dict, var_dict_to_value_dict(var_dict))

        var_dict = {'x/y': x/y, 'long_expression': (x+y)**2/(x*y)}
        result_dict = {'x/y': first_val/second_val, 'long_expression': (first_val+second_val)**2/(first_val*second_val)}
        self.assertEqual(result_dict, var_dict_to_value_dict(var_dict))

        var_dict = {'x-x': x-x}
        result_dict = {'x-x': mpf(0.0)}
        self.assertEqual(result_dict, var_dict_to_value_dict(var_dict))

        with self.assertRaises(AttributeError):
            var_dict = {'x': 0}
            var_dict_to_value_dict(var_dict)

        logging.info('\nAll var_dict_to_value_dict tests passed\n')


# Necessary for unittest class to work properly
if __name__ == '__main__':
    unittest.main()
