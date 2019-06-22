
import unittest
import logging
from sys import version_info

# TODO: Change level based on desired amount of output.
# ERROR -> Outputs minimal information -- only when there's an error
# INFO -> Outputs when starting and finishing a module, as well as everything in ERROR
# DEBUG -> Displays all pairs of values being compared, as well as everything in INFO
# NOTSET -> Displays symbolic dictionary for all modules, as well as everything in DEBUG
logging.basicConfig(level=logging.INFO)


class TestFunctions(unittest.TestCase):

    # If Python version 2.x or 3.y, y < 4
    if version_info[0] == 2 or version_info[1] < 4:
        def test_calc_error(self):
            from UnitTesting.calc_error import calc_error
            from testfixtures import LogCapture
            from mpmath import mpf, mp
            from datetime import date
            from UnitTesting.trusted_values_dict import trusted_values_dict

            mp.dps = trusted_values_dict['precision']

            mod = 'TestModule'

            calculated_dict = {}
            trusted_dict = {}

            self.assertEqual(True, calc_error(mod, calculated_dict, trusted_dict))

            calculated_dict = {'a': mpf('1.0')}
            trusted_dict =    {}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'ERROR', '\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.'),
                ('root', 'ERROR', '\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['a']")
            )

            calculated_dict = {'a': mpf('1.0'), 'b': mpf('2.0')}
            trusted_dict =    {}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'ERROR', '\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.'),
                ('root', 'ERROR', '\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['a', 'b']")
            )

            calculated_dict = {'a': mpf('1.0')}
            trusted_dict =    {'b': mpf('2.0')}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'ERROR', '\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.'),
                ('root', 'ERROR', '\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['a']"),
                ('root', 'ERROR', '\n\tTrusted Dictionary variables not in Calculated Dictionary: \n\t' + "['b']")
            )

            calculated_dict = {'a': mpf('2.0'), 'b': mpf('3.0')}
            trusted_dict =    {'c': mpf('1.0')}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'ERROR', '\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.'),
                ('root', 'ERROR', '\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['a', 'b']"),
                ('root', 'ERROR', '\n\tTrusted Dictionary variables not in Calculated Dictionary: \n\t' + "['c']")
            )

            calculated_dict = {'a': mpf('2.0'), 'b': mpf('3.0')}
            trusted_dict =    {'a': mpf('1.0'), 'c': mpf('4.0')}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'ERROR', '\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.'),
                ('root', 'ERROR', '\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['b']"),
                ('root', 'ERROR', '\n\tTrusted Dictionary variables not in Calculated Dictionary: \n\t' + "['c']")
            )

            calculated_dict = {'a': mpf('1.0'), 'b': mpf('2.0')}
            trusted_dict =    {'c': mpf('3.0'), 'd': mpf('4.0')}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'ERROR', '\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.'),
                ('root', 'ERROR', '\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['a', 'b']"),
                ('root', 'ERROR', '\n\tTrusted Dictionary variables not in Calculated Dictionary: \n\t' + "['c', 'd']")
            )

            calculated_dict = {'b': mpf('1.0'), 'a': mpf('2.0')}
            trusted_dict =    {'c': mpf('3.0')}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'ERROR', '\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.'),
                ('root', 'ERROR', '\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['a', 'b']"),
                ('root', 'ERROR', '\n\tTrusted Dictionary variables not in Calculated Dictionary: \n\t' + "['c']")
            )

            calculated_dict = {'a': mpf('1.0')}
            trusted_dict =    {'a': mpf('1.0')}

            with LogCapture() as log:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a' +
                 ': Trusted:    ' + '1.0' + '\n')
            )

            calculated_dict = {'a': mpf('1.0'), 'b': mpf('2.0')}
            trusted_dict = {'a': mpf('1.0'), 'b': mpf('2.0')}

            with LogCapture() as log:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a' +
                 ': Trusted:    ' + '1.0' + '\n'),
                ('root', 'DEBUG', '\nTestModule: ' + 'b' + ': Calculated: ' + '2.0' + '\nTestModule: ' + 'b' +
                 ': Trusted:    ' + '2.0' + '\n')
            )

            calculated_dict = {'a': mpf('1.0')}
            trusted_dict =    {'a': mpf('2.0')}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a' +
                 ': Trusted:    ' + '2.0' + '\n'),
                ('root', 'INFO', '\n\nVariable ' + "'a'" + ' in module TestModule failed. Please check values.\n\nIf you '
                 + 'are confident that the newly calculated values are correct, comment out the old trusted values for ' +
                 "'TestModuleGlobals' in trusted_values_dict and copy the following code between the ##### into " +
                 'trusted_values_dict. Make sure to fill out the TODO comment describing why the values had to be changed.'
                 + ' Then re-run test script.\n#####\n\n# Generated on: ' + str(date.today()) + '\n# Reason for changing' +
                 " values: TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'a'" + ": mpf('" + '1.0' + "')}\n\n" +
                 '#####')
            )

            calculated_dict = {'a': mpf('0.0')}
            trusted_dict =    {'a': mpf('0.0')}

            with LogCapture() as log:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' + '\nTestModule: ' + 'a'
                 + ': Trusted:    ' + '0.0' + '\n')
            )

            calculated_dict = {'a': mpf('1.0')}
            trusted_dict =    {'a': mpf('1.00000000000000010000000000000')}

            with LogCapture() as log:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a'
                 + ': Trusted:    ' + '1.0000000000000001' + '\n')
            )

            calculated_dict = {'a': mpf('0.0')}
            trusted_dict =    {'a': mpf('0.0000000000000001')}

            with LogCapture() as log:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' + '\nTestModule: ' + 'a'
                 + ': Trusted:    ' + '1.0e-16' + '\n')
            )

            calculated_dict = {'b': mpf('0.0')}
            trusted_dict =    {'b': mpf('0.000000000000001')}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'b' + ': Calculated: ' + '0.0' + '\nTestModule: ' + 'b'
                 + ': Trusted:    ' + '1.0e-15' + '\n'),
                ('root', 'INFO', '\n\nVariable ' + "'b'" + ' in module TestModule failed. Please check values.\n\nIf you '
                 + 'are confident that the newly calculated values are correct, comment out the old trusted values for ' +
                 "'TestModuleGlobals' in trusted_values_dict and copy the following code between the ##### into " +
                 'trusted_values_dict. Make sure to fill out the TODO comment describing why the values had to be changed.'
                 + ' Then re-run test script.\n#####\n\n# Generated on: ' + str(date.today()) + '\n# Reason for changing' +
                 " values: TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'b'" + ": mpf('" + '0.0' + "')}\n\n" +
                 '#####')
            )

            calculated_dict = {'a': mpf('0.0000000000000001')}
            trusted_dict =    {'a': mpf('0.0')}

            with LogCapture() as log:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0e-16' + '\nTestModule: ' + 'a'
                 + ': Trusted:    ' + '0.0' + '\n')
            )

            calculated_dict = {'alpha': mpf('0.000000000000001')}
            trusted_dict =    {'alpha': mpf('0.0')}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'alpha' + ': Calculated: ' + '1.0e-15' + '\nTestModule: ' + 'alpha'
                 + ': Trusted:    ' + '0.0' + '\n'),
                ('root',
                 'INFO', '\n\nVariable ' + "'alpha'" + ' in module TestModule failed. Please check values.\n\nIf you '
                 + 'are confident that the newly calculated values are correct, comment out the old trusted values for ' +
                 "'TestModuleGlobals' in trusted_values_dict and copy the following code between the ##### into " +
                 'trusted_values_dict. Make sure to fill out the TODO comment describing why the values had to be changed.'
                 + ' Then re-run test script.\n#####\n\n# Generated on: ' + str(date.today()) + '\n# Reason for changing' +
                 " values: TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'alpha'" + ": mpf('" + '1.0e-15' +
                 "')}\n\n#####")
            )

            calculated_dict = {'f': mpf('123.012345678901234567890123456789')}
            trusted_dict =    {'f': mpf('123.012345678901234567890123456789')}

            with LogCapture() as log:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'f' + ': Calculated: ' + '123.012345678901234567890123457' +
                 '\nTestModule: ' + 'f' + ': Trusted:    ' + '123.012345678901234567890123457' + '\n')
            )

            calculated_dict = {'f': mpf('123.012345678101234567890123457')}
            trusted_dict =    {'f': mpf('123.012345678901234567890123457')}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'f' + ': Calculated: ' + '123.012345678101234567890123457' +
                 '\nTestModule: ' + 'f' + ': Trusted:    ' + '123.012345678901234567890123457' + '\n'),
                ('root',
                 'INFO', '\n\nVariable ' + "'f'" + ' in module TestModule failed. Please check values.\n\nIf you '
                 + 'are confident that the newly calculated values are correct, comment out the old trusted values for ' +
                 "'TestModuleGlobals' in trusted_values_dict and copy the following code between the ##### into " +
                 'trusted_values_dict. Make sure to fill out the TODO comment describing why the values had to be changed.'
                 + ' Then re-run test script.\n#####\n\n# Generated on: ' + str(date.today()) + '\n# Reason for changing' +
                 " values: TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'f'" + ": mpf('" +
                 '123.012345678101234567890123457' + "')}\n\n#####")
            )

            calculated_dict = {'f': mpf('123.012345678901234567890123456')}
            trusted_dict =    {'f': mpf('123.012345678901234567890123457')}

            with LogCapture() as log:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'f' + ': Calculated: ' + '123.012345678901234567890123456' +
                 '\nTestModule: ' + 'f' + ': Trusted:    ' + '123.012345678901234567890123457' + '\n'),
            )

            calculated_dict = {'a': mpf('0.0'), 'b': mpf('1.0')}
            trusted_dict =    {'a': mpf('0.0'), 'b': mpf('1.0')}

            with LogCapture() as log:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' +
                 '\nTestModule: ' + 'a' + ': Trusted:    ' + '0.0' + '\n'),
                ('root', 'DEBUG', '\nTestModule: ' + 'b' + ': Calculated: ' + '1.0' +
                 '\nTestModule: ' + 'b' + ': Trusted:    ' + '1.0' + '\n'),
            )

            calculated_dict = {'a': mpf('0.0'), 'b': mpf('1.0')}
            trusted_dict =    {'a': mpf('0.1'), 'b': mpf('1.0')}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' +
                 '\nTestModule: ' + 'a' + ': Trusted:    ' + '0.1' + '\n'),
                ('root',
                 'INFO', '\n\nVariable ' + "'a'" + ' in module TestModule failed. Please check values.\n\nIf you '
                 + 'are confident that the newly calculated values are correct, comment out the old trusted values for ' +
                 "'TestModuleGlobals' in trusted_values_dict and copy the following code between the ##### into " +
                 'trusted_values_dict. Make sure to fill out the TODO comment describing why the values had to be changed.'
                 + ' Then re-run test script.\n#####\n\n# Generated on: ' + str(date.today()) + '\n# Reason for changing' +
                 " values: TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'a'" + ": mpf('" + '0.0' + "')" + ', ' +
                 "'b'" + ": mpf('1.0')" + "}\n\n#####")
            )

            calculated_dict = {'a': mpf('0.0'), 'b': mpf('1.1')}
            trusted_dict =    {'a': mpf('0.1'), 'b': mpf('1.0')}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' +
                 '\nTestModule: ' + 'a' + ': Trusted:    ' + '0.1' + '\n'),
                ('root',
                 'INFO', '\n\nVariable ' + "'a'" + ' in module TestModule failed. Please check values.\n\nIf you '
                 + 'are confident that the newly calculated values are correct, comment out the old trusted values for ' +
                 "'TestModuleGlobals' in trusted_values_dict and copy the following code between the ##### into " +
                 'trusted_values_dict. Make sure to fill out the TODO comment describing why the values had to be changed.'
                 + ' Then re-run test script.\n#####\n\n# Generated on: ' + str(date.today()) + '\n# Reason for changing' +
                 " values: TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'a'" + ": mpf('" + '0.0' + "')" + ', ' +
                 "'b'" + ": mpf('1.1')" + "}\n\n#####")
            )

            calculated_dict = {'a': mpf('0.0'), 'b': mpf('1.1')}
            trusted_dict =    {'a': mpf('0.0'), 'b': mpf('1.0')}

            with LogCapture() as log:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            log.check(
                ('root', 'DEBUG', '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' +
                 '\nTestModule: ' + 'a' + ': Trusted:    ' + '0.0' + '\n'),
                ('root', 'DEBUG', '\nTestModule: ' + 'b' + ': Calculated: ' + '1.1' +
                 '\nTestModule: ' + 'b' + ': Trusted:    ' + '1.0' + '\n'),
                ('root',
                 'INFO', '\n\nVariable ' + "'b'" + ' in module TestModule failed. Please check values.\n\nIf you '
                 + 'are confident that the newly calculated values are correct, comment out the old trusted values for ' +
                 "'TestModuleGlobals' in trusted_values_dict and copy the following code between the ##### into " +
                 'trusted_values_dict. Make sure to fill out the TODO comment describing why the values had to be changed.'
                 + ' Then re-run test script.\n#####\n\n# Generated on: ' + str(date.today()) + '\n# Reason for changing' +
                 " values: TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'a'" + ": mpf('" + '0.0' + "')" + ', ' +
                 "'b'" + ": mpf('1.1')" + "}\n\n#####")
            )

            logging.info('\nAll calc_error tests passed.\n')
    else:
        def test_calc_error(self):
            from UnitTesting.calc_error import calc_error
            from mpmath import mpf, mp
            from datetime import date
            from UnitTesting.trusted_values_dict import trusted_values_dict

            mp.dps = trusted_values_dict['precision']

            mod = 'TestModule'

            calculated_dict = {}
            trusted_dict = {}

            self.assertEqual(True, calc_error(mod, calculated_dict, trusted_dict))

            calculated_dict = {'a': mpf('1.0')}
            trusted_dict = {}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                ['ERROR:root:\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.',
                 'ERROR:root:\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['a']"])

            calculated_dict = {'a': mpf('1.0'), 'b': mpf('2.0')}
            trusted_dict = {}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                 'ERROR:root:\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.',
                 'ERROR:root:\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['a', 'b']"
                ])

            calculated_dict = {'a': mpf('1.0')}
            trusted_dict = {'b': mpf('2.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'ERROR:root:\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.',
                    'ERROR:root:\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['a']",
                    'ERROR:root:\n\tTrusted Dictionary variables not in Calculated Dictionary: \n\t' + "['b']"
                ]
            )

            calculated_dict = {'a': mpf('2.0'), 'b': mpf('3.0')}
            trusted_dict = {'c': mpf('1.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'ERROR:root:\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.',
                    'ERROR:root:\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['a', 'b']",
                    'ERROR:root:\n\tTrusted Dictionary variables not in Calculated Dictionary: \n\t' + "['c']"
                ]
            )

            calculated_dict = {'a': mpf('2.0'), 'b': mpf('3.0')}
            trusted_dict = {'a': mpf('1.0'), 'c': mpf('4.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'ERROR:root:\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.',
                    'ERROR:root:\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['b']",
                    'ERROR:root:\n\tTrusted Dictionary variables not in Calculated Dictionary: \n\t' + "['c']"
                ]
            )

            calculated_dict = {'a': mpf('1.0'), 'b': mpf('2.0')}
            trusted_dict = {'c': mpf('3.0'), 'd': mpf('4.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'ERROR:root:\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.',
                    'ERROR:root:\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['a', 'b']",
                    'ERROR:root:\n\tTrusted Dictionary variables not in Calculated Dictionary: \n\t' + "['c', 'd']"
                ]
            )

            calculated_dict = {'b': mpf('1.0'), 'a': mpf('2.0')}
            trusted_dict = {'c': mpf('3.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'ERROR:root:\n\tTestModule: Calculated dictionary and trusted dictionary have different variables.',
                    'ERROR:root:\n\tCalculated Dictionary variables not in Trusted Dictionary: \n\t' + "['a', 'b']",
                    'ERROR:root:\n\tTrusted Dictionary variables not in Calculated Dictionary: \n\t' + "['c']"
                ]
            )

            calculated_dict = {'a': mpf('1.0')}
            trusted_dict = {'a': mpf('1.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                 'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a' +
                 ': Trusted:    ' + '1.0' + '\n'
                ]
            )

            calculated_dict = {'a': mpf('1.0'), 'b': mpf('2.0')}
            trusted_dict = {'a': mpf('1.0'), 'b': mpf('2.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                 'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a' +
                 ': Trusted:    ' + '1.0' + '\n',
                 'DEBUG:root:' + '\nTestModule: ' + 'b' + ': Calculated: ' + '2.0' + '\nTestModule: ' + 'b' +
                 ': Trusted:    ' + '2.0' + '\n'
                ]
            )

            calculated_dict = {'a': mpf('1.0')}
            trusted_dict = {'a': mpf('2.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a' +
                    ': Trusted:    ' + '2.0' + '\n',
                    'INFO:root:' +
                    '\n\nVariable ' + "'a'" + ' in module TestModule failed. Please check values.\n\nIf you are ' +
                    'confident that the newly calculated values are correct, comment out the old trusted values for ' +
                    "'TestModuleGlobals' in trusted_values_dict and copy the following code between the ##### into " +
                    'trusted_values_dict. Make sure to fill out the TODO comment describing why the values had to be ' +
                    'changed. Then re-run test script.\n#####\n\n# Generated on: ' + str(date.today()) +
                    '\n# Reason for changing' + " values: TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'a'" +
                    ": mpf('" + '1.0' + "')}\n\n#####"
                ]
            )

            calculated_dict = {'a': mpf('0.0')}
            trusted_dict = {'a': mpf('0.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' + '\nTestModule: ' + 'a' +
                    ': Trusted:    ' + '0.0' + '\n'
                ]
            )

            calculated_dict = {'a': mpf('1.0')}
            trusted_dict = {'a': mpf('1.00000000000000010000000000000')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0' + '\nTestModule: ' + 'a' +
                    ': Trusted:    ' + '1.0000000000000001' + '\n'
                ]
            )

            calculated_dict = {'a': mpf('0.0')}
            trusted_dict = {'a': mpf('0.0000000000000001')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' + '\nTestModule: ' + 'a' +
                    ': Trusted:    ' + '1.0e-16' + '\n'
                ]
            )

            calculated_dict = {'b': mpf('0.0')}
            trusted_dict = {'b': mpf('0.000000000000001')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'b' + ': Calculated: ' + '0.0' + '\nTestModule: ' + 'b'
                    + ': Trusted:    ' + '1.0e-15' + '\n',
                    'INFO:root:' + '\n\nVariable ' + "'b'" + ' in module TestModule failed. Please check values.\n\nIf'
                    + ' you are confident that the newly calculated values are correct, comment out the old trusted ' +
                    'values for ' + "'TestModuleGlobals' in trusted_values_dict and copy the following code between " +
                    "the ##### into " + 'trusted_values_dict. Make sure to fill out the TODO comment describing why ' +
                    'the values had to be changed. Then re-run test script.\n#####\n\n# Generated on: ' +
                    str(date.today()) + '\n# Reason for changing' +
                    " values: TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'b'" + ": mpf('" + '0.0' +
                    "')}\n\n#####"
                ]
            )

            calculated_dict = {'a': mpf('0.0000000000000001')}
            trusted_dict = {'a': mpf('0.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '1.0e-16' + '\nTestModule: ' + 'a'
                    + ': Trusted:    ' + '0.0' + '\n'
                ]
            )

            calculated_dict = {'alpha': mpf('0.000000000000001')}
            trusted_dict = {'alpha': mpf('0.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'alpha' + ': Calculated: ' + '1.0e-15' + '\nTestModule: ' +
                    'alpha' + ': Trusted:    ' + '0.0' + '\n',
                    'INFO:root:' + '\n\nVariable ' + "'alpha'" + ' in module TestModule failed. Please check values.' +
                    '\n\nIf you are confident that the newly calculated values are correct, comment out the old ' +
                    'trusted values for ' + "'TestModuleGlobals' in trusted_values_dict and copy the following code" +
                    " between the ##### into " + 'trusted_values_dict. Make sure to fill out the TODO comment' +
                    ' describing why the values had to be changed. Then re-run test script.\n#####\n\n# Generated on: '
                    + str(date.today()) + '\n# Reason for changing values:' +
                    " TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'alpha'" + ": mpf('" + '1.0e-15' +
                    "')}\n\n#####"
                ]
            )

            calculated_dict = {'f': mpf('123.012345678901234567890123456789')}
            trusted_dict = {'f': mpf('123.012345678901234567890123456789')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'f' + ': Calculated: ' + '123.012345678901234567890123457' +
                    '\nTestModule: ' + 'f' + ': Trusted:    ' + '123.012345678901234567890123457' + '\n'
                ]
            )

            calculated_dict = {'f': mpf('123.012345678101234567890123457')}
            trusted_dict = {'f': mpf('123.012345678901234567890123457')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'f' + ': Calculated: ' + '123.012345678101234567890123457' +
                    '\nTestModule: ' + 'f' + ': Trusted:    ' + '123.012345678901234567890123457' + '\n',
                    'INFO:root:' + '\n\nVariable ' + "'f'" + ' in module TestModule failed. Please check values.\n\n' +
                    'If you are confident that the newly calculated values are correct, comment out the old trusted ' +
                    'values for ' + "'TestModuleGlobals' in trusted_values_dict and copy the following code between" +
                    ' the ##### into trusted_values_dict. Make sure to fill out the TODO comment describing why the ' +
                    'values had to be changed. Then re-run test script.\n#####\n\n# Generated on: '
                    + str(date.today()) + '\n# Reason for changing values: ' +
                    "TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'f'" + ": mpf('" +
                    '123.012345678101234567890123457' + "')}\n\n#####"
                ]
            )

            calculated_dict = {'f': mpf('123.012345678901234567890123456')}
            trusted_dict = {'f': mpf('123.012345678901234567890123457')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'f' + ': Calculated: ' + '123.012345678901234567890123456' +
                    '\nTestModule: ' + 'f' + ': Trusted:    ' + '123.012345678901234567890123457' + '\n'
                ]
            )

            calculated_dict = {'a': mpf('0.0'), 'b': mpf('1.0')}
            trusted_dict = {'a': mpf('0.0'), 'b': mpf('1.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertTrue(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' +
                    '\nTestModule: ' + 'a' + ': Trusted:    ' + '0.0' + '\n',
                    'DEBUG:root:' + '\nTestModule: ' + 'b' + ': Calculated: ' + '1.0' +
                    '\nTestModule: ' + 'b' + ': Trusted:    ' + '1.0' + '\n'
                ]
            )

            calculated_dict = {'a': mpf('0.0'), 'b': mpf('1.0')}
            trusted_dict = {'a': mpf('0.1'), 'b': mpf('1.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' +
                    '\nTestModule: ' + 'a' + ': Trusted:    ' + '0.1' + '\n',
                    'INFO:root:' + '\n\nVariable ' + "'a'" + ' in module TestModule failed. Please check values.\n\n' +
                    'If you are confident that the newly calculated values are correct, comment out the old trusted ' +
                    'values for ' + "'TestModuleGlobals' in trusted_values_dict and copy the following code between" +
                    ' the ##### into trusted_values_dict. Make sure to fill out the TODO comment describing why the' +
                    ' values had to be changed. Then re-run test script.\n#####\n\n# Generated on: ' +
                    str(date.today()) + '\n# Reason for changing values: ' +
                    "TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'a'" + ": mpf('" + '0.0' + "')" + ', ' +
                    "'b'" + ": mpf('1.0')" + "}\n\n#####"
                ]
            )

            calculated_dict = {'a': mpf('0.0'), 'b': mpf('1.1')}
            trusted_dict = {'a': mpf('0.1'), 'b': mpf('1.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' +
                    '\nTestModule: ' + 'a' + ': Trusted:    ' + '0.1' + '\n',
                    'INFO:root:' + '\n\nVariable ' + "'a'" + ' in module TestModule failed. Please check values.\n\n' +
                    'If you are confident that the newly calculated values are correct, comment out the old trusted ' +
                    'values for ' + "'TestModuleGlobals' in trusted_values_dict and copy the following code between" +
                    ' the ##### into trusted_values_dict. Make sure to fill out the TODO comment describing why the' +
                    ' values had to be changed. Then re-run test script.\n#####\n\n# Generated on: ' +
                    str(date.today()) + '\n# Reason for changing values: ' +
                    "TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'a'" + ": mpf('" + '0.0' + "')" + ', ' +
                    "'b'" + ": mpf('1.1')" + "}\n\n#####"
                ]
            )

            calculated_dict = {'a': mpf('0.0'), 'b': mpf('1.1')}
            trusted_dict = {'a': mpf('0.0'), 'b': mpf('1.0')}

            with self.assertLogs(level='DEBUG') as logger:
                self.assertFalse(calc_error(mod, calculated_dict, trusted_dict))
            self.assertEqual(logger.output,
                [
                    'DEBUG:root:' + '\nTestModule: ' + 'a' + ': Calculated: ' + '0.0' +
                    '\nTestModule: ' + 'a' + ': Trusted:    ' + '0.0' + '\n',
                    'DEBUG:root:' + '\nTestModule: ' + 'b' + ': Calculated: ' + '1.1' +
                    '\nTestModule: ' + 'b' + ': Trusted:    ' + '1.0' + '\n',
                    'INFO:root:' + '\n\nVariable ' + "'b'" + ' in module TestModule failed. Please check values.\n\n' +
                    'If you are confident that the newly calculated values are correct, comment out the old trusted ' +
                    'values for ' + "'TestModuleGlobals' in trusted_values_dict and copy the following code between" +
                    ' the ##### into trusted_values_dict. Make sure to fill out the TODO comment describing why the' +
                    ' values had to be changed. Then re-run test script.\n#####\n\n# Generated on: ' +
                    str(date.today()) + '\n# Reason for changing values: ' +
                    "TODO\ntrusted_values_dict['TestModuleGlobals'] = {" + "'a'" + ": mpf('" + '0.0' + "')" + ', ' +
                    "'b'" + ": mpf('1.1')" + "}\n\n#####"
                ]
            )

            logging.info('\nAll calc_error tests passed.\n')

    def test_create_trusted_globals_dict(self):
        from UnitTesting.create_trusted_globals_dict import create_trusted_globals_dict
        from mpmath import mpf
        from UnitTesting.trusted_values_dict import trusted_values_dict

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

        from UnitTesting.calc_error import calc_error

        mod_dict = {'BrillLindquist': ['foo', 'bar']}
        first_times = {'BrillLindquist': False}
        self.assertTrue(calc_error('BrillLindquist', {'alphaCart': mpf('0.122483331574515176153136610247876'), 'betaCartU[0]': mpf('0.0'), 'betaCartU[1]': mpf('0.0'), 'betaCartU[2]': mpf('0.0'), 'BCartU[0]': mpf('0.0'), 'BCartU[1]': mpf('0.0'), 'BCartU[2]': mpf('0.0'), 'gammaCartDD[0][0]': mpf('66.6570391079152319165851690987334'), 'gammaCartDD[0][1]': mpf('0.0'), 'gammaCartDD[0][2]': mpf('0.0'), 'gammaCartDD[1][0]': mpf('0.0'), 'gammaCartDD[1][1]': mpf('66.6570391079152319165851690987334'), 'gammaCartDD[1][2]': mpf('0.0'), 'gammaCartDD[2][0]': mpf('0.0'), 'gammaCartDD[2][1]': mpf('0.0'), 'gammaCartDD[2][2]': mpf('66.6570391079152319165851690987334'), 'KCartDD[0][0]': mpf('0.0'), 'KCartDD[0][1]': mpf('0.0'), 'KCartDD[0][2]': mpf('0.0'), 'KCartDD[1][0]': mpf('0.0'), 'KCartDD[1][1]': mpf('0.0'), 'KCartDD[1][2]': mpf('0.0'), 'KCartDD[2][0]': mpf('0.0'), 'KCartDD[2][1]': mpf('0.0'), 'KCartDD[2][2]': mpf('0.0')}
                        , create_trusted_globals_dict(mod_dict, first_times)['BrillLindquist']))

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
        from UnitTesting.evaluate_globals import evaluate_globals
        from UnitTesting.functions_and_globals import functions_and_globals
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
        from UnitTesting.expand_variable_dict import expand_variable_dict

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
        from UnitTesting.functions_and_globals import functions_and_globals

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
        from UnitTesting.get_variable_dimension import get_variable_dimension

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
        from UnitTesting.is_first_time import is_first_time

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
        from UnitTesting.run_test import run_test

        mod_dict = {}
        with self.assertRaises(AssertionError):
            run_test(self, mod_dict, locals())

    def test_var_dict_to_value_dict(self):
        from UnitTesting.var_dict_to_value_dict import var_dict_to_value_dict
        from mpmath import mpf, sqrt, mp
        from random import random, seed
        from UnitTesting.trusted_values_dict import trusted_values_dict
        from UnitTesting.calc_error import calc_error
        from sympy.abc import x, y, z

        seed(trusted_values_dict['seed'])
        mp.dps = trusted_values_dict['precision']

        var_dict = dict()
        self.assertEqual(dict(), var_dict_to_value_dict(var_dict))

        first_val = mpf(sqrt(random()))
        second_val = mpf(sqrt(random()))
        third_val = mpf(sqrt(random()))

        var_dict = {'x': x}
        trusted_dict = {'x': first_val}
        calculated_dict = var_dict_to_value_dict(var_dict)
        self.assertTrue(calc_error('mod', calculated_dict, trusted_dict, output=False))

        var_dict = {'y': y}
        trusted_dict = {'y': first_val}
        calculated_dict = var_dict_to_value_dict(var_dict)
        self.assertTrue(calc_error('mod', calculated_dict, trusted_dict, output=False))

        var_dict = {'x': x, 'y': y}
        trusted_dict = {'x': first_val, 'y': second_val}
        calculated_dict = var_dict_to_value_dict(var_dict)
        self.assertTrue(calc_error('mod', calculated_dict, trusted_dict, output=False))

        var_dict = {'y': y, 'x': x}
        trusted_dict = {'y': second_val, 'x': first_val}
        calculated_dict = var_dict_to_value_dict(var_dict)
        self.assertTrue(calc_error('mod', calculated_dict, trusted_dict, output=False))

        var_dict = {'x': x, 'y': y, 'z': z}
        trusted_dict = {'x': first_val, 'y': second_val, 'z': third_val}
        calculated_dict = var_dict_to_value_dict(var_dict)
        self.assertTrue(calc_error('mod', calculated_dict, trusted_dict, output=False))

        var_dict = {'x+y': x+y}
        trusted_dict = {'x+y': first_val+second_val}
        calculated_dict = var_dict_to_value_dict(var_dict)
        self.assertTrue(calc_error('mod', calculated_dict, trusted_dict, output=False))

        var_dict = {'x/y': x/y, 'long_expression': (x+y)**2/(x*y)}
        trusted_dict = {'x/y': first_val/second_val, 'long_expression': (first_val+second_val)**2/(first_val*second_val)}
        calculated_dict = var_dict_to_value_dict(var_dict)
        self.assertTrue(calc_error('mod', calculated_dict, trusted_dict, output=False))

        var_dict = {'x-x': x-x}
        trusted_dict = {'x-x': mpf(0.0)}
        calculated_dict = var_dict_to_value_dict(var_dict)
        self.assertTrue(calc_error('mod', calculated_dict, trusted_dict, output=False))

        with self.assertRaises(AttributeError):
            var_dict = {'x': 0}
            calculated_dict = var_dict_to_value_dict(var_dict)

        var_dict = {'t1': x**2 + y**2, 't2': x**2/y, 't3': x**4, 't4':x**2 + x**2*y + x**2/y + x**2*z}
        trusted_dict = {'t1': first_val**2 + second_val**2, 't2': first_val**2/second_val, 't3': first_val**4,
                        't4': first_val**2 + first_val**2*second_val + first_val**2/second_val + first_val**2*third_val}
        calculated_dict = var_dict_to_value_dict(var_dict)
        self.assertTrue(calc_error('mod', calculated_dict, trusted_dict, output=False))

        from sympy import symbols

        M_SQRT1_2, M_PI = symbols('M_SQRT1_2 M_PI')

        var_dict = {'pi': M_PI, 'sqrt': M_SQRT1_2}
        trusted_dict = {'pi': mpf('3.14159265358979323846264338327933'),
                        'sqrt': mpf('0.707106781186547524400844362104785')}
        calculated_dict = var_dict_to_value_dict(var_dict)
        for var in trusted_dict:
            self.assertTrue(calc_error('Constants', calculated_dict, trusted_dict, output=False))

        logging.info('\nAll var_dict_to_value_dict tests passed\n')


# Necessary for unittest class to work properly
if __name__ == '__main__':
    unittest.main()
