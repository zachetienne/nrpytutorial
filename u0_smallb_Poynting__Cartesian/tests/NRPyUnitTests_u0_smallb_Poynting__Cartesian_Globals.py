# Necessary imports for unit testing framework
import unittest
import logging
from UnitTesting.run_test import run_test
from UnitTesting.functions_and_globals import functions_and_globals
from UnitTesting.RepeatedTimer import RepeatedTimer
from trusted_values_dict import trusted_values_dict

# TODO: Change level based on desired amount of output.
# ERROR -> Outputs minimal information -- only when there's an error
# INFO -> Outputs when starting and finishing a module, as well as everything in ERROR
# DEBUG -> Displays all pairs of values being compared, as well as everything in INFO
# NOTSET -> Displays symbolic dictionary for all modules, as well as everything in DEBUG
logging.basicConfig(level=logging.INFO)


# https://stackoverflow.com/questions/3393612/run-certain-code-every-n-seconds/13151299
# Creates a threaded timer object that prints to the console every 5 minutes
Timer = RepeatedTimer(300, logging.info, "\nPrinting every 5 minutes to prevent timeouts.\n")


# Python unittest class
class TestGlobals(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        Timer.stop()

    # Testing globals
    def test_globals(self):

        # TODO: Import modules to be tested
        # Note: Even though it says the modules are unused, these imports are vital for run_test to work properly.
        # Their information gets passed into run_test through locals()
        import u0_smallb_Poynting__Cartesian.u0_smallb_Poynting__Cartesian as u0sbPoyn

        # TODO: Create lists of globals to calculate
        global_list = ['g4DD', 'u0', 'uD', 'uBcontraction', 'uU', 'smallb4U', 'g4UU',
                       'smallb4D', 'smallb2etk', 'PoynSU']

        # TODO: Create Module dictionary based on imported modules, functions to initialize the modules, and globals
        # Note that the name of the modules in mod_dict MUST have the same name as the imported module.
        # Example: If you say 'import My_Modules.Module1 as M1', then mod_dict should have the entry 'M1' as a string.
        mod_dict = {'u0sbPoyn': functions_and_globals(['compute_u0_smallb_Poynting__Cartesian()'], global_list)}

        # TODO: Call run_test with arguments (self, mod_dict, locals())
        run_test(self, mod_dict, trusted_values_dict, locals())


# Necessary for unittest class to work properly
if __name__ == '__main__':
    unittest.main()
