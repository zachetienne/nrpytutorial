# Necessary imports for unit testing framework
import unittest
import logging
from UnitTesting.run_test import run_test
from UnitTesting.functions_and_globals import functions_and_globals
from UnitTesting.RepeatedTimer import RepeatedTimer

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

    # Testing globals for reference_metric
    def test_reference_metric_globals(self):

        import NRPy_param_funcs as par

        # Globals used by all coordinate systems
        basic_global_list = ['UnitVectors', 'ReU', 'ReDD', 'ghatDD', 'ghatUU', 'detgammahat', 'detgammahatdD',
                             'detgammahatdDD', 'ReUdD', 'ReUdDD', 'ReDDdD', 'ReDDdDD', 'ghatDDdD', 'ghatDDdDD',
                             'GammahatUDD', 'GammahatUDDdD']

        # Dictionary of coordinate systems and their additional globals
        coord_dict = {'Spherical': ('True', ['xxmin', 'xxmax']), 'SinhSpherical': ('True', []),
                      'SinhSphericalv2': ('True', []), 'NobleSphericalThetaOptionOne': ('False', []),
                      'NobleSphericalThetaOptionTwo': ('False', []), 'Cylindrical': ('True', []),
                      'SinhCylindrical': ('True', []), 'SinhCylindricalv2': ('True', []), 'SymTP': ('True', []),
                      'SinhSymTP': ('True', []), 'Cartesian': ('True', [])}

        # For each module and its respective boolean and additional globals
        for coord, bool_global_tuple in coord_dict.items():

            # Import module based on string
            exec('import reference_metric as rfm_' + coord)

            # Set reference metric
            par.set_parval_from_str("reference_metric::CoordSystem", coord)

            # Create mod_dict
            global_list = bool_global_tuple[1] + basic_global_list
            mod_dict = {('rfm_' + coord): functions_and_globals(['reference_metric(' + bool_global_tuple[0] + ')'],
                                                                global_list)}

            # Run test and delete old entry
            run_test(self, mod_dict, locals())
            del locals()['rfm_' + coord]


# Necessary for unittest class to work properly
if __name__ == '__main__':
    unittest.main()
