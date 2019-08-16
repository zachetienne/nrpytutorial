from UnitTesting.create_test import create_test


def test_Spherical():

    module = 'reference_metric'

    module_name = 'rfm_Spherical'

    function_and_global_dict = {'reference_metric(True)': ['xxmin', 'xxmax', 'UnitVectors', 'ReU', 'ReDD', 'ghatDD', 'ghatUU', 'detgammahat',
                       'detgammahatdD', 'detgammahatdDD', 'ReUdD', 'ReUdDD', 'ReDDdD', 'ReDDdDD', 'ghatDDdD',
                       'ghatDDdDD', 'GammahatUDD', 'GammahatUDDdD', 'Cart_to_xx','xxCart','xxSph','scalefactor_orthog']}

    initialization_string_dict = {'reference_metric(True)': '''
import NRPy_param_funcs as par
par.set_parval_from_str("reference_metric::CoordSystem", "Spherical")
'''}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)

def test_SinhSpherical():

    module = 'reference_metric'

    module_name = 'rfm_SinhSpherical'

    function_and_global_dict = {'reference_metric(True)': ['xxmin', 'xxmax', 'UnitVectors', 'ReU', 'ReDD', 'ghatDD', 'ghatUU', 'detgammahat',
                       'detgammahatdD', 'detgammahatdDD', 'ReUdD', 'ReUdDD', 'ReDDdD', 'ReDDdDD', 'ghatDDdD',
                       'ghatDDdDD', 'GammahatUDD', 'GammahatUDDdD', 'Cart_to_xx','xxCart','xxSph','scalefactor_orthog']}

    initialization_string_dict = {'reference_metric(True)': '''
import NRPy_param_funcs as par
par.set_parval_from_str("reference_metric::CoordSystem", "SinhSpherical")
'''}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)


def test_SinhSphericalv2():

    module = 'reference_metric'

    module_name = 'rfm_SinhSphericalv2'

    function_and_global_dict = {'reference_metric(True)': ['xxmin', 'xxmax', 'UnitVectors', 'ReU', 'ReDD', 'ghatDD', 'ghatUU', 'detgammahat',
                       'detgammahatdD', 'detgammahatdDD', 'ReUdD', 'ReUdDD', 'ReDDdD', 'ReDDdDD', 'ghatDDdD',
                       'ghatDDdDD', 'GammahatUDD', 'GammahatUDDdD', 'Cart_to_xx','xxCart','xxSph','scalefactor_orthog']}

    initialization_string_dict = {'reference_metric(True)': '''
import NRPy_param_funcs as par
par.set_parval_from_str("reference_metric::CoordSystem", "SinhSphericalv2")
'''}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)

def test_NobleSphericalThetaOptionOne():

    module = 'reference_metric'

    module_name = 'rfm_NobleSphericalThetaOptionOne'

    function_and_global_dict = {'reference_metric(False)': ['UnitVectors', 'ReU', 'ReDD', 'ghatDD', 'ghatUU', 'detgammahat',
                       'detgammahatdD', 'detgammahatdDD', 'ReUdD', 'ReUdDD', 'ReDDdD', 'ReDDdDD', 'ghatDDdD',
                       'ghatDDdDD', 'GammahatUDD', 'GammahatUDDdD', 'Cart_to_xx','xxCart','xxSph','scalefactor_orthog']}

    initialization_string_dict = {'reference_metric(False)': '''
import NRPy_param_funcs as par
par.set_parval_from_str("reference_metric::CoordSystem", "NobleSphericalThetaOptionOne")
'''}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)

def test_NobleSphericalThetaOptionTwo():

    module = 'reference_metric'

    module_name = 'rfm_NobleSphericalThetaOptionTwo'

    function_and_global_dict = {'reference_metric(False)': ['UnitVectors', 'ReU', 'ReDD', 'ghatDD', 'ghatUU', 'detgammahat',
                       'detgammahatdD', 'detgammahatdDD', 'ReUdD', 'ReUdDD', 'ReDDdD', 'ReDDdDD', 'ghatDDdD',
                       'ghatDDdDD', 'GammahatUDD', 'GammahatUDDdD', 'Cart_to_xx','xxCart','xxSph','scalefactor_orthog']}

    initialization_string_dict = {'reference_metric(False)': '''
import NRPy_param_funcs as par
par.set_parval_from_str("reference_metric::CoordSystem", "NobleSphericalThetaOptionTwo")
'''}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)

def test_Cylindrical():

    module = 'reference_metric'

    module_name = 'rfm_Cylindrical'

    function_and_global_dict = {'reference_metric(True)': ['xxmin', 'xxmax', 'UnitVectors', 'ReU', 'ReDD', 'ghatDD', 'ghatUU', 'detgammahat',
                       'detgammahatdD', 'detgammahatdDD', 'ReUdD', 'ReUdDD', 'ReDDdD', 'ReDDdDD', 'ghatDDdD',
                       'ghatDDdDD', 'GammahatUDD', 'GammahatUDDdD', 'Cart_to_xx','xxCart','xxSph','scalefactor_orthog']}

    initialization_string_dict = {'reference_metric(True)': '''
import NRPy_param_funcs as par
par.set_parval_from_str("reference_metric::CoordSystem", "Cylindrical")
'''}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)

def test_SinhCylindrical():

    module = 'reference_metric'

    module_name = 'rfm_SinhCylindrical'

    function_and_global_dict = {'reference_metric(True)': ['xxmin', 'xxmax', 'UnitVectors', 'ReU', 'ReDD', 'ghatDD', 'ghatUU', 'detgammahat',
                       'detgammahatdD', 'detgammahatdDD', 'ReUdD', 'ReUdDD', 'ReDDdD', 'ReDDdDD', 'ghatDDdD',
                       'ghatDDdDD', 'GammahatUDD', 'GammahatUDDdD', 'Cart_to_xx','xxCart','xxSph','scalefactor_orthog']}

    initialization_string_dict = {'reference_metric(True)': '''
import NRPy_param_funcs as par
par.set_parval_from_str("reference_metric::CoordSystem", "SinhCylindrical")
'''}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)

def test_SinhCylindricalv2():

    module = 'reference_metric'

    module_name = 'rfm_SinhCylindricalv2'

    function_and_global_dict = {'reference_metric(True)': ['xxmin', 'xxmax', 'UnitVectors', 'ReU', 'ReDD', 'ghatDD', 'ghatUU', 'detgammahat',
                       'detgammahatdD', 'detgammahatdDD', 'ReUdD', 'ReUdDD', 'ReDDdD', 'ReDDdDD', 'ghatDDdD',
                       'ghatDDdDD', 'GammahatUDD', 'GammahatUDDdD', 'Cart_to_xx','xxCart','xxSph','scalefactor_orthog']}

    initialization_string_dict = {'reference_metric(True)': '''
import NRPy_param_funcs as par
par.set_parval_from_str("reference_metric::CoordSystem", "SinhCylindricalv2")
'''}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)

def test_SymTP():

    module = 'reference_metric'

    module_name = 'rfm_SymTP'

    function_and_global_dict = {'reference_metric(True)': ['xxmin', 'xxmax', 'UnitVectors', 'ReU', 'ReDD', 'ghatDD', 'ghatUU', 'detgammahat',
                       'detgammahatdD', 'detgammahatdDD', 'ReUdD', 'ReUdDD', 'ReDDdD', 'ReDDdDD', 'ghatDDdD',
                       'ghatDDdDD', 'GammahatUDD', 'GammahatUDDdD', 'Cart_to_xx','xxCart','xxSph','scalefactor_orthog']}

    initialization_string_dict = {'reference_metric(True)': '''
import NRPy_param_funcs as par
par.set_parval_from_str("reference_metric::CoordSystem", "SymTP")
'''}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)

def test_SinhSymTP():

    module = 'reference_metric'

    module_name = 'rfm_SinhSymTP'

    function_and_global_dict = {'reference_metric(True)': ['xxmin', 'xxmax', 'UnitVectors', 'ReU', 'ReDD', 'ghatDD', 'ghatUU', 'detgammahat',
                       'detgammahatdD', 'detgammahatdDD', 'ReUdD', 'ReUdDD', 'ReDDdD', 'ReDDdDD', 'ghatDDdD',
                       'ghatDDdDD', 'GammahatUDD', 'GammahatUDDdD', 'Cart_to_xx','xxCart','xxSph','scalefactor_orthog']}

    initialization_string_dict = {'reference_metric(True)': '''
import NRPy_param_funcs as par
par.set_parval_from_str("reference_metric::CoordSystem", "SinhSymTP")
'''}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)

def test_Cartesian():

    module = 'reference_metric'

    module_name = 'rfm_Cartesian'

    function_and_global_dict = {'reference_metric(True)': ['UnitVectors', 'ReU', 'ReDD', 'ghatDD', 'ghatUU', 'detgammahat',
                       'detgammahatdD', 'detgammahatdDD', 'ReUdD', 'ReUdDD', 'ReDDdD', 'ReDDdDD', 'ghatDDdD',
                       'ghatDDdDD', 'GammahatUDD', 'GammahatUDDdD', 'Cart_to_xx','xxCart','xxSph','scalefactor_orthog']}

    initialization_string_dict = {'reference_metric(True)': '''
import NRPy_param_funcs as par
par.set_parval_from_str("reference_metric::CoordSystem", "Cartesian")
'''}

    create_test(module, module_name, function_and_global_dict, initialization_string_dict=initialization_string_dict)


if __name__ == '__main__':
    import sys

    if len(sys.argv) <= 3:
        failed_functions = []
        for fun in dir():
            if fun[0:5] == 'test_':
                print('\nTesting ' + str(fun) + '...\n')
                try:
                    exec(fun + '()')
                except SystemExit:
                    failed_functions.append(fun)

        if failed_functions != []:
            import sys, os
            with open(os.path.join('UnitTesting', 'failed_tests.txt'), 'a') as file:
                for function in failed_functions:
                    file.write(sys.argv[0] + ': ' + str(function) + '\n')
            exit(1)

    else:
        globals()[sys.argv[4]]()
