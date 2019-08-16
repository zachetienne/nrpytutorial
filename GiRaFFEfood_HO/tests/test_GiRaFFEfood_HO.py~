from UnitTesting.create_test import create_test


def test_GiRaFFEfood_HO():

    module = 'GiRaFFEfood_HO.GiRaFFEfood_HO'

    module_name = 'GiRaFFEfood_HO'

    function_and_global_dict = {'GiRaFFEfood_HO()': ['AD', 'ValenciavU']}

    create_test(module, module_name, function_and_global_dict)

def test_GiRaFFEfood_HO_Aligned_Rotator():

    module = 'GiRaFFEfood_HO.GiRaFFEfood_HO_Aligned_Rotator'

    module_name = 'GiRaFFEfood_HO_Aligned_Rotator'

    function_and_global_dict = {'GiRaFFEfood_HO_Aligned_Rotator()': ['AD', 'ValenciavU']}

    create_test(module, module_name, function_and_global_dict)

def test_GiRaFFEfood_HO_1D_tests():

    module = 'GiRaFFEfood_HO.GiRaFFEfood_HO_1D_tests'

    module_name = 'GiRaFFEfood_HO_1D_tests'

    function_and_global_dict = {'GiRaFFEfood_HO_1D_tests()': ['AleftD', 'AcenterD', 'ArightD', 'ValenciavleftU', 'ValenciavcenterU', 'ValenciavrightU']}

    create_test(module, module_name, function_and_global_dict)


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
