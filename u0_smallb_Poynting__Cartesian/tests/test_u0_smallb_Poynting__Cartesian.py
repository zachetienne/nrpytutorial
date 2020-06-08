from UnitTesting.create_test import create_test


def test_u0_smallb_Poynting__Cartesian():

    module = 'u0_smallb_Poynting__Cartesian.u0_smallb_Poynting__Cartesian'

    module_name = 'u0sbPoyn'

    function_and_global_dict = {'compute_u0_smallb_Poynting__Cartesian()': ['u0', 'uD', 'uBcontraction', 'uU',
                                                                            'smallb4U', 'smallb4D',
                                                                            'smallb2etk', 'PoynSU']}

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
            sys.exit(1)

    else:
        globals()[sys.argv[4]]()

