from UnitTesting.create_test import create_test


def test_GiRaFFE_HO():

    module = 'GiRaFFE_HO.GiRaFFE_Higher_Order'

    module_name = 'GiRaFFE_HO'

    function_and_global_dict = {'GiRaFFE_Higher_Order()': ['uD', 'uU', 'gammaUU', 'gammadet', 'u0alpha', 'alpsqrtgam',
                                                             'Stilde_rhsD', 'AevolParen', 'PevolParenU', 'A_rhsD',
                                                             'psi6Phi_rhs']}

    create_test(module, module_name, function_and_global_dict)

def test_GiRaFFE_HO_v2():

    module = 'GiRaFFE_HO.GiRaFFE_Higher_Order_v2'

    module_name = 'GiRaFFE_HO_v2'

    function_and_global_dict = {'GiRaFFE_Higher_Order_v2()': ['gammaUU', 'gammadet', 'SevolParenUD', 'Stilde_rhsD',
                                                              'AevolParen', 'PevolParenU', 'A_rhsD', 'psi6Phi_rhs']}

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
