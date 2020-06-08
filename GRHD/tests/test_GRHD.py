from UnitTesting.create_test import create_test


def test_equations():

    module = 'GRHD.equations'

    module_name = 'GRHD'

    function_and_global_dict = {'generate_everything_for_UnitTesting()': ['h','T4UU','T4UD','sqrtgammaDET','rho_star','tau_tilde','S_tildeD','vU','rho_star_fluxU','tau_tilde_fluxU','S_tilde_fluxUD',
                                                                          's_source_term','g4DD_zerotimederiv_dD','S_tilde_source_termD','rescaledValenciavU','u4U_ito_ValenciavU','rescaledvU','u4U_ito_vU']}

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
