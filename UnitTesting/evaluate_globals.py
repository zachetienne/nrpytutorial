import logging
from importlib import import_module

# evaluate_globals sequentially imports self.module, runs self.initialization_string, calls self.function on
# self.module, and gets an expression for each global in self.global_list. It returns a dictionary whose keys represent
# the names of the globals and whose values represent the expressions calculated for each respective global.

# Called by run_test

# Uses self.module, self.module_name, self.initialization_string, self.function, self.global_list


def evaluate_globals(self):

    logging.debug(' Importing ' + self.module + '...')

    # Try to import self.module
    try:
        imported_module = import_module(self.module)
    # If user supplied an incorrect module, error
    except ImportError:
        logging.error(" Attribute 'module' for " + self.module_name + " does not exist as a module. This attribute "
                      "should be what you would type if you were importing 'module' in your own file.\n")
        self.assertTrue(False)
    except ValueError:
        logging.error(" Attribute 'module' for " + self.module_name + " is empty -- it must have a value.")
        self.assertTrue(False)

    logging.debug(' ...Success: Imported module.')

    logging.debug(' Executing initialization_string...')

    logging.debug('initialization_string: ' + self.initialization_string)

    # Execute self.initialization_string
    exec(self.initialization_string)

    logging.debug(' ...Successfully executed.')

    # If the user supplied a function, add it to string_exec
    if self.function != '':
        string_exec = self.module_name + '.' + self.function + '\n'
    else:
        string_exec = ''

    # Add each global to string_exec
    for glob in self.global_list:
        string_exec += glob + '=' + self.module_name + '.' + glob + '\n'

    # Initializing location that result will be stored in
    var_dict = {}

    logging.debug(' Executing function call and global assignment...')

    # Execute string_exec with imported_module as environment and store result in var_dict
    exec(string_exec, {self.module_name: imported_module}, var_dict)

    logging.debug(' ...Successfully executed.')

    return var_dict
