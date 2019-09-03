# Deprecated: See jupyter notebook

# NRPy Unit Testing: Functions In-Depth 
If you have not already, please read through the Jupyter notebook
tutorial for unit testing `Tutorial-UnitTesting.ipynb` This will contain
in-depth information on all functions used for unit testing, not a
high-level user tutorial. For examples, the default module that will be
used is `UnitTesting/Test_UnitTesting/test_module.py`.

## failed_tests:
`failed_tests.txt` is a simple text file that keeps track of which tests
failed. Line 1 is by default 'Failures: '. The subsequent lines tell the
user which test functions in which test files failed in the following
format: `[test file path]: [test function]`

Example:

Say that the test function `test_module_for_testing_no_gamma()` failed.
Then we'd expect `failed_tests.txt` to be the following:

```
Failures:

UnitTesting/Test_UnitTesting/test_module.py: test_module_for_testing_no_gamma
```

## standard_constants:
`standard_constants.py` stores test-wide information that the user can
modify to impact the numerical result for their globals. It currently
only has one field, `precision`, which determines how precise the values
for the globals are. It is by default set to `30`, which we've
determined to be a reasonable amount. This file has the ability to be
expanded upon in the future, but it is currently minimal.

## run_NRPy_UnitTests:
`run_NRPy_UnitTests.sh` is a bash script that acts as the hub for
running tests -- it is where the user specifies the tests they'd like to
be run. It keeps track of which tests failed by interacting with
[failed_tests](#failed_tests), giving the user easily readable output
from the file. It also has the option to automatically rerun the tests
that failed in `DEBUG` mode if the boolean `rerun_if_fail` is `true`.

The script is run with the following syntax:

```
./UnitTesting/run_NRPy_UnitTests.sh [python interpreter]
```

This of course assumes that the user is in the nrpy directory; the user
simply has to specify the path from their current directory to the bash
file.
 
Examples of `python interpreter` are `python` and `python3`.

The script first lets the user know if they forgot to pass a python
interpreter. Then if they didn't, it prints some baseline information
about Python variables: `PYTHONPATH`, `PYTHONEXEC`, and `PYTHONEXEC`
version.

`failed_tests.txt` is then overwritten with the default information.
This makes it so that each subsequent test call has a unique list of the
tests that passed; it wouldn't make sense to store this information.

The user can then change the boolean `rerun_if_fail` if need be. Next,
the user can add tests using the `add_test` function. The syntax is as
follows: `add_test [path to test file]`

Example:

```
add_test UnitTesting/Test_UnitTesting/test_module.py
```

Finally, the bash script will read any failures from `failed_tests.txt`
and, if `rerun_if_fail` is `true`, rerun those tests. It lastly prints
which tests failed in the same format as `failed_tests.txt`, and if no
tests failed, a success message.

## create_test:
create_test is a function that takes the following user-supplied
information: a module to test `module`, the name of the module
`module_name`, and a dictionary whose keys are functions and whose
values are lists of globals `function_and_global_dict`. It uses this
information to generate a test file that is automatically run as a bash
script; this test file does all the heavy lifting in calling the
function, getting expressions for all the globals, evaluating the
expressions to numerical values, and storing the values in the proper
trusted_values_dict.

create_test additionally takes optional arguments `logging_level` and
`initialization_string_dict`, which respectively determine the desired
level of output (think verbosity) and run some python code prior calling
the specified function. Usage is as following:

```
module = 'BSSN.BrillLindquist'

module_name = 'BrillLindquist'

function_and_global_dict = {'BrillLindquist(ComputeADMGlobalsOnly = True)': ['alphaCart', 'betaCartU', 'BCartU', 'gammaCartDD', 'KCartDD']}

create_test(module, module_name, function_and_global_dict)
```

The way to think of this is that the module to be tested is
BSSN.BrillLindquist. The module_name is how you refer to this module --
it's a bit arbitrary, so whether you prefer BrillLindquist or bl, it
won't change the computation. The function_and_global_dict contains
entry 'BrillLindquist(ComputeADMGlobalsOnly = True)', which is the
function that gets called in the module. It's value in the dictionary is
a list of globals that get created when this function gets called.

Now let's add the optional arguments into the same example:

```
module = 'BSSN.BrillLindquist'

module_name = 'BrillLindquist'

function_and_global_dict = {'BrillLindquist(ComputeADMGlobalsOnly = True)': ['alphaCart', 'betaCartU', 'BCartU', 'gammaCartDD', 'KCartDD']}

logging_level = 'DEBUG'

initialization_string_dict = {'BrillLindquist(ComputeADMGlobalsOnly = True)': 'print("example")\nprint("Hello world!")'}

create_test(module, module_name, function_and_global_dict, logging_level=logging_level, initialization_string_dict=initialization_string_dict)
```

Now when create_test runs, the user will be given much more output due
to the logging_level; additionally, the user-specified print will occur
due to initialization_string_dict.

You may now be wondering why we use dictionaries to store this data
instead of simply having separate variables `function`, `global_list`,
and `initialization_string`. This is where some of the power of this
testing method lies: we can test multiple functions and their globals
with ease! In other words, function_and_global_dict can contain multiple
entries, each a specific function call with its own associated list of
globals. Since not every function being tested must have an associated
initialization_string, we make an entry for each function optional. An
example is as follows:

```
module = 'BSSN.BrillLindquist'

module_name = 'BrillLindquist'

function_and_global_dict = {'BrillLindquist(ComputeADMGlobalsOnly = True)': ['alphaCart', 'betaCartU', 'BCartU', 'gammaCartDD', 'KCartDD'],
                            'BrillLindquist(ComputeADMGlobalsOnly = False)': ['alphaCart', 'betaCartU', 'BCartU', 'gammaCartDD', 'KCartDD']}

logging_level = 'DEBUG'

initialization_string_dict = {'BrillLindquist(ComputeADMGlobalsOnly = True)': 'print("example")\nprint("Hello world!")'}

create_test(module, module_name, function_and_global_dict, logging_level=logging_level, initialization_string_dict=initialization_string_dict)
```

Both instances will be called separately, with their own globals. The
print statements will only be called in the first function, since there
is no associated initialization_string for the second function as well.

An important note when using `create_test` is that all arguments are
**strings**. This includes the module, module_name, function, each
global in the list of globals, logging level, and initialization_string.
The reason for making these fields strings is that when setting
module_name, for example, there doesn't exist anything in Python with
the name BrillLindquist. So, we wrap it in a string. This is true of
every input. Be careful with the dicts and lists, however: their
arguments are strings, they aren't themselves strings.

## setup_trusted_values_dict:
`setup_trusted_values_dict` takes in a path to a test directory `path`,
and checks whether or not a `trusted_values_dict.py` exists in the test
directory. If it does exist, the function does nothing. If it doesn't
exist, `setup_trusted_values_dict` creates the file
`trusted_values_dict.py` in the test directory. In then writes the
following default code into the file:

```
from mpmath import mpf, mp, mpc
from UnitTesting.standard_constants import precision

mp.dps = precision
trusted_values_dict = {}

```

The default code allows the unit test to properly interact with and
write to the file.

## run_test:
`run_test` acts as the hub for an individual unit test. It takes in all
the module-wide information, and goes through all the steps of
determining whether the test passed or failed by calling many
sub-functions. A fundamentally important part of `run_test` is the
notion of `self`; `self` stores a test's information (i.e. `module`,
`module_name`, etc.) to be able to easily pass information to and make
assertions in sub-functions. When `self` is referenced, simply think
"information storage".

`run_test` begins by importing the `trusted_values_dict` of the current
module being tested; since `setup_trusted_values_dict` is called before
`run_test`, we know it exists.

`run_test` then determines if the current function/module is being done
for the first time based off the existence of the proper entry in
`trusted_values_dict`. 

[evaluate_globals](#evaluate_globals) is then run in order to generate
the SymPy expressions for each global being tested. 

Next,
[cse_simplify_and_evaluate_sympy_expressions](#cse_simplify_and_evaluate_sympy_expressions)
is 



## evaluate_globals:

## cse_simplify_and_evaluate_sympy_expressions:

## create_dict_string:

## first_time_print:

## calc_error:
