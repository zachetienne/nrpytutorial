## NRPy Unit Testing Globals: An In-Depth Guide

### Motivation:

What is the purpose of unit testing, and why should you do it? To begin thinking about that, consider what subtleties 
can occur with your code that are almost unnoticeable to the eye, but wind up giving you a very incorrect result. You
could make a small optimization, and observe that nothing changes about your result. However, maybe the optimization 
you made only works on Python 3 and not Python 2, or it changes a value by some tiny amount -- too small to be obviously
noticeable, but enough to make a difference.

This is where unit testing comes in. By initially calculating values for the globals of your modules in a **trusted**
version of your code and storing those values in a dictionary, you can then easily check if something stopped working
by comparing your newly calculated values to the ones you've stored. On the frontend, there are four modules 
essential to get your unit tests up and running: `trusted_values_dict`, `functions_and_globals`, `run_test`, and your 
testing module (which we'll simply reference as `Your_Test_File`). The usage of each of these modules is outlined in 
the  **Interactive Modules** section. There is also some important prerequisite knowledge that may be helpful to grasp
before beginning your testing. There are many functions at play in the backend as well, all of which will 
be described in detail below in the **Functions** section. Mastery of these functions may not be essential to get your 
tests up-and-running, but some basic understanding of these modules with undoubtedly streamline the testing process and 
how to potentially create your own, different types of tests.

An important caveat is that the unit testing does not test the **correctness** of your code or your variables. The 
unit tests function as a protective measure to ensure that nothing was broken; it gets its values by running _your_ 
code, so if something starts out incorrect, it will be stored as incorrect in the system. There are measures against 
this, but it relies on the user's knowledge of what versions of their code are correct.

### Important Prerequisite Knowledge:

#### Unit Testing:
TODO

#### Logging:
`logging` is a Python module with the purpose of creating meaningful and flexible event logging. For the purpose of
NRPy Unit Testing, we use `logging` to display information, errors, debug statements, etc. to the console in real-time.
Unlike simple printing, which requires sloppy if-statements to determine whether to print or not and only prints
once a function call is complete, `logging` has a built-in `level` parameter that makes it easy to set the level of
output for all functions and an ability to print output during a function call. The following levels and their 
descriptions are what we'll be using in the NRPy Unit Testing context.

###### Levels:
````
ERROR -> Outputs minimal information about processes. Only prints when necessary due to an error.
INFO -> Outputs when a module begins and finishes being tested, as well as everything above.
DEBUG -> Outputs all pairs of calculated and trusted values being compared, as well as everything above.
NOTSET -> Outputs symbolic expressions for each variable, as well as everything above.
````

The recommended level for `logging` is `INFO`, as it gives the user an understanding of how long each module takes while
not being overbearing. The logging level should be set at the top of the user's Unit Testing file with the following 
command:

````
logging.basicConfig(level=logging.INFO)
```` 

The level can be changed by simply changing `INFO` to any of `ERROR`, `INFO`, etc. 

#### Your_Test_File.py:
TODO

### Interactive Modules:

####trusted_values_dict:
`trusted_values_dict` acts as the storage hub for NRPy unit tests. It is a module that stores trusted values of the 
globals that you calculate in an easily accessible dictionary. It also stores the precision `precision` that is used in 
comparing your trusted values and calculated values, as well as a seed `seed` that is used by some functions below. A 
good default value for `precision` is `30`, and a standard `seed` is `1234`. 


#### functions_and_globals:
`functions_and_globals` is a simple function that takes in a list of functions and a list of globals, and returns a
dictionary containing those functions and globals in a format that is readable by the other functions. The purpose of
this function is to be called to create the value for your module in your module dictionary. Examples are shown below:

````
my_function_list = ['function_1()', 'function_2(True)', 'function_3(1, 2, True)']
my_global_list = ['x', 'y', 'z']
functions_and_globals(function_list, global_list) -> {'function_list': ['function_1()', 
'function_2(True)', 'function_3(1, 2, True)'], 'global_list': ['x', 'y', 'z']}
````

Once in this format, the resulting dictionary can easily be passed into your module dictionary without having to 
concern yourself with formatting the dictionary correctly. Example usage is as follows:

````
mod_dict = {'myMod': functions_and_globals(function_list, global_list)}
````

`functions_and_globals` throws an `assertion_error` if either argument passed in isn't a list, or if any entry 
in each list passed in isn't a string.

#### run_test basic:
`run_test` is the culmination of all the functions outlined below. It does all the heavy lifting in calculating values,
comparing them to trusted values, throwing errors if necessary, printing the desired output, etc. `run_test` takes in 
`self`, which simply allows it to use the assert functions of `unittest`, `mod_dict`, which is the user-created module 
dictionary containing the modules they're testing, the necessary functions to run for each module, and the globals to 
evaluate and compare for each module, and `locs` (almost always a simple call to the built-in function 
`locals`), which allows the current local variables to be passed into `run_test`.

At the most basic level, an understanding of exactly what `run_test` does isn't necessary to start creating your own
tests, and as such a simple usage tutorial will be given here, while a more in-depth tutorial will be given 
[below](#run_test) in the **Functions** section. Assuming you've created your module dictionary `mod_dict` 
correctly and imported its necessary modules, running `run_test` is as simple as calling the following:

````
run_test(self, mod_dict, locals())
````

### Functions:

#### calc_error:
`calc_error` takes in a module `mod`, a calculated dictionary `calculated_dict`, a trusted dictionary `trusted_dict`, 
and a symbolic dictionary `symbolic_dict`. It computes the error between each respective value in `result_dict` and 
`trusted_dict`, and returns `False` if the error between any values in the two dictionaries is too large, `True` if
all values in both dictionaries are close enough in value. It prints output to the screen according to the level of
`logging`. A more in-depth explanation is given after the following example:

````
mod = 'myMod'
calculated_dict = {'alpha': mpf('0.12759172659175684561176895'), 'beta': mpf('60.561850917509181234786')}
trusted_dict = {'alpha': mpf('0.12759172659175684561176904'), 'beta': mpf('60.5618509175149189740873')}
symbolic_dict = {'alpha': sqrt(x**2 + y**2), 'beta': x**3/y**5}
calc_error(mod, calculated_dict, trusted_dict, symbolic_dict) -> False

Console output:

ERROR:root:
Variable beta in module myMod failed. Please check values.
If you are confident that the newly calculated values are correct, comment out the old trusted values 
for 'myModGlobals' in trustedValuesDict and copy the following code between the ##### into 
trustedValuesDict. Make sure to fill out the TODO comment describing why the values had to be changed. 
Then re-run test script.
#####
# Generated on: 2019-06-14 14:03:20.771968
# Reason for changing values: TODO
trustedValuesDict['myModGlobals'] = {'alpha': mpf('0.127591726591756854380932395542914'), 
'beta': mpf('60.5618509175091830343262699898332')}
#####

````

As we can see above, variable `alpha` was within the precision value, and as such was found to be correct. However, 
variable `beta` was too far off to be considered correct, so `calc_error` returned `False` and output to the console
with instructions.

Now for a more in-depth look into what calc_error does, how it does it, and what it outputs to the console in what case.
First, `calc_error` takes `calculated_dict` and `trusted_dict` and creates their sets (`calculated_set` and 
`trusted_set`, respectively). This makes it super easy to compare the keys in each dictionary. If the sets are not 
equal, that means at least one of the sets has a value that's not in the other set. For example, this could mean that 
`calculated_dict` has variables `alpha` and `beta`, while `trusted_dict` only has variable `alpha`. Note that there's 
nothing special about creating sets out of the dictionaries, other than for ease of comparison. Then, the fact that 
they differ is printed to the screen under the `Error` level of `logging`, along with the specific variables in each
set that aren't in the other set:

````
ERROR:root:
	myMod: Calculated dictionary and trusted dictionary have different variables.
ERROR:root:
    Calculated Dictionary variables not in Trusted Dictionary: 
	{'beta'}
````

Since there are no variables in `trusted_dict` that _aren't_ in `calculated_dict`, there is no output of this type. 
However, if both dictionaries had variables that weren't in the other dictionary, we'd get output from both of them.
After this printing occurs, `calc_error` returns `False`; it doesn't bother comparing the variables that were in both
dictionaries, since it'll have to be rerun anyway and that would be a waste of processing power.

In the case that both dictionaries have the same variables, `calc_error` then compares the values for each variable.
If the `logging` level is `Debug` or lower, both the calculated and trusted values will be printed to the console. 
Additionally, if the `logging` level is `Notset` (which is strongly not recommended), the symbolic expression will be
printed as well. Finally, `calc_error` will calculate the error. 

The error is calculated by looking at how many digits of the two values being compared are identical. It does this by 
taking the `log10` of the error between the trusted value and the calculated value divided by the trusted value.
A keen eye may notice that this will give us a division by zero in the case where the trusted value is zero; in this
we simply look at distance from zero of the calculated value.

After we get this number, we compare it with our `precision` value as set in `trusted_values_dict`. Since it'd be 
unreasonable to expect the calculated and trusted values to be identical, we simply say that the the calculated and 
trusted values must agree to at least `precision/2` digits. If the error meets this constraint, we easily move on
to the next variable to be compared. However, if the error is too large, we print to the console at the `Error`
level of `logging` that two values differed by too much, as well as a potentially correct new entry for 
`trusted_values_dict` if you're confident that the newly calculated values are correct. `calc_error` then returns 
`False`. In the case that the error for every value meets the constraint, `True` is returned and `calc_error` is
finished.

#### create_trusted_globals_dict:
`create_trusted_globals_dict` takes in a module dictionary `mod_dict` and a boolean list `first_times`  representing 
if it's the first time running each module in `mod_dict`. For each module, if `first_time` is `True`, then that module 
gets assigned an empty dictionary. This ensures that, once the result of `create_trusted_globals_dict` is passed into
`calc_error`, that an error will arise. If `first_time` is False, then that module gets assigned its entry in 
`trusted_values_dict` along the same naming convention we defined above in `trusted_values_dict`'s tutorial.. If 
`mod_dict` and `first_times` have different lengths, then an `AssertionError` is thrown. If a module has `first_time` 
as `False` but it doesn't have an entry in `trusted_values_dict` (this should never happen in the default code), 
then a `KeyError` is thrown. Example usage is shown below:
 
````
mod_dict = {'myOldMod': 'arbitrary functions and globals', 
            'myNewMod': 'arbitrary functions and globals'}
first_times = [False, True]
trusted_dict =  create_trusted_globals_dict(mod_dict, first_times)
````

Then `trusted_dict` would look like the following:

````
trusted_dict -> {'myOldMod': *trusted_values_dict['myOldModGlobals']*,
                 'myNewMod': *empty dictionary*}
````

As we can see, `trusted_dict` fetches the existing trusted values from `trusted_values_dict` and assigns them to their
respective modules if `first_time` is `True`, and assigns an empty dictionary to their respective modules if `first_time`
is `False`.

#### var_dict_to_value_dict:
`var_dict_to_value_dict` takes in a dictionary of sympy expressions `var_dict` and returns `value_dict`, which is 
equivalent to `var_dict` except each expression is evaluated by assigning every variable a pseudorandom number. Since 
the pseudorandom number is assigned according to the `seed` value in `trusted_values_dict`, we expect that over multiple
iterations, the same pseudorandom numbers will be generated in the same order. This ensures that our trusted and 
calculated `mpf` values for every variable is using the same number for each of their expressions. Say you had a dict of 
sympy expressions `var_dict` as follows,  and you want to calculate values for them according to the algorithm in 
`var_dict_to_value_dict`. Then it's as simple as doing the following:

````
var_list = {'alpha': a + b / b, 'betaU[0]': r**2, 'betaU[1]': r*cos(t) + a, 'betaU[2]': r*sin(t) + b]
value_list = list_to_value_list(var_list)
````
This will result in output of `mpf` values along the lines of:

````
value_dict -> {'alpha': mpf('2.00029108574223017786936345469015'), 
               'betaU[0]': mpf('0.00749147005858719072790563586749991'), 
               'betaU[1]': mpf('1.03311642840535379455820753746684'),
               'betaU[2]': mpf('0.734504099569315555167726889705803')}
````
In the backend, what var_dict_to_value_dict does first is get all the `free_symbols` in the values of `var_dict` and 
sort it by variable name. Then each variable is assigned a pseudorandom number using Python's `random` module. The 
standardized `seed` in `trusted_values_dict` and the sorting based on variable name act to ensure that the same variable 
gets the same number every time the code is run. Next, Sympy's `cse` (common subexpression elimination) algorithm is 
used to optimize the calculation. For example, if the term `a**2` appears multiple times in `var_list`, it is more 
efficient to store `a**2` as its own variable and replace all instances of `a**2` with that new variable. Finally, each 
variable in the optimized `var_dict` is replaced with its `mpf` value, and the resulting value dict is returned.

#### evaluate_globals:
`evaluate_globals` takes in a module dictionary `mod_dict` and a `locals()` call from where the respective 
modules were imported `old_locals`, and returns a dictionary `result_dict` that contains all modules in `mod_dict` as 
its keys and each respective module's globals as its values. Note that `locals()` is a built-in Python function, 
so it must be passed in with the parenthesis. Example usage:

````
import myModule as myMod

mod_dict = {'myMod': functions_and_globals(['myModuleInit()'], ['alpha', 'betaU'])}

result_dict = evaluate_globals(mod_dict, locals())

result_dict -> {'myMod': {'alpha': sqrt(x**2 + y**2),
                          'betaU': [1, cos(x)*sin(y), x/y]}
````

`evaluate_globals` works by using Python's `exec()` function, which allows the user to create Python code as a string, 
which we'll refer to as `stringexec`, and then pass it into `exec()` along with any globals and locals that need to be 
accessed. This allows us to call the functions and globals from `mod_dict` on our module even though they are passed
in as strings. First, `string_exec` is created by calling all the functions from `mod_dict` on their respective module
`mod`, then assigning all globals from `mod_dict` to their own variables. `stringexec` is then passed into `exec()` 
along with `old_locals`, which imports the modules that we initially imported into `exec()`. Then `exec()` runs 
`stringexec` in this environment, and we extract the resulting variables and store them in `result_dict`.

#### first_time_print:
`first_time_print` takes in a module `mod` and a value dictionary `value_dict`, and prints the code that needs to be 
copied into `trusted_values_dict` assuming the entries in `value_dict` correspond to the module `mod`. <br />
Example Usage:
````
mod = 'myModule'
value_dict = {'x': mpf('0.122483331574515176153136610247876'), 'y': mpf('0.0'), 
'z': mpf('66.6570391079152319165851690987334')}

first_time_print(mod, value_dict)
````
Output:
````
Module: myModule
Please copy the following code between the ##### and paste it into your trusted_values_dict.py file:
#####

# Generated on: 2019-06-11 12:25:19.221572
trusted_values_dict['myModuleGlobals'] = {'x': mpf('0.122483331574515176153136610247876'), 
'y': mpf('0.0'), 'z': mpf('66.6570391079152319165851690987334')}`

#####
````
Note that `first_time_print` does not check if it _should_ be called at any given time based on the existence of 
`mod` in `trusted_values_dict`. It is up to the user to determine when the correct time to call the function is. 
If `run_test` is used without any modifications (as recommended), `first_time_print` will run as determined by the 
boolean result from `is_first_time`.

#### get_variable_dimension:
`get_variable_dimension` takes in a tensor `tensor` and returns a tuple containing the rank of the tensor `dim` and the 
length of the tensor `length`. `dim` is defined as the number of dimensions of `tensor`. For example, `dim` of a scalar
is `0`, `dim` of a vector is `1`, etc. `length` is defined as the number of variables in each `dim` of `tensor`. For 
example, `length` of `[1,2,3]` is `3`. It is assumed that any tensor being passed into the function is 'square'. This
means that a `d`-dimensional tensor being passed in is made up of `n` tensors of dimension `d-1`, then each `d-1` 
dimensional tensors must also be made up of `n` tensors of dimension `d-2`, etc. `get_variable_dimension` of an empty 
list `[]` throws an `IndexError`. Example usage is shown below:

````
scalar = 2
get_variable_dimension(scalar) -> 0, 1

vector = [1, 2, 3]
get_variable_dimension(vector) -> 1, 3

long_vector = [2, 1, 3, 3, 4, 10, 12]
get_variable_dimension(long_vector) -> 1, 7

basic_tensor = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
get_variable_dimension(basic_tensor) -> 2, 3

high_dim_tensor = [[[[1]]]]
get_variable_dimension(high_dim_tensor) -> 4

empty_list = []
get_variable_dimension(empty_list) -> Raises IndexError
````

#### expand_variable_dict:
`expand_variable_dict` takes in a variable dictionary `variable_dict` and returns a dictionary containing every variable
in `variable_dict` The dict is created such that each tensor in `variable_dict` is broken down into each of its 
scalars and then put into the resulting dictionary, named using the traditional Python list syntax. 
Example usage is shown below:

````
variable_dict = {'alpha' : r / (M + r), 'betaDD': [[r * cos(theta), r * sin(theta)], 
[r * tan(theta), 0]]}

expand_variable_dict(variable_dict) -> {'alpha': r / (M + r), 'betaDD[0][0]': r * cos(theta), 
'betaDD[0][1]': r * sin(theta), 'betaDD[1][0]': r * tan(theta), 'betaDD[1][1]': 0}
````

#### is_first_time:
`is_first_time` takes in a module dictionary `mod_dict` and returns a dictionary containing corresponding booleans for 
each module in `mod_dict`. The boolean for each module `mod` is `True` if `trusted_values_dict` contains a dictionary 
entry for `mod` according to the naming convention defined in `create_trusted_globals_dict`, `False` otherwise. 
Say we have `trusted_values_dict` with the following keys (and corresponding values which aren't listed):

````
trusted_values_dict = {'Module1Globals', 'module1Globals', 'mod1globals', 'Mod2Globals', 'Module3Globs'}
````

Then we'd get the following results with the following module dictionaries (without their corresponding values):

````
mod_dict_1 = {'Module1', 'Module2', 'Module3'}
is_first_time(mod_dict_1) -> [True, False, False]

mod_dict_2 = {'Mod1', 'Mod2', 'Mod3'}
is_first_time(mod_dict_2) -> [False, True, False]

mod_dict_3 = {'module1', 'Module3', 'Mod2', 'mod1'}
is_first_time(mod_dict_3) -> [True, False, True, False]
````

An important note is that the order of keys in the module dictionary is the same order that `is_first_time` generates
its boolean list. This ensures that the resulting boolean list properly corresponds with the input module dictionary.

#### run_test:
TODO

### Example Usage:

Say you have a module `myModule` which has globals `x, y, z` and an initialization function `myModuleInit()`. How would
you go about testing these globals using the NRPy unit testing globals infrastructure? <br />
The first step is to create a unittest file. An example name is `NRPyUnitTests_Globals_Tests.py` <br />
Then it's important to import the necessary modules. Here we import unittest, which is Python's built-in testing platform, 
logging, which allows us to easily specify the level of desired output, and run_test and functions_and_globals,
as described above. <br />
````
import unittest
import logging
from run_test import run_test
from functions_and_globals import functions_and_globals
````

We then set the logging level according to the desired level of output. A good default level is INFO.
````
# ERROR -> Outputs minimal information -- only when there's an error
# INFO -> Outputs when starting and finishing a module, as well as everything in ERROR
# DEBUG -> Displays all pairs of values being compared, as well as everything in INFO
# NOTSET -> Displays symbolic dictionary for all modules, as well as everything in DEBUG

logging.basicConfig(level=logging.INFO)
````
Now that all the initialization is done, it's time to begin testing. The first step in testing is to create the unittest
class and a function for your specific test as follows. Note that the function MUST begin with the word 'test'.
````
class TestMyGlobals(unittest.TestCase):

    def test_my_module(self):
````

The next steps are to import the modules to be tested and any required pre-initialization for the modules, create a list
 of globals, and create a list of functions. We'll use `myModule` as described above.
 
````
import myModule as myMod
function_list = ['myModuleInit()']
global_list = ['x', 'y', 'z']
````

Note that the globals and functions are listed as strings! This is intentional and will cause an  error otherwise.

Next is to create our module dictionary `mod_dict`, which stores all information for our modules. This is where we use 
`functions_and_globals`, the function described above. We pass into it first our list of functions then our list of 
globals.

````
mod_dict = {'myMod': functions_and_globals(function_list, global_list)}
````

What this specifically does is create a dictionary with key `myMod` and value `functions_and_globals(function_list, 
global)list)`, which we know to be a dictionary as well from our functions tutorials.<br />
IMPORTANT NOTE: The key of our entry in `mod_dict` is `'myMod'`. It MUST be a string, and it MUST have the same name as the imported
module. Since we said `import myModule as myMod` above, its name must be `myMod`. If we had instead said
`import myModule as exampleMod` then in `ModDict` we would have to have the key as `'exampleMod'`. This is vital
for the way `runTest` functions.

The next step is to call `run_test` in order to calculate the globals for our modules. 

````
run_test(self, mod_dict, locals())
````

Finally, we must put the following if statement at the bottom of our file in order for everything to communicate 
properly.

````
if __name__ == '__main__':
unittest.main()
````

Our resulting file should look as such:

````
import unittest
import logging
from run_test import run_test
from functions_and_globals import functions_and_globals

logging.basicConfig(level=logging.INFO)

class TestMyGlobals(unittest.TestCase):

    def test_my_module(self):

        import myModule as myMod
        function_list = ['myModuleInit()']
        global_list = ['x', 'y', 'z']

        mod_dict = {'myMod': functions_and_globals(function_list, global_list)}

        run_test(self, mod_dict, locals())

if __name__ == '__main__':
    unittest.main()
````

Now that we have completed our test file, it's time to run it to calculate our trusted values for 
the first time. Run the test file, and you should see that text has been written to the console. Follow
the instructions given by copying the code in between the `#####` into your trustedValuesDict.py file.
Example output is as follows:

````
Module: myMod
Please copy the following code between the ##### and paste it into your trusted_values_dict.py file:
#####

# Generated on: 2019-06-11 12:25:19.221572
trusted_values_dict['myModGlobals'] = {'x': mpf('0.122483331574515176153136610247876'), 
'y': mpf('0.0'), 'z': mpf('66.6570391079152319165851690987334')}`

#####
````

After following the instructions, run the code again and you should see that the test passes with no errors.

Now let's say you have another module that's similar to `myModule` -- say `myModule2`, with initialization function 
`myMod2Init()` and globals `x2 y2`-- and you want to test it as well
without having to add an entirely new test file, copying all the code, etc. With some simile variable renaming, we can
add `myModule2` to our `test_my_module` function to easily test it as well. One possible way of doing this is to 
adjust our test file as follows:

````
import unittest
import logging
from run_test import run_test
from functions_and_globals import functions_and_globals

logging.basicConfig(level=logging.INFO)

class TestMyGlobals(unittest.TestCase):

    def test_my_modules(self):

        import myModule as myMod
        function_list_1 = ['myModuleInit()']
        global_list_1 = ['x', 'y', 'z']
        
        import myModule2 as myMod2
        function_list_2 = ['myMod2Init()']
        global_list_2 = ['x2', 'y2']


        mod_dict = {'myMod':  functions_and_globals(function_list_1, global_list_1),
                    'myMod2': functions_and_globals(function_list_2, global_list_2)}

        run_test(self, mod_dict, locals())

if __name__ == '__main__':
    unittest.main()
````

It's as simple as that! All we did was change the function name from `test_my_module` to `test_my_modules` (note that
this isn't necessary, but just makes things more understandable), we renamed `function_list` and `global_list` by
adding `_1` onto the end, and we imported and created our new modules, function list, and global list. Then 
by simply adding `myMod2` to the dictionary with the proper function and global list, we are now easily testing 
`myMod2` as well as the original `myMod`. Of course, after running the code for the first time you'll have to copy
the calculated values for `trusted_values_dict` into place, but you _don't_ have to redo anything to do with `myMod`.
Its tests are unaffected by the introduction of `myModule2`, and as such nothing about them will change.

