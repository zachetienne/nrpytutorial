#!/bin/bash

# Set proper path
export PYTHONPATH=`pwd`:`pwd`/UnitTesting:`pwd`/in_progress

# Make sure the Python interpreter to use was passed in
if [ -z "$1" ]
then
    echo "ERROR: Was expecting parameter."
    echo " Usage: ./run_NRPy_UnitTests.sh [Python interpreter; e.g., python]"
    exit
fi

PYTHONEXEC=$1

echo "########################################"
echo $PYTHONPATH
echo Using $PYTHONEXEC as Python interpreter.
echo $PYTHONEXEC version info:
$PYTHONEXEC --version
echo "########################################"

# Overwrite failed_tests.txt
failed_tests_file=UnitTesting/failed_tests.txt
:> $failed_tests_file
printf "Failures:\n\n" > $failed_tests_file

# Usage: add_test [path to test file]
add_test () {
  $PYTHONEXEC $1 $PYTHONEXEC $failed_tests_file $2 $3
}

# Change boolean to true/false depending on if you want the tests that fail to automatically re-run with
# logging_level='DEBUG'
# By default, this boolean is set to true. A motivation for changing it to false would be if every test is failing, and
# the sheer amount of output given with the logging level is too much to process. Having it set to true is most useful
# when only one or two modules fail.
rerun_if_fail=false

# $PYTHONEXEC UnitTesting/Test_UnitTesting/test_functions.py

add_test CarpetX/tests/test_WaveToy.py
add_test UnitTesting/Test_UnitTesting/test_module.py
add_test BSSN/tests/test_BSSN.py
add_test GRHD/tests/test_GRHD.py
add_test GRFFE/tests/test_GRFFE.py
add_test GRMHD/tests/test_GRMHD.py
add_test FishboneMoncriefID/tests/test_FishboneMoncriefID.py
#add_test in_progress-Maxwell/Maxwell/tests/test_Maxwell.py
add_test ScalarWave/tests/test_ScalarWave.py
add_test tests/test_reference_metric.py
add_test u0_smallb_Poynting__Cartesian/tests/test_u0_smallb_Poynting__Cartesian.py
add_test WeylScal4NRPy/tests/test_WeylScal4NRPy.py

# TODO: add your tests here
echo "Starting doctest unit tests!"
failed_unittest=0
for file in expr_tree.py indexedexp.py loop.py functional.py finite_difference_helpers.py assert_equal.py sugar.py; do
    echo Running doctest on file: $file
    $PYTHONEXEC -m doctest $file
    if [ $? == 1 ]
    then
        failed_unittest=1
    fi
    echo Doctest of $file finished.
done
$PYTHONEXEC -c "import sys, sympy; major, minor = int(sympy.__version__.split('.')[0]), int(sympy.__version__.split('.')[1]); sys.exit(major == 1 and minor <= 3)"
if [ $? == 0 ]
then
    echo Running doctest on file: cse_helpers.py
    $PYTHONEXEC -m doctest cse_helpers.py
    if [ $? == 1 ]
    then
        failed_unittest=1
    fi
    echo Doctest of cse_helpers.py finished.
fi
for file in tests/test_parse_BSSN.py; do
    echo Running unittest on file: $file
    $PYTHONEXEC $file
    if [ $? == 1 ]
    then
        failed_unittest=1
    fi
    echo Unittest of $file finished.
done

# Checking failed_tests.txt to see what failed
contents=$(<$failed_tests_file)

if [ "$contents" == $"Failures:" ] && [ $failed_unittest == 0 ]
then
  printf "All tests passed!\n\n"
  exit 0
else
  printf "Tests failed!\n\n"
  printf "$contents \n\n"
  printf '%s\n' '----------------------------------------------------------------------'

  # If tests failed and rerun_if_fail is true, the rerun failed tests with logging level DEBUG
  if $rerun_if_fail
  then
    printf "Re-running failed tests with logging_level=DEBUG:\n\n"
    while IFS=': ' read -r col1 col2
    do
      if [ "$col1" != "Failures" ] && [ "$col1" != "" ]
      then
        add_test "$col1" "DEBUG" "$col2"
      fi
    done <UnitTesting/failed_tests.txt

    printf "Completed by re-running following tests with logging_level=DEBUG\n\n"
    printf "$contents \n\n"
    printf '%s\n' '----------------------------------------------------------------------'
  fi

  exit 1
fi
