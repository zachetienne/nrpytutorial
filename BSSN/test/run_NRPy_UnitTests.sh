#!/bin/bash

export PYTHONPATH=$PYTHONPATH:`pwd`

PYTHONEXEC=python
if [ -x "$(command -v pypy)" ]; then
    # pypy is NOT installed
    PYTHONEXEC=pypy
fi

echo "########################################"
echo Using $PYTHONEXEC as Python interpreter.
echo $PYTHONEXEC version info:
$PYTHONEXEC --version
echo "########################################"

$PYTHONEXEC BSSN/test/NRPyUnitTests_Function_Tests.py &&
$PYTHONEXEC BSSN/test/NRPyUnitTests_Globals_Tests.py
