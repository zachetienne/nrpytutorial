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

$PYTHONEXEC UnitTesting/NRPyUnitTests_Functions.py #&&
$PYTHONEXEC BSSN/tests/NRPyUnitTests_BSSN_Globals.py &&
$PYTHONEXEC tests/NRPyUnitTests_Reference_Metric_Globals.py &&
$PYTHONEXEC FishboneMoncriefID/tests/NRPyUnitTests_FishboneMoncriefID_Globals.py &&
$PYTHONEXEC GiRaFFE_HO/tests/NRPyUnitTests_GiRaFFE_HO_Globals.py &&
$PYTHONEXEC GiRaFFEfood_HO/tests/NRPyUnitTests_GiRaFFEfood_HO_Globals.py &&
$PYTHONEXEC ScalarWave/tests/NRPyUnitTests_ScalarWave_Globals.py &&
$PYTHONEXEC ScalarWaveCurvilinear/tests/NRPyUnitTests_ScalarWaveCurvilinear_Globals.py
