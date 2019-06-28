#!/bin/bash

export PYTHONPATH=`pwd`:`pwd`/UnitTesting


PYTHONEXEC=$1

echo "########################################"
echo $PYTHONPATH
echo Using $PYTHONEXEC as Python interpreter.
echo $PYTHONEXEC version info:
$PYTHONEXEC --version
echo "########################################"

$PYTHONEXEC UnitTesting/NRPyUnitTests_Functions.py &&
$PYTHONEXEC BSSN/tests/NRPyUnitTests_BSSN_Globals.py &&
$PYTHONEXEC tests/NRPyUnitTests_Reference_Metric_Globals.py &&
$PYTHONEXEC FishboneMoncriefID/tests/NRPyUnitTests_FishboneMoncriefID_Globals.py &&
$PYTHONEXEC GiRaFFE_HO/tests/NRPyUnitTests_GiRaFFE_HO_Globals.py &&
$PYTHONEXEC GiRaFFEfood_HO/tests/NRPyUnitTests_GiRaFFEfood_HO_Globals.py &&
$PYTHONEXEC ScalarWave/tests/NRPyUnitTests_ScalarWave_Globals.py &&
$PYTHONEXEC ScalarWaveCurvilinear/tests/NRPyUnitTests_ScalarWaveCurvilinear_Globals.py &&
$PYTHONEXEC Maxwell/tests/NRPyUnitTests_Maxwell_Globals.py &&
$PYTHONEXEC u0_smallb_Poynting__Cartesian/tests/NRPyUnitTests_u0_smallb_Poynting__Cartesian_Globals.py
