#!/bin/bash

export PYTHONPATH=$PYTHONPATH:`pwd`

pypy BSSN/test/BSSN_ID_test.py
pypy BSSN/test/BSSN_RHS_test.py
