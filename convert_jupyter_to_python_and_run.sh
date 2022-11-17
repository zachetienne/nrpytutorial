#!/bin/bash

jupyter nbconvert --to python $1 --output=blah
ipython3 blah.py
rm -f blah.py
