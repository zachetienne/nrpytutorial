#!/bin/bash

# Remove trailing whitespace
sed -i 's/[ \t]*$//' $1

# Add newline at end of file
sed -i '$a\' $1
