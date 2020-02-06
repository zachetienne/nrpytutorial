#!/usr/bin/env python

import os

# Replace \r\n with \n:
for file in os.listdir("."):
   # `file` can be a file or a directory. Check that it's in fact a file:
   if os.path.isfile(file):
      with open(file,"r") as readfile:
         with open("Tut.ipynb","w") as writefile:
            for line in readfile.readlines():
               writefile.write(line.replace(r'\r\n',r'\n'))
      os.rename("Tut.ipynb", file)
