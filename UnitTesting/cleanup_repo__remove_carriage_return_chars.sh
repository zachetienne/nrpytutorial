#!/usr/bin/env python

# Replace \r\n with \n:
for file in os.listdir("."):
   with open(file,"r") as readfile: 
      with open("Tut.ipynb","w") as writefile: 
         for line in readfile.readlines(): 
            writefile.write(line.replace(r'\r\n',r'\n')) 
   os.rename("Tut.ipynb", file)
