#import NRPy_param_funcs as par
#from outputC import *

#import sympy as sp

# Step 1: Initialize free parameters for this module:
# modulename here will be set to "scalarwave", based on the filename.
thismodule = __name__

#def boundary_conditions():
    # Step 3a: Register gridfunctions that are needed as input.
#    DIM = par.parval_from_str("DIM")

#import numpy as np


def ternary(x, D):
    nums = []
    for i in range(D):
        x, r = divmod(x, 3)
        nums.append(r - 1)
    return nums

D = 5
Ngz = 3**D

########################################################################################################################

idx_list = "i0"
idx_g_list = "i0_g"
for dim in range(1, D):
    idx_list += ", i" + str(dim)
    idx_g_list += ", i" + str(dim) + "_g"

idx_string = "(i" + str(D-1) + ")"
for dim in range(D - 1, 0, -1):
    idx_string = "( (i" + str(dim - 1) + ") + Npts[" + str(dim - 1) + "] * " + idx_string + " )"
print "#define IDX" + str(D) + "(" + idx_list + ") " + idx_string

########################################################################################################################

rgz = []
for i in range(Ngz):
    sum = 0
    num = ternary(i, D)
    for j in num:
        sum += j * j
    rgz.append([num, sum])

rgz_sorted = sorted(rgz, key=lambda x: x[1])


for i in range(Ngz):
    print "// Mapping region " + str(rgz_sorted[i][0])

    for dim in range(D - 1, -1, -1):
        if rgz_sorted[i][0][dim] == 0:
            print "for(int i" + str(dim) + "_g = NGHOSTS; i" + str(dim) + "_g < Npts[" + str(dim) + "] - NGHOSTS; i" + str(dim) + "_g++) {"
        else:
            print "for(int i" + str(dim) + " = 0; i" + str(dim) + " < NGHOSTS; i" + str(dim) + "++) {"
            if rgz_sorted[i][0][dim] == -1:
                print "const int i" + str(dim) + "_g = NGHOSTS - 1 - i" + str(dim) + ";\n"
            elif rgz_sorted[i][0][dim] == 1:
                print "const int i" + str(dim) + "_g = Npts[" + str(dim) + "] - NGHOSTS + i" + str(dim) + ";\n"

    print "Map_Ghost_To_Available_Point(IDX" + str(D) + "(" + idx_g_list + "), Cartxyz, ListOfSetPoints, &NumberOfSetPoints, SymmetryMap);"

    for dim in range(D):
        print "}"

    print ""

########################################################################################################################
########################################################################################################################

loop_list = ""
for dim in range(D):
    loop_list = "for(int i" + str(dim) + " = 0; i" + str(dim) + " < Npts[" + str(dim) + "]; i" + str(dim) + "++) " + loop_list

print "#define LOOP_GZFILL(" + idx_list + ") _Pragma(\"omp parallel for\") \\\n  " + loop_list

########################################################################################################################

loop_list = ""
for dim in range(D):
    loop_list = "for(int i" + str(dim) + " = NGHOSTS; i" + str(dim) + " < Npts[" + str(dim) + "] - NGHOSTS; i" + str(dim) + "++) " + loop_list

print "#define LOOP_NOGZFILL(" + idx_list + ") _Pragma(\"omp parallel for\") \\\n  " + loop_list

########################################################################################################################

loop_list = ""
for dim in range(D):
    loop_list = "  for(int i" + str(dim) + " = NGHOSTS; i" + str(dim) + " < Npts[" + str(dim) + "] - NGHOSTS; i" + str(dim) + "++) { \\\n    const REAL x" + str(dim) + " = x" + str(dim) + "G[i" + str(dim) + "]; \\\n" + loop_list
    if dim == 0:
        loop_list = "  _Pragma(\"ivdep\") \\\n  _Pragma(\"vector always\") \\\n" + loop_list + "    const int idx = IDX" + str(D) + "(" + idx_list + "); \\"

print "#define START_LOOP_NOGZFILL(" + idx_list + ") _Pragma(\"omp parallel for\") \\\n" + loop_list

########################################################################################################################

loop_list = ""
for dim in range(D):
    loop_list = "  for(int i" + str(dim) + " = 0; i" + str(dim) + " < Npts[" + str(dim) + "]; i" + str(dim) + "++) { \\\n    const REAL x" + str(dim) + " = x" + str(dim) + "G[i" + str(dim) + "]; \\\n" + loop_list
    if dim == 0:
        loop_list = "  _Pragma(\"ivdep\") \\\n  _Pragma(\"vector always\") \\\n" + loop_list + "    const int idx = IDX" + str(D) + "(" + idx_list + "); \\"

print "#define START_LOOP_GZFILL(" + idx_list + ") _Pragma(\"omp parallel for\") \\\n" + loop_list

########################################################################################################################

loop_list = ""
for dim in range(D):
    loop_list += "} "

print "#define END_LOOP_GZFILL " + loop_list

print "#define END_LOOP_NOGZFILL " + loop_list

########################################################################################################################

loop_list = ""
for dim in range(D):
    if dim == 0:
        loop_list = "  for(int i0 = 0; i0 < Npts[0]; i0 += 4) { \\\n    const int idx = IDX" + str(D) + "(" + idx_list + "); \\\n    const __m256d x0 = _mm256_loadu_pd(&x0G[i0]); \\"
    else:
        loop_list = "  for(int i" + str(dim) + " = 0; i" + str(dim) + " < Npts[" + str(dim) + "]; i" + str(dim) + "++) { \\\n    const __m256d x" + str(dim) + " = _mm256_loadu_pd(&x" + str(dim) + "G[i" + str(dim) + "]); \\\n" + loop_list

print "#define START_LOOP_NOGZFILL_SIMD(" + idx_list + ") \\\n" + loop_list

########################################################################################################################

loop_list = ""
for dim in range(D):
    if dim == 0:
        loop_list = "  for(int i0 = 0; i0 < Npts[0]; i0 += 8) { \\\n    const int idx = IDX" + str(D) + "(" + idx_list + "); \\\n    const __m512d x0 = _mm512_loadu_pd(&x0G[i0]); \\"
    else:
        loop_list = "  for(int i" + str(dim) + " = 0; i" + str(dim) + " < Npts[" + str(dim) + "]; i" + str(dim) + "++) { \\\n    const __m512d x" + str(dim) + " = _mm512_loadu_pd(&x" + str(dim) + "G[i" + str(dim) + "]); \\\n" + loop_list

print "#define START_LOOP_NOGZFILL_SIMD(" + idx_list + ") \\\n" + loop_list
