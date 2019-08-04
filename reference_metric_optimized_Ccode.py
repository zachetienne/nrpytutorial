import reference_metric as rfm
import NRPy_param_funcs as par
from outputC import *
import sympy as sp
import re

par.set_parval_from_str("reference_metric::CoordSystem","Spherical")
rfm.reference_metric()

print(rfm.f0_of_xx0)

rfm.ref_metric__hatted_quantities(scalefactor_input="generic_functions")

list_of_hatted_exprs = []
DIM=3
# First rename all functions and derivs
def make_replacements(input):
    for i in ["0","1","2"]:
        inputnew = input.replace(", xx"+i+", xx"+i+", xx"+i+")", "__dDDD"+i+i+i).\
                         replace(", xx"+i+", xx"+i+")", "__dDD"+i+i).\
                         replace(", xx"+i+")", "__dD"+i).\
                         replace("Derivative(", "").\
                         replace("f"+i+"_of_xx"+i+"_funcform(xx"+i+")", "f"+i+"_of_xx"+i+"_funcform")
        input = inputnew
    return input

rfm.detgammahat = sp.sympify(make_replacements(str(rfm.detgammahat)))
for i in range(DIM):
    rfm.ReU[i] = sp.sympify(make_replacements(str(rfm.ReU[i])))
    rfm.detgammahatdD[i] = sp.sympify(make_replacements(str(rfm.detgammahatdD[i])))
    for j in range(DIM):
        rfm.ReDD[i][j]   = sp.sympify(make_replacements(str(  rfm.ReDD[i][j])))
        rfm.ReUdD[i][j]  = sp.sympify(make_replacements(str( rfm.ReUdD[i][j])))
        rfm.ghatDD[i][j] = sp.sympify(make_replacements(str(rfm.ghatDD[i][j])))
        rfm.ghatUU[i][j] = sp.sympify(make_replacements(str(rfm.ghatUU[i][j])))
        rfm.detgammahatdDD[i][j] = sp.sympify(make_replacements(str(rfm.detgammahatdDD[i][j])))
        for k in range(DIM):
            rfm.ReDDdD[i][j][k] = sp.sympify(make_replacements(str(rfm.ReDDdD[i][j][k])))
            rfm.ReUdDD[i][j][k] = sp.sympify(make_replacements(str(rfm.ReUdDD[i][j][k])))
            rfm.ghatDDdD[i][j][k] = sp.sympify(make_replacements(str(rfm.ghatDDdD[i][j][k])))
            rfm.GammahatUDD[i][j][k] = sp.sympify(make_replacements(str(rfm.GammahatUDD[i][j][k])))
            for l in range(DIM):
                rfm.ReDDdDD[i][j][k][l] = sp.sympify(make_replacements(str(rfm.ReDDdDD[i][j][k][l])))
                rfm.ghatDDdDD[i][j][k][l] = sp.sympify(make_replacements(str(rfm.ghatDDdDD[i][j][k][l])))
                rfm.GammahatUDDdD[i][j][k][l] = sp.sympify(make_replacements(str(rfm.GammahatUDDdD[i][j][k][l])))

freevars = []
freevars.extend(rfm.detgammahat.free_symbols)
for i in range(DIM):
    freevars.extend(rfm.ReU[i].free_symbols)
    freevars.extend(rfm.detgammahatdD[i].free_symbols)
    for j in range(DIM):
        freevars.extend(rfm.ReDD[i][j].free_symbols)
        freevars.extend(rfm.ReUdD[i][j].free_symbols)
        freevars.extend(rfm.ghatDD[i][j].free_symbols)
        freevars.extend(rfm.ghatUU[i][j].free_symbols)
        freevars.extend(rfm.detgammahatdDD[i][j].free_symbols)
        for k in range(DIM):
            freevars.extend(rfm.ReDDdD[i][j][k].free_symbols)
            freevars.extend(rfm.ReUdDD[i][j][k].free_symbols)
            freevars.extend(rfm.ghatDDdD[i][j][k].free_symbols)
            freevars.extend(rfm.GammahatUDD[i][j][k].free_symbols)
            for l in range(DIM):
                freevars.extend(rfm.ReDDdDD[i][j][k][l].free_symbols)
                freevars.extend(rfm.ghatDDdDD[i][j][k][l].free_symbols)
                freevars.extend(rfm.GammahatUDDdD[i][j][k][l].free_symbols)

freevars_uniq = superfast_uniq(freevars)
print(freevars_uniq)

for var in freevars_uniq:
    basename = str(var).split("__")[0].replace("_funcform","")
    derivatv = ""
    if "__" in str(var):
        derivatv = str(var).split("__")[1].replace("_funcform","")
    if basename == "f0_of_xx0":
        basefunc = rfm.f0_of_xx0
    elif basename == "f1_of_xx1":
        basefunc = rfm.f1_of_xx1
    else:
        print("Error: function inside "+var+" undefined.")
        sys.exit(1)
    if derivatv == "":
        if basefunc == sp.sympify(0):
            var=0
    else:
        derivorder = derivatv.replace("d","").replace("D","").replace("0","0 ").replace("1","1 ").replace("2","2 ").split(" ")
        order = len(derivorder)-1
        diff_result = basefunc
        for derivdirn in derivorder:
            if derivdirn != "":
                derivwrt = rfm.xx[int(derivdirn)]
                diff_result = sp.diff(diff_result,derivwrt)
                if diff_result == sp.sympify(0):
                    var=0
        print(var,diff_result)

basemodule = sp
print("ggg",rfm.detgammahat,basemodule.diff(rfm.detgammahat,rfm.xx[0]))

freevars.extend(rfm.detgammahat.free_symbols)
for i in range(DIM):
    freevars.extend(rfm.ReU[i].free_symbols)
    freevars.extend(rfm.detgammahatdD[i].free_symbols)
    for j in range(DIM):
        freevars.extend(rfm.ReDD[i][j].free_symbols)
        freevars.extend(rfm.ReUdD[i][j].free_symbols)
        freevars.extend(rfm.ghatDD[i][j].free_symbols)
        freevars.extend(rfm.ghatUU[i][j].free_symbols)
        freevars.extend(rfm.detgammahatdDD[i][j].free_symbols)
        for k in range(DIM):
            freevars.extend(rfm.ReDDdD[i][j][k].free_symbols)
            freevars.extend(rfm.ReUdDD[i][j][k].free_symbols)
            freevars.extend(rfm.ghatDDdD[i][j][k].free_symbols)
            freevars.extend(rfm.GammahatUDD[i][j][k].free_symbols)
            for l in range(DIM):
                freevars.extend(rfm.ReDDdDD[i][j][k][l].free_symbols)
                freevars.extend(rfm.ghatDDdDD[i][j][k][l].free_symbols)
                freevars.extend(rfm.GammahatUDDdD[i][j][k][l].free_symbols)

    # if basename=="f0_of_xx0":
    #     pass
    # print(basename,derivatv)

# print(rfm.f0_of_xx0)

#.replace("(x0)","").replace(", x0","prm").replace("Derivative",""))
# newexpr = sp.sympify(str(rfm.ghatUU[0][0]).replace("f_of_xx0(xx0)","f_of_xx0").replace("Derivative","").replace(", xx0","__dD0"))
# print(newexpr)

# for i in range(DIM):
#     for j in range(DIM):
#         list_of_hatted_exprs.extend(rfm.ghatUU[i][j].free_symbols)
# print(list_of_hatted_exprs)
