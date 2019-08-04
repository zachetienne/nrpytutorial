import reference_metric as rfm
import NRPy_param_funcs as par
from outputC import *
import sympy as sp
import indexedexp as ixp

def reference_metric_optimized_Ccode(SymPySimplifyExpressions=True):

    # Step 0: Set dimension to zero
    DIM=3

    # Step 2: Declare all globals
    global ReU,ReDD,ghatDD,ghatUU,detgammahat
    ReU    = ixp.zerorank1()
    ReDD   = ixp.zerorank2()
    ghatDD = ixp.zerorank2()
    ghatUU = ixp.zerorank2()

    global detgammahatdD, detgammahatdDD
    detgammahatdD = ixp.zerorank1(DIM)
    detgammahatdDD = ixp.zerorank2(DIM)

    global ReUdD,ReUdDD
    ReUdD  = ixp.zerorank2(DIM)
    ReUdDD = ixp.zerorank3(DIM)

    global ReDDdD,ReDDdDD
    ReDDdD = ixp.zerorank3(DIM)
    ReDDdDD = ixp.zerorank4(DIM)

    global ghatDDdD, ghatDDdDD
    ghatDDdD = ixp.zerorank3(DIM)
    ghatDDdDD = ixp.zerorank4(DIM)

    global GammahatUDD
    GammahatUDD = ixp.zerorank3(DIM)

    global GammahatUDDdD
    GammahatUDDdD = ixp.zerorank4(DIM)

    # Step 3: Define all hatted quantities in terms of generic SymPy functions for scale factors,
    #         which in turn were defined in rfm.reference_metric()
    # Step 1: Set up reference metric.
    rfm.reference_metric(SymPySimplifyExpressions)
    rfm.ref_metric__hatted_quantities(SymPySimplifyExpressions,scalefactor_input="generic_functions")
    print(rfm.detgammahat)

    # Step 4: Now that all hatted quantities are written in terms of generic SymPy functions,
    #         we will now replace SymPy functions with simple variables using rigid NRPy+ syntax,
    #         and store these variables to globals defined above.
    def make_replacements(input):
        for i in ["0","1","2"]:
            inputnew = input.replace(", xx"+i+", xx"+i+", xx"+i+")", "__dDDD"+i+i+i).\
                             replace(", xx"+i+", xx"+i+")", "__dDD"+i+i).\
                             replace(", xx"+i+")", "__dD"+i).\
                             replace("Derivative(", "").\
                             replace("f"+i+"_of_xx"+i+"_funcform(xx"+i+")", "f"+i+"_of_xx"+i+"_funcform")
            input = inputnew
        return input
    #
    detgammahat = sp.sympify(make_replacements(str(rfm.detgammahat)))
    for i in range(DIM):
        ReU[i] = sp.sympify(make_replacements(str(rfm.ReU[i])))
        detgammahatdD[i] = sp.sympify(make_replacements(str(rfm.detgammahatdD[i])))
        for j in range(DIM):
            ReDD  [i][j] = sp.sympify(make_replacements(str(rfm.  ReDD[i][j])))
            ReUdD [i][j] = sp.sympify(make_replacements(str(rfm. ReUdD[i][j])))
            ghatDD[i][j] = sp.sympify(make_replacements(str(rfm.ghatDD[i][j])))
            ghatUU[i][j] = sp.sympify(make_replacements(str(rfm.ghatUU[i][j])))
            detgammahatdDD[i][j] = sp.sympify(make_replacements(str(rfm.detgammahatdDD[i][j])))
            for k in range(DIM):
                ReDDdD     [i][j][k] = sp.sympify(make_replacements(str(rfm.ReDDdD[i][j][k])))
                ReUdDD     [i][j][k] = sp.sympify(make_replacements(str(rfm.ReUdDD[i][j][k])))
                ghatDDdD   [i][j][k] = sp.sympify(make_replacements(str(rfm.ghatDDdD[i][j][k])))
                GammahatUDD[i][j][k] = sp.sympify(make_replacements(str(rfm.GammahatUDD[i][j][k])))
                for l in range(DIM):
                    ReDDdDD[i][j][k][l] = sp.sympify(make_replacements(str(rfm.ReDDdDD[i][j][k][l])))
                    ghatDDdDD[i][j][k][l] = sp.sympify(make_replacements(str(rfm.ghatDDdDD[i][j][k][l])))
                    GammahatUDDdD[i][j][k][l] = sp.sympify(make_replacements(str(rfm.GammahatUDDdD[i][j][k][l])))

    # At this point, each expression is written in terms of the generic functions of
    #    xx0, xx1, and/or xx2 and their derivatives. Depending on the functions, some
    #    of these derivatives may be zero. In Step 5 we'll evaluate the function
    #    derivatives exactly and set the expressions to zero. Otherwise in the C code
    #    we'd be storing performing arithmetic with zeros -- wasteful!

    # Step 5.a: Construct the full list of *unique* NRPy+ variables representing the
    #           SymPy functions and derivatives, so that all zero derivatives can be
    #           computed.
    freevars = []
    freevars.extend(detgammahat.free_symbols)
    for i in range(DIM):
        freevars.extend(ReU[i].free_symbols)
        freevars.extend(detgammahatdD[i].free_symbols)
        for j in range(DIM):
            freevars.extend(ReDD[i][j].free_symbols)
            freevars.extend(ReUdD[i][j].free_symbols)
            freevars.extend(ghatDD[i][j].free_symbols)
            freevars.extend(ghatUU[i][j].free_symbols)
            freevars.extend(detgammahatdDD[i][j].free_symbols)
            for k in range(DIM):
                freevars.extend(ReDDdD[i][j][k].free_symbols)
                freevars.extend(ReUdDD[i][j][k].free_symbols)
                freevars.extend(ghatDDdD[i][j][k].free_symbols)
                freevars.extend(GammahatUDD[i][j][k].free_symbols)
                for l in range(DIM):
                    freevars.extend(ReDDdDD[i][j][k][l].free_symbols)
                    freevars.extend(ghatDDdDD[i][j][k][l].free_symbols)
                    freevars.extend(GammahatUDDdD[i][j][k][l].free_symbols)

    freevars_uniq = superfast_uniq(freevars)

    print(freevars_uniq)
    # Step 5.b: Using the expressions rfm.f?_of_xx? set in rfm.reference_metric(),
    #           evaluate each needed derivative and, in the case it is zero,
    #           set the corresponding "freevar" variable to zero.
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
            print("Error: function inside "+str(var)+" undefined.")
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
                        print("Attempting to set "+str(var)+" to zero")
                        var=0
            print(var,diff_result)

    # Step 5.c:
    print(freevars_uniq)

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

par.set_parval_from_str("reference_metric::CoordSystem","Spherical")
reference_metric_optimized_Ccode()
