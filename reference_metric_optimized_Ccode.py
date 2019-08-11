# reference_metric_optimized_Ccode.py: Generate C code for
#     precomputing and storing in memory complicated expressions
#     related to the reference metric (a.k.a., "hatted quantities")

# The precomputed "hatted quantity" expressions will be stored in
#    a C struct called rfmstruct. As these expressions generally
#    involve computationally expensive transcendental functions
#    of xx0,xx1,or xx2, and xx0,xx1, and xx2 remain fixed across
#    most (if not all) of a given simulation, setting up the
#    rfmstruct can greatly improve performance.

# The core challenge in setting up the rfmstruct is collecting
#    all the information needed to automatically generate it.
# To this end, this module implements the following algorithm:
# Step 0: Set dimension, DIM, to three
# Step 1: Declare all globals
# Step 2: Set up reference metric.
# Step 3: Define all hatted quantities in terms of generic SymPy functions for scale factors,
#         which in turn were defined in rfm.reference_metric()
# Step 4: Now that all hatted quantities are written in terms of generic SymPy functions,
#         we will now replace SymPy functions with simple variables using rigid NRPy+ syntax,
#         and store these variables to globals defined above.
# Step 5: At this point, each expression is written in terms of the generic functions
#         of xx0, xx1, and/or xx2 and their derivatives. Depending on the functions, some
#         of these derivatives may be zero. In Step 5 we'll evaluate the function
#         derivatives exactly and set the expressions to zero. Otherwise in the C code
#         we'd be storing performing arithmetic with zeros -- wasteful!
# Step 5.a: Construct the full list of *unique* NRPy+ variables representing the
#           SymPy functions and derivatives, so that all zero derivatives can be
#           computed.
# Step 5.b: Using the expressions rfm.f?_of_xx? set in rfm.reference_metric(),
#           evaluate each needed derivative and, in the case it is zero,
#           set the corresponding "freevar" variable to zero.
# Step 5.c: Substitute zero for all functions & derivatives that evaluate to zero.
# Step 6: Construct needed C code for declaring rfmstruct, allocating storage for
#         rfmstruct arrays, defining each element in each array, reading the
#         rfmstruct data from memory (both with and without SIMD enabled), and
#         freeing allocated memory for the rfmstrcut arrays.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

import reference_metric as rfm
import NRPy_param_funcs as par
from outputC import *
import sympy as sp
import indexedexp as ixp
import grid as gri

def reference_metric_optimized_Ccode(directory, SymPySimplifyExpressions=True):

    # Step 0: Set dimension to three
    DIM=3

    # Step 1: Declare all globals
    global ReU,ReDD,ghatDD,ghatUU,detgammahat
    ReU    = ixp.zerorank1()
    ReDD   = ixp.zerorank2()
    ghatDD = ixp.zerorank2()
    ghatUU = ixp.zerorank2()

    global detgammahatdD, detgammahatdDD
    detgammahatdD  = ixp.zerorank1(DIM)
    detgammahatdDD = ixp.zerorank2(DIM)

    global ReUdD,ReUdDD
    ReUdD  = ixp.zerorank2(DIM)
    ReUdDD = ixp.zerorank3(DIM)

    global ReDDdD,ReDDdDD
    ReDDdD  = ixp.zerorank3(DIM)
    ReDDdDD = ixp.zerorank4(DIM)

    global ghatDDdD, ghatDDdDD
    ghatDDdD = ixp.zerorank3(DIM)
    ghatDDdDD = ixp.zerorank4(DIM)

    global GammahatUDD
    GammahatUDD = ixp.zerorank3(DIM)

    global GammahatUDDdD
    GammahatUDDdD = ixp.zerorank4(DIM)

    # Step 2: Set up reference metric.
    rfm.reference_metric(SymPySimplifyExpressions)

    # Step 3: Define all hatted quantities in terms of generic SymPy functions for scale factors,
    #         which in turn were defined in rfm.reference_metric()
    rfm.ref_metric__hatted_quantities(SymPySimplifyExpressions,scalefactor_input="generic_functions")

    # Step 4: Now that all hatted quantities are written in terms of generic SymPy functions,
    #         we will now replace SymPy functions with simple variables using rigid NRPy+ syntax,
    #         and store these variables to globals defined above.
    def make_replacements(input):
        for i in ["0", "1", "2"]:
            inputnew = input.replace("Derivative(f"+i+"_of_xx"+i+"_funcform(xx"+i+"), (xx"+i+", 3))",
                                     "Derivative(f"+i+"_of_xx"+i+"_funcform(xx"+i+"), xx"+i+", xx"+i+", xx"+i+")"). \
                             replace("Derivative(f"+i+"_of_xx"+i+"_funcform(xx"+i+"), (xx"+i+", 2))",
                                     "Derivative(f"+i+"_of_xx"+i+"_funcform(xx"+i+"), xx"+i+", xx"+i+")"). \
                             replace(", xx"+i+", xx"+i+", xx"+i+")", "__DDD"+i+i+i).\
                             replace(", xx"+i+", xx"+i+")", "__DD"+i+i).\
                             replace(", xx"+i+")", "__D"+i).\
                             replace("f"+i+"_of_xx"+i+"_funcform(xx"+i+")", "f"+i+"_of_xx"+i)
            input = inputnew
        inputnew = input.replace("Derivative(", "")
        if "Derivative" in inputnew:
            print("Error: ",inputnew)
            sys.exit(1)
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

    # Step 5: At this point, each expression is written in terms of the generic functions
    #         of xx0, xx1, and/or xx2 and their derivatives. Depending on the functions, some
    #         of these derivatives may be zero. In Step 5 we'll evaluate the function
    #         derivatives exactly and set the expressions to zero. Otherwise in the C code
    #         we'd be storing performing arithmetic with zeros -- wasteful!

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

    freevars_uniq_zeroed = []
    for i in range(len(freevars_uniq)):
        freevars_uniq_zeroed.append(freevars_uniq[i])

    # Step 5.b: Using the expressions rfm.f?_of_xx? set in rfm.reference_metric(),
    #           evaluate each needed derivative and, in the case it is zero,
    #           set the corresponding "freevar" variable to zero.
    freevars_uniq_vals = []
    for i in range(len(freevars_uniq)):
        var = freevars_uniq[i]
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
        diff_result = basefunc
        if derivatv == "":
            pass
        else:
            derivorder = derivatv.replace("d","").replace("D","").replace("0","0 ").replace("1","1 ").replace("2","2 ").split(" ")
            for derivdirn in derivorder:
                if derivdirn != "":
                    derivwrt = rfm.xx[int(derivdirn)]
                    diff_result = sp.diff(diff_result,derivwrt)
        freevars_uniq_vals.append(diff_result)
        if diff_result == sp.sympify(0):
            freevars_uniq_zeroed[i]  = 0

    # Step 5.c: Finally, substitute zero for all functions & derivatives that evaluate to zero.
    for varidx in range(len(freevars_uniq)):
        detgammahat = detgammahat.subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
        for i in range(DIM):
            ReU[i] = ReU[i].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
            detgammahatdD[i] = detgammahatdD[i].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
            for j in range(DIM):
                ReDD  [i][j] = ReDD  [i][j].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
                ReUdD [i][j] = ReUdD [i][j].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
                ghatDD[i][j] = ghatDD[i][j].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
                ghatUU[i][j] = ghatUU[i][j].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
                # print(varidx,i,j,freevars_uniq[varidx],freevars_uniq_zeroed[varidx],detgammahatdDD[i][j])
                detgammahatdDD[i][j] = detgammahatdDD[i][j].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
                # print(varidx,i,j,freevars_uniq[varidx],freevars_uniq_zeroed[varidx],detgammahatdDD[i][j])
                for k in range(DIM):
                    ReDDdD     [i][j][k] = ReDDdD     [i][j][k].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
                    ReUdDD     [i][j][k] = ReUdDD     [i][j][k].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
                    ghatDDdD   [i][j][k] = ghatDDdD   [i][j][k].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
                    GammahatUDD[i][j][k] = GammahatUDD[i][j][k].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
                    for l in range(DIM):
                        ReDDdDD      [i][j][k][l] = ReDDdDD      [i][j][k][l].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
                        ghatDDdDD    [i][j][k][l] = ghatDDdDD    [i][j][k][l].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])
                        GammahatUDDdD[i][j][k][l] = GammahatUDDdD[i][j][k][l].subs(freevars_uniq[varidx],freevars_uniq_zeroed[varidx])

    # Step 6: Construct needed C code for declaring rfmstruct, allocating storage for
    #         rfmstruct arrays, defining each element in each array, reading the
    #         rfmstruct data from memory (both with and without SIMD enabled), and
    #         freeing allocated memory for the rfmstrcut arrays.
    # struct_str: String that declares the rfmstruct struct.
    struct_str = "typedef struct __rfmstruct__ {\n"
    # rfmstruct stores pointers to (so far) 1D arrays. The malloc_str string allocates space for the arrays.
    malloc_str = "rfm_struct rfmstruct;\n"
    # define_str sets the arrays to appropriate values. Note that elements of
    #    these arrays will generally be transcendental functions of xx0,xx1,or xx2.
    #    Since xx0,xx1, and xx2 remain fixed across many (if not all) iterations,
    #    and these transcendental functions are quite expensive, setting up this
    #    struct can greatly improve performance.
    define_str = ""
    # readvr_str reads the arrays from memory as needed
    readvr_str = ["","",""]
    readvr_SIMD_outer_str = ["","",""]
    readvr_SIMD_inner_str = ["","",""]
    freemm_str = ""
    for dirn in [0,1,2]:
        malloc_size = gri.Nxx_plus_2NGHOSTS[dirn]
#        malloc_size = "Nxx_plus_2NGHOSTS["+str(dirn)+"]"

        numvars = 0
        for varidx in range(len(freevars_uniq)):
            if "xx"+str(dirn) in str(freevars_uniq_zeroed[varidx]):
                numvars = numvars + 1
                struct_str += "\tREAL *restrict "+str(freevars_uniq_zeroed[varidx])+";\n"
                malloc_str += "rfmstruct."+str(freevars_uniq_zeroed[varidx])+" = (REAL *)malloc(sizeof(REAL)*"+str(malloc_size)+");\n"
                freemm_str += "free(rfmstruct."+str(freevars_uniq_zeroed[varidx])+");\n"

                define_str += """
for(int ii=0;ii<"""+str(malloc_size)+""";ii++) {
    const REAL xx"""+str(dirn)+""" = xx["""+str(dirn)+"""][ii];
    rfmstruct."""+str(freevars_uniq_zeroed[varidx])+"""[ii] = """+str(sp.ccode(freevars_uniq_vals[varidx]))+""";
}"""
                readvr_str[dirn] += "const REAL "+str(freevars_uniq_zeroed[varidx])+" = rfmstruct->"+str(freevars_uniq_zeroed[varidx])+"[i"+str(dirn)+"];\n"
                readvr_SIMD_outer_str[dirn] += "const double NOSIMD"+str(freevars_uniq_zeroed[varidx])+" = rfmstruct->"+str(freevars_uniq_zeroed[varidx])+"[i"+str(dirn)+"]; "
                readvr_SIMD_outer_str[dirn] += "const REAL_SIMD_ARRAY "+str(freevars_uniq_zeroed[varidx])+" = ConstSIMD(NOSIMD"+str(freevars_uniq_zeroed[varidx])+");\n"
                readvr_SIMD_inner_str[dirn] += "const REAL_SIMD_ARRAY "+str(freevars_uniq_zeroed[varidx])+" = ReadSIMD(&rfmstruct->"+str(freevars_uniq_zeroed[varidx])+"[i"+str(dirn)+"]);\n"

    struct_str += "} rfm_struct;\n\n"

    # Step 7: Output needed C code to files
    with open(directory+"/rfm_struct__declare.h", "w") as file:
        file.write(struct_str)
    with open(directory+"/rfm_struct__malloc.h", "w") as file:
        file.write(malloc_str)
    with open(directory+"/rfm_struct__define.h", "w") as file:
        file.write(define_str)
    for i in range(3):
        with open(directory + "/rfm_struct__read" + str(i) + ".h", "w") as file:
            file.write(readvr_str[i])
        with open(directory + "/rfm_struct__SIMD_outer_read" + str(i) + ".h", "w") as file:
            file.write(readvr_SIMD_outer_str[i])
        with open(directory + "/rfm_struct__SIMD_inner_read" + str(i) + ".h", "w") as file:
            file.write(readvr_SIMD_inner_str[i])
    with open(directory+"/rfm_struct__freemem.h", "w") as file:
        file.write(freemm_str)


    # print(struct_str)
    # print(malloc_str)
    # print(define_str)
    # print(readvr_str[0])
    # print(readvr_str[1])
    # print(readvr_str[2])
    # print(freemm_str)

# print(rfm.f0_of_xx0)

#.replace("(x0)","").replace(", x0","prm").replace("Derivative",""))
# newexpr = sp.sympify(str(rfm.ghatUU[0][0]).replace("f_of_xx0(xx0)","f_of_xx0").replace("Derivative","").replace(", xx0","__dD0"))
# print(newexpr)

# for i in range(DIM):
#     for j in range(DIM):
#         list_of_hatted_exprs.extend(rfm.ghatUU[i][j].free_symbols)
# print(list_of_hatted_exprs)

# par.set_parval_from_str("reference_metric::CoordSystem","SinhSpherical")
# reference_metric_optimized_Ccode("ScalarWaveCurvilinear")
