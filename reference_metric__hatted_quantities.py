import sympy as sp

import NRPy_param_funcs as par
# from outputC import *
import indexedexp as ixp

import reference_metric as rfm

# Step 0a: Initialize parameters
thismodule = __name__
# par.initialize_param(par.glb_param("char", thismodule, "CoordSystem", "Cartesian"))

def reference_metric__hatted_quantities(call_reference_metric=False):
    if call_reference_metric:
        rfm.reference_metric()
    
    # Step 0: Set dimension DIM
    DIM = par.parval_from_str("grid::DIM")

    global ReU,ReDD,ghatDD
    ReU    = ixp.zerorank1()
    ReDD   = ixp.zerorank2()
    ghatDD = ixp.zerorank2()

    # Step 1: Compute ghatDD (reference metric), as well as 
    #         rescaling vector ReU & rescaling matrix ReDD
    for i in range(DIM):
        rfm.scalefactor_orthog[i] = sp.sympify(rfm.scalefactor_orthog[i])
        ghatDD[i][i] = rfm.scalefactor_orthog[i]**2
        ReU[i] = 1/rfm.scalefactor_orthog[i]
        for j in range(DIM):
            ReDD[i][j] = rfm.scalefactor_orthog[i]*rfm.scalefactor_orthog[j]
    # Step 1b: Compute ghatUU
    ghatUU, detgammahat = ixp.symm_matrix_inverter3x3(ghatDD)

    # Step 1c: Sanity check: verify that ReDD, ghatDD, 
    #          and ghatUU are all symmetric rank-2: 
    for i in range(DIM):
        for j in range(DIM):
            if ReDD[i][j] != ReDD[j][i]:
                print("Error: ReDD["+ str(i) + "][" + str(j) + "] != ReDD["+ str(j) + "][" + str(i) + ": " + str(ReDD[i][j]) + "!=" + str(ReDD[j][i]))
                exit(1)
            if ghatDD[i][j] != ghatDD[j][i]:
                print("Error: ghatDD["+ str(i) + "][" + str(j) + "] != ghatDD["+ str(j) + "][" + str(i) + ": " + str(ghatDD[i][j]) + "!=" + str(ghatDD[j][i]))
                exit(1)
            if ghatUU[i][j] != ghatUU[j][i]:
                print("Error: ghatUU["+ str(i) + "][" + str(j) + "] != ghatUU["+ str(j) + "][" + str(i) + ": " + str(ghatUU[i][j]) + "!=" + str(ghatUU[j][i]))
                exit(1)

    # Step 2: Compute det(ghat) and its 1st & 2nd derivatives
    global detgammahatdD,detgammahatdDD
    detgammahatdD  = ixp.zerorank1(DIM)
    detgammahatdDD = ixp.zerorank2(DIM)
    for i in range(DIM):
        detgammahatdD[i] = (sp.diff(detgammahat, rfm.xx[i]))
        for j in range(DIM):
            detgammahatdDD[i][j] = sp.diff(detgammahatdD[i], rfm.xx[j])

    # Step 3a: Compute 1st & 2nd derivatives of rescaling vector.
    #          (E.g., needed in BSSN for betaUdDD computation)
    global ReUdD,ReUdDD
    ReUdD  = ixp.zerorank2(DIM)
    ReUdDD = ixp.zerorank3(DIM)
    for i in range(DIM):
        for j in range(DIM):
            ReUdD[i][j] = sp.diff(ReU[i], rfm.xx[j])
            for k in range(DIM):
                ReUdDD[i][j][k] = sp.diff(ReUdD[i][j], rfm.xx[k])

    # Step 3b: Compute 1st & 2nd derivatives of rescaling matrix.
    global ReDDdD,ReDDdDD
    ReDDdD = ixp.zerorank3(DIM)
    ReDDdDD = ixp.zerorank4(DIM)
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                ReDDdD[i][j][k] = (sp.diff(ReDD[i][j],rfm.xx[k]))
                for l in range(DIM):
                    # Simplifying this doesn't appear to help overall NRPy run time.
                    ReDDdDD[i][j][k][l] = sp.diff(ReDDdD[i][j][k],rfm.xx[l])

    # Step 3c: Compute 1st & 2nd derivatives of reference metric.
    global ghatDDdD,ghatDDdDD
    ghatDDdD = ixp.zerorank3(DIM)
    ghatDDdDD = ixp.zerorank4(DIM)
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                ghatDDdD[i][j][k] = sp.simplify(sp.diff(ghatDD[i][j],rfm.xx[k])) # FIXME: BAD: MUST BE SIMPLIFIED OR ANSWER IS INCORRECT! Must be some bug in sympy...
                for l in range(DIM):
                    ghatDDdDD[i][j][k][l] = (sp.diff(ghatDDdD[i][j][k],rfm.xx[l]))

    # Step 4a: Compute Christoffel symbols of reference metric.
    global GammahatUDD
    GammahatUDD = ixp.zerorank3(DIM)
    for i in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                for m in range(DIM):
                    GammahatUDD[i][k][l] += (sp.Rational(1,2))*ghatUU[i][m]*\
                                            (ghatDDdD[m][k][l] + ghatDDdD[m][l][k] - ghatDDdD[k][l][m])

    # Step 4b: Compute derivs of Christoffel symbols of reference metric.
    global GammahatUDDdD
    GammahatUDDdD = ixp.zerorank4(DIM)
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    GammahatUDDdD[i][j][k][l] = (sp.diff(GammahatUDD[i][j][k],rfm.xx[l]))