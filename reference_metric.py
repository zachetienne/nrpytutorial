# reference_metric.py: In terms of uniform
#     coordinate (xx[0],xx[1],xx[2]), define:
#     1) scalefactor_orthog:
#       orthogonal coordinate scale factor
#       (positive root of diagonal metric components),
#     2) xxCart[]: Cartesian coordinate (x,y,z)
#     3) xxSph[]: Spherical coordinate (r,theta,phi)
#

import time
import sympy as sp

import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri

# Step 0a: Initialize parameters
thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "CoordSystem", "Spherical"))

# Step 0b: Declare global variables
xx = gri.xx
xxCart = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s
Cart_to_xx = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s
Cartx,Carty,Cartz = sp.symbols("Cartx Carty Cartz", real=True)
Cart = [Cartx,Carty,Cartz]
xxSph  = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s
scalefactor_orthog = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s
have_already_called_reference_metric_function = False

def reference_metric(SymPySimplifyExpressions=True):
    global have_already_called_reference_metric_function # setting to global enables other modules to see updated value.
    have_already_called_reference_metric_function = True

    CoordSystem = par.parval_from_str("reference_metric::CoordSystem")
    M_PI = par.Cparameters("REAL",thismodule,"M_PI")

    global UnitVectors
    UnitVectors = ixp.zerorank2(3)

    # Set up hatted metric tensor, rescaling matrix, and rescaling vector
    if CoordSystem == "Spherical" or CoordSystem == "SinhSpherical" or CoordSystem == "SinhSphericalv2":
        xx[0] = sp.symbols("xx0", real=True)
        xx[1] = sp.symbols("xx1", real=True)

        r  = xx[0]
        th = xx[1]
        ph = xx[2]

        if CoordSystem == "Spherical":
            RMAX = par.Cparameters("REAL", thismodule, ["RMAX"])
            global xxmin
            global xxmax
            xxmin = [sp.sympify(0), sp.sympify(0), -M_PI]
            xxmax = [         RMAX,          M_PI,  M_PI]

            Cart_to_xx[0] = sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2)
            Cart_to_xx[1] = sp.acos(Cartz / Cart_to_xx[0])
            Cart_to_xx[2] = sp.atan2(Carty, Cartx)
        elif CoordSystem == "SinhSpherical" or CoordSystem == "SinhSphericalv2":
            AMPL, SINHW = par.Cparameters("REAL",thismodule,["AMPL","SINHW"])
    
            xxmin = [sp.sympify(0), sp.sympify(0), -M_PI]
            xxmax = [sp.sympify(1),          M_PI,  M_PI]
    
            # Set SinhSpherical radial coordinate by default; overwrite later if CoordSystem == "SinhSphericalv2".
            r = AMPL * (sp.exp(xx[0] / SINHW) - sp.exp(-xx[0] / SINHW)) / \
                       (sp.exp(1 / SINHW) - sp.exp(-1 / SINHW))

            Cart_to_xx[0] = SINHW*sp.asinh(sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2)*sp.sinh(1/SINHW)/AMPL)
            Cart_to_xx[1] = sp.acos(Cartz / sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2))
            Cart_to_xx[2] = sp.atan2(Carty, Cartx)

            # SinhSphericalv2 adds the parameter "const_dr", which allows for a region near xx[0]=0 to have
            # constant radial resolution of const_dr, provided the sinh() term does not dominate near xx[0]=0.
            if CoordSystem == "SinhSphericalv2":
                const_dr = par.Cparameters("REAL",thismodule,["const_dr"])
                r = AMPL*( const_dr*xx[0] + (sp.exp(xx[0] / SINHW) - sp.exp(-xx[0] / SINHW)) /
                           (sp.exp(1 / SINHW) - sp.exp(-1 / SINHW)) )

                Cart_to_xx[0] = "NewtonRaphson"
                Cart_to_xx[1] = sp.acos(Cartz / sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2))
                Cart_to_xx[2] = sp.atan2(Carty, Cartx)

        xxSph[0] = r
        xxSph[1] = th
        xxSph[2] = ph

        # Now define xCart, yCart, and zCart in terms of x0,xx[1],xx[2].
        #   Note that the relation between r and x0 is not necessarily trivial in SinhSpherical coordinates. See above.
        xxCart[0] = xxSph[0]*sp.sin(xxSph[1])*sp.cos(xxSph[2])
        xxCart[1] = xxSph[0]*sp.sin(xxSph[1])*sp.sin(xxSph[2])
        xxCart[2] = xxSph[0]*sp.cos(xxSph[1])

        scalefactor_orthog[0] = sp.diff(xxSph[0],xx[0])
        scalefactor_orthog[1] = xxSph[0]
        scalefactor_orthog[2] = xxSph[0]*sp.sin(xxSph[1])

        # Set the unit vectors
        UnitVectors = [[ sp.sin(xxSph[1])*sp.cos(xxSph[2]), sp.sin(xxSph[1])*sp.sin(xxSph[2]),  sp.cos(xxSph[1])],
                       [ sp.cos(xxSph[1])*sp.cos(xxSph[2]), sp.cos(xxSph[1])*sp.sin(xxSph[2]), -sp.sin(xxSph[1])],
                       [                 -sp.sin(xxSph[2]),                  sp.cos(xxSph[2]),  sp.sympify(0)   ]]
        
    elif CoordSystem == "NobleSphericalThetaOptionOne" or CoordSystem == "NobleSphericalThetaOptionTwo":
        R0,x0beg = par.Cparameters("REAL", thismodule, ["R0","x0beg"])
        xx[0] = sp.symbols("xx0", real=True)
        r  = R0 + sp.exp(x0beg + xx[0])

        th_c,xi,x1beg = par.Cparameters("REAL", thismodule, ["th_c","xi","x1beg"])
        xx[1] = sp.symbols("xx1", real=True)
        x1j = x1beg + xx[1]
        if CoordSystem == "NobleSphericalThetaOptionOne":
            th = th_c + (M_PI - 2*th_c)*x1j + xi*sp.sin(2*M_PI)*x1j
        elif CoordSystem == "NobleSphericalThetaOptionTwo":
            x1_n_exponent = par.Cparameters("REAL", thismodule, ["x1_n_exponent"])
            th = M_PI/2 * ( 1 + (1 - xi)*(2*x1j - 1) + (xi - 2*th_c/M_PI)*(2*x1j - 1)**x1_n_exponent )

        xx[2] = sp.symbols("xx2", real=True)
        ph = xx[2]

        # These DO NOT MATTER for interp_to_Sph_grids.
#         global xxmin
#         global xxmax
#         xxmin = [sp.sympify(0), sp.sympify(0), -M_PI]
#         xxmax = [         RMAX,          M_PI,  M_PI]

        xxSph[0] = r
        xxSph[1] = th
        xxSph[2] = ph

        Cart_to_xx[0] = sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2)
        Cart_to_xx[1] = sp.acos(Cartz / Cart_to_xx[0])
        Cart_to_xx[2] = sp.atan2(Carty, Cartx)
        
        # Now define xCart, yCart, and zCart in terms of x0,xx[1],xx[2].
        #   Note that the relation between r and x0 is not necessarily trivial in SinhSpherical coordinates. See above.
        xxCart[0] = xxSph[0]*sp.sin(xxSph[1])*sp.cos(xxSph[2])
        xxCart[1] = xxSph[0]*sp.sin(xxSph[1])*sp.sin(xxSph[2])
        xxCart[2] = xxSph[0]*sp.cos(xxSph[1])

        scalefactor_orthog[0] = sp.diff(xxSph[0],xx[0])
        scalefactor_orthog[1] = xxSph[0]
        scalefactor_orthog[2] = xxSph[0]*sp.sin(xxSph[1])

        # Set the unit vectors
        UnitVectors = [[ sp.sin(xxSph[1])*sp.cos(xxSph[2]), sp.sin(xxSph[1])*sp.sin(xxSph[2]),  sp.cos(xxSph[1])],
                       [ sp.cos(xxSph[1])*sp.cos(xxSph[2]), sp.cos(xxSph[1])*sp.sin(xxSph[2]), -sp.sin(xxSph[1])],
                       [                 -sp.sin(xxSph[2]),                  sp.cos(xxSph[2]),  sp.sympify(0)   ]]

    elif CoordSystem == "Cylindrical" or CoordSystem == "SinhCylindrical" or CoordSystem == "SinhCylindricalv2":
        # Assuming the cylindrical radial coordinate
        #   is positive makes nice simplifications of
        #   unit vectors possible.
        xx[0] = sp.symbols("xx0", real=True)

        RHOCYL = xx[0]
        PHICYL = xx[1]
        ZCYL   = xx[2]

        if CoordSystem == "Cylindrical":
            RHOMAX,ZMIN,ZMAX = par.Cparameters("REAL",thismodule,["RHOMAX","ZMIN","ZMAX"])
            xxmin = [sp.sympify(0), -M_PI, ZMIN]
            xxmax = [       RHOMAX,  M_PI, ZMAX]

            Cart_to_xx[0] = sp.sqrt(Cartx ** 2 + Carty ** 2)
            Cart_to_xx[1] = sp.atan2(Carty, Cartx)
            Cart_to_xx[2] = Cartz

        elif CoordSystem == "SinhCylindrical" or CoordSystem == "SinhCylindricalv2":
            AMPLRHO, SINHWRHO, AMPLZ, SINHWZ = par.Cparameters("REAL",thismodule,["AMPLRHO","SINHWRHO","AMPLZ","SINHWZ"])
    
            # Set SinhCylindrical radial & z coordinates by default; overwrite later if CoordSystem == "SinhCylindricalv2".
            RHOCYL = AMPLRHO * (sp.exp(xx[0] / SINHWRHO) - sp.exp(-xx[0] / SINHWRHO)) / (sp.exp(1 / SINHWRHO) - sp.exp(-1 / SINHWRHO))
            # phi coordinate remains unchanged.
            PHICYL = xx[1]
            ZCYL   = AMPLZ   * (sp.exp(xx[2] / SINHWZ)   - sp.exp(-xx[2] / SINHWZ))   / (sp.exp(1 / SINHWZ)   - sp.exp(-1 / SINHWZ))

            # SinhCylindricalv2 adds the parameters "const_drho", "const_dz", which allows for regions near xx[0]=0
            # and xx[2]=0 to have constant rho and z resolution of const_drho and const_dz, provided the sinh() terms
            # do not dominate near xx[0]=0 and xx[2]=0.
            if CoordSystem == "SinhCylindricalv2":
                const_drho, const_dz = par.Cparameters("REAL",thismodule,["const_drho","const_dz"])
 
                RHOCYL = AMPLRHO * ( const_drho*xx[0] + (sp.exp(xx[0] / SINHWRHO) - sp.exp(-xx[0] / SINHWRHO)) / (sp.exp(1 / SINHWRHO) - sp.exp(-1 / SINHWRHO)) )
                ZCYL   = AMPLZ   * ( const_dz  *xx[2] + (sp.exp(xx[2] / SINHWZ  ) - sp.exp(-xx[2] / SINHWZ  )) / (sp.exp(1 / SINHWZ  ) - sp.exp(-1 / SINHWZ  )) )
    
            xxmin = [sp.sympify(0), -M_PI, sp.sympify(-1)]
            xxmax = [sp.sympify(1),  M_PI, sp.sympify(+1)]

        xxCart[0] = RHOCYL*sp.cos(PHICYL)
        xxCart[1] = RHOCYL*sp.sin(PHICYL)
        xxCart[2] = ZCYL
    
        xxSph[0] = sp.sqrt(RHOCYL**2 + ZCYL**2)
        xxSph[1] = sp.acos(ZCYL / xxSph[0])
        xxSph[2] = PHICYL
    
        scalefactor_orthog[0] = sp.diff(RHOCYL,xx[0])
        scalefactor_orthog[1] = RHOCYL
        scalefactor_orthog[2] = sp.diff(ZCYL,xx[2])

        # Set the unit vectors
        UnitVectors = [[ sp.cos(PHICYL), sp.sin(PHICYL), sp.sympify(0)],
                       [-sp.sin(PHICYL), sp.cos(PHICYL), sp.sympify(0)],
                       [ sp.sympify(0),  sp.sympify(0),  sp.sympify(1)]]
        
    elif CoordSystem == "SymTP" or CoordSystem == "SinhSymTP":
        var1, var2= sp.symbols('var1 var2',real=True)
        bScale, AW, AA, AMAX, RHOMAX, ZMIN, ZMAX = par.Cparameters("REAL",thismodule,["bScale","AW","AA","AMAX","RHOMAX","ZMIN","ZMAX"])

        # Assuming xx0, xx1, and bScale
        #   are positive makes nice simplifications of
        #   unit vectors possible.
        xx[0],xx[1],bScale = sp.symbols("xx0 xx1 bScale", real=True)

        xxmin = ["0.0","0.0","0.0"]
        xxmax = ["params.AMAX","M_PI","2.0*M_PI"]
    
        AA = xx[0]
    
        if CoordSystem == "SinhSymTP":
            AA = (sp.exp(xx[0]/AW)-sp.exp(-xx[0]/AW))/2
    
        var1 = sp.sqrt(AA**2 + (bScale * sp.sin(xx[1]))**2)
        var2 = sp.sqrt(AA**2 + bScale**2)
    
        RHOSYMTP = AA*sp.sin(xx[1])
        PHSYMTP = xx[2]
        ZSYMTP = var2*sp.cos(xx[1])
    
        xxCart[0] = AA  *sp.sin(xx[1])*sp.cos(xx[2])
        xxCart[1] = AA  *sp.sin(xx[1])*sp.sin(xx[2])
        xxCart[2] = var2*sp.cos(xx[1])
    
        xxSph[0] = sp.sqrt(RHOSYMTP**2 + ZSYMTP**2)
        xxSph[1] = sp.acos(ZSYMTP / xxSph[0])
        xxSph[2] = PHSYMTP
    
        scalefactor_orthog[0] = sp.diff(AA,xx[0]) * var1 / var2
        scalefactor_orthog[1] = var1
        scalefactor_orthog[2] = AA * sp.sin(xx[1])

        # Set the transpose of the matrix of unit vectors
        UnitVectors = [[sp.sin(xx[1]) * sp.cos(xx[2]) * var2 / var1,
                        sp.sin(xx[1]) * sp.sin(xx[2]) * var2 / var1,
                        AA * sp.cos(xx[1]) / var1],
                       [AA * sp.cos(xx[1]) * sp.cos(xx[2]) / var1,
                        AA * sp.cos(xx[1]) * sp.sin(xx[2]) / var1,
                            -sp.sin(xx[1]) * var2 / var1],
                       [-sp.sin(xx[2]), sp.cos(xx[2]), 0]]

    elif CoordSystem == "Cartesian":
        xmin, xmax, ymin, ymax, zmin, zmax = par.Cparameters("REAL",thismodule,["xmin","xmax","ymin","ymax","zmin","zmax"])
        xxmin = ["xmin", "ymin", "zmin"]
        xxmax = ["xmax", "ymax", "zmax"]
    
        xxCart[0] = xx[0]
        xxCart[1] = xx[1]
        xxCart[2] = xx[2]

        xxSph[0] = sp.sqrt(xx[0] ** 2 + xx[1] ** 2 + xx[2] ** 2)
        xxSph[1] = sp.acos(xx[2] / xxSph[0])
        xxSph[2] = sp.atan2(xx[1], xx[0])

        scalefactor_orthog[0] = sp.sympify(1)
        scalefactor_orthog[1] = sp.sympify(1)
        scalefactor_orthog[2] = sp.sympify(1)

        # Set the transpose of the matrix of unit vectors
        UnitVectors = [[sp.sympify(1), sp.sympify(0), sp.sympify(0)],
                       [sp.sympify(0), sp.sympify(1), sp.sympify(0)],
                       [sp.sympify(0), sp.sympify(0), sp.sympify(1)]]
    else:
        print("CoordSystem == " + CoordSystem + " is not supported.")
        exit(1)

    # Finally, call ref_metric__hatted_quantities()
    #  to construct hatted metric, derivs of hatted
    #  metric, and Christoffel symbols
    ref_metric__hatted_quantities(SymPySimplifyExpressions)
    
def ref_metric__hatted_quantities(SymPySimplifyExpressions=True):
    # Step 0: Set dimension DIM
    DIM = par.parval_from_str("grid::DIM")

    global ReU,ReDD,ghatDD,ghatUU,detgammahat
    ReU    = ixp.zerorank1()
    ReDD   = ixp.zerorank2()
    ghatDD = ixp.zerorank2()

    # Step 1: Compute ghatDD (reference metric), ghatUU
    #         (inverse reference metric), as well as 
    #         rescaling vector ReU & rescaling matrix ReDD
    for i in range(DIM):
        scalefactor_orthog[i] = sp.sympify(scalefactor_orthog[i])
        ghatDD[i][i] = scalefactor_orthog[i]**2
        ReU[i] = 1/scalefactor_orthog[i]
        for j in range(DIM):
            ReDD[i][j] = scalefactor_orthog[i]*scalefactor_orthog[j]
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
        detgammahatdD[i] = (sp.diff(detgammahat, xx[i]))
        for j in range(DIM):
            detgammahatdDD[i][j] = sp.diff(detgammahatdD[i], xx[j])

    # Step 3a: Compute 1st & 2nd derivatives of rescaling vector.
    #          (E.g., needed in BSSN for betaUdDD computation)
    global ReUdD,ReUdDD
    ReUdD  = ixp.zerorank2(DIM)
    ReUdDD = ixp.zerorank3(DIM)
    for i in range(DIM):
        for j in range(DIM):
            ReUdD[i][j] = sp.diff(ReU[i], xx[j])
            for k in range(DIM):
                ReUdDD[i][j][k] = sp.diff(ReUdD[i][j], xx[k])

    # Step 3b: Compute 1st & 2nd derivatives of rescaling matrix.
    global ReDDdD,ReDDdDD
    ReDDdD = ixp.zerorank3(DIM)
    ReDDdDD = ixp.zerorank4(DIM)
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                ReDDdD[i][j][k] = (sp.diff(ReDD[i][j],xx[k]))
                for l in range(DIM):
                    # Simplifying this doesn't appear to help overall NRPy run time.
                    ReDDdDD[i][j][k][l] = sp.diff(ReDDdD[i][j][k],xx[l])

    # Step 3c: Compute 1st & 2nd derivatives of reference metric.
    global ghatDDdD,ghatDDdDD
    ghatDDdD = ixp.zerorank3(DIM)
    ghatDDdDD = ixp.zerorank4(DIM)
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                if SymPySimplifyExpressions==True:
                    ghatDDdD[i][j][k] = sp.simplify(sp.diff(ghatDD[i][j],xx[k])) # FIXME: BAD: MUST BE SIMPLIFIED OR ANSWER IS INCORRECT! Must be some bug in sympy...
                else:
                    ghatDDdD[i][j][k] = (sp.diff(ghatDD[i][j],xx[k])) # FIXME: BAD: MUST BE SIMPLIFIED OR ANSWER IS INCORRECT! Must be some bug in sympy...
                for l in range(DIM):
                    ghatDDdDD[i][j][k][l] = (sp.diff(ghatDDdD[i][j][k],xx[l]))

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
                    GammahatUDDdD[i][j][k][l] = (sp.diff(GammahatUDD[i][j][k],xx[l]))

# Compute proper distance in all 3 directions. Used to find the appropriate timestep for the CFL condition.
def ds_dirn(delxx):
    ds_dirn = ixp.zerorank1(3)
    for i in range(3):
        ds_dirn[i] = delxx[i]*scalefactor_orthog[i]
    return ds_dirn
