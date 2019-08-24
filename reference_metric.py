# reference_metric.py: Define all needed quantities
#     for a reference metric. 
# Given uniform (reference metric) coordinate 
#    (xx[0],xx[1],xx[2]), you must define:
#     1) xxmin[3],xxmax[3]: Valid ranges for each
#       uniform coordinate xx0,xx1,xx2
#     2) xxSph[3]: Spherical coordinate (r,theta,phi),
#       in terms of uniform coordinate xx0,xx1,xx2
#     3) xxCart[3]: Cartesian coordinate (x,y,z),
#       in terms of uniform coordinate xx0,xx1,xx2
#     4) scalefactor_orthog:
#       orthogonal coordinate scale factor
#       (positive root of diagonal reference metric 
#       components)
#     5) Cart_to_xx[3]: Inverse of xxCart:
#       xx0,xx1,xx2 as functions of (x,y,z). 
#       In the case that there exists no closed-form
#       expression, then a root finder might be needed
#     6) UnitVectors[3][3]: Unit vectors of reference 
#       metric.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

import time
import sys
import sympy as sp

import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
from outputC import superfast_uniq # contains superfast_uniq()

# Step 0a: Initialize parameters
thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "CoordSystem", "Spherical"))
par.initialize_param(par.glb_param("char", thismodule, "enable_rfm_precompute", "False"))
par.initialize_param(par.glb_param("char", thismodule, "rfm_precompute_Ccode_outdir", "Ccode"))

# Step 0b: Declare global variables
xx = gri.xx
xxCart = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s
Cart_to_xx = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s
Cartx,Carty,Cartz = sp.symbols("Cartx Carty Cartz", real=True)
Cart = [Cartx,Carty,Cartz]
xxSph  = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s
scalefactor_orthog = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s

scalefactor_orthog_funcform = ixp.zerorank1(DIM=4) # Must be set in terms of generic functions of xx[]s

have_already_called_reference_metric_function = False

def reference_metric(SymPySimplifyExpressions=True):
    global f0_of_xx0_funcform, f1_of_xx1_funcform, f2_of_xx0_xx1_funcform, f3_of_xx0_funcform
    global f0_of_xx0, f1_of_xx1, f2_of_xx0_xx1, f3_of_xx0
    f0_of_xx0_funcform     = sp.Function('f0_of_xx0_funcform')(xx[0])
    f1_of_xx1_funcform     = sp.Function('f1_of_xx1_funcform')(xx[1])
    f2_of_xx0_xx1_funcform = sp.Function('f2_of_xx0_xx1_funcform')(xx[0], xx[1])
    f3_of_xx0_funcform     = sp.Function('f3_of_xx0_funcform')(xx[0])
    f0_of_xx0, f1_of_xx1, f2_of_xx0_xx1, f3_of_xx0 = par.Cparameters("REAL", thismodule,
                                                          ["f0_of_xx0", "f1_of_xx1", "f2_of_xx0_xx1", "f3_of_xx0"], 1e300)
    # FIXME: Hack
    f0_of_xx0__D0, f0_of_xx0__DD00, f0_of_xx0__DDD000 = par.Cparameters("REAL", thismodule,
                                                                        ["f0_of_xx0__D0", "f0_of_xx0__DD00",
                                                                         "f0_of_xx0__DDD000"], 1e300)
    f1_of_xx1__D1, f1_of_xx1__DD11, f1_of_xx1__DDD111 = par.Cparameters("REAL", thismodule,
                                                                        ["f1_of_xx1__D1", "f1_of_xx1__DD11",
                                                                         "f1_of_xx1__DDD111"], 1e300)

    global have_already_called_reference_metric_function # setting to global enables other modules to see updated value.
    have_already_called_reference_metric_function = True

    CoordSystem = par.parval_from_str("reference_metric::CoordSystem")
    M_PI,M_SQRT1_2 = par.Cparameters("#define",thismodule,["M_PI","M_SQRT1_2"],"")

    global xxmin
    global xxmax

    global UnitVectors
    UnitVectors = ixp.zerorank2(DIM=3)

    # Set up hatted metric tensor, rescaling matrix, and rescaling vector
    
    #####################################################################
    # SPHERICAL-LIKE COORDINATE SYSTEMS WITH & WITHOUT RADIAL RESCALING #
    #####################################################################
    if CoordSystem == "Spherical" or CoordSystem == "SinhSpherical" or CoordSystem == "SinhSphericalv2":

        # Adding assumption real=True can help simplify expressions involving xx[0] & xx[1] below.
        xx[0] = sp.symbols("xx0", real=True)
        xx[1] = sp.symbols("xx1", real=True)

        if CoordSystem == "Spherical":
            RMAX = par.Cparameters("REAL", thismodule, ["RMAX"],10.0)
            xxmin = [sp.sympify(0), sp.sympify(0), -M_PI]
            xxmax = [         RMAX,          M_PI,  M_PI]

            r  = xx[0]
            th = xx[1]
            ph = xx[2]

            Cart_to_xx[0] = sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2)
            Cart_to_xx[1] = sp.acos(Cartz / Cart_to_xx[0])
            Cart_to_xx[2] = sp.atan2(Carty, Cartx)

        elif CoordSystem == "SinhSpherical":
            xxmin = [sp.sympify(0), sp.sympify(0), -M_PI]
            xxmax = [sp.sympify(1),          M_PI,  M_PI]
            
            AMPL, SINHW = par.Cparameters("REAL",thismodule,["AMPL","SINHW"],[10.0,0.2])
            # Set SinhSpherical radial coordinate by default; overwrite later if CoordSystem == "SinhSphericalv2".
            r = AMPL * (sp.exp(xx[0] / SINHW) - sp.exp(-xx[0] / SINHW)) / \
                       (sp.exp(1 / SINHW) - sp.exp(-1 / SINHW))
            th = xx[1]
            ph = xx[2]

            Cart_to_xx[0] = SINHW*sp.asinh(sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2)*sp.sinh(1/SINHW)/AMPL)
            Cart_to_xx[1] = sp.acos(Cartz / sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2))
            Cart_to_xx[2] = sp.atan2(Carty, Cartx)

        # SinhSphericalv2 adds the parameter "const_dr", which allows for a region near xx[0]=0 to have
        # constant radial resolution of const_dr, provided the sinh() term does not dominate near xx[0]=0.
        elif CoordSystem == "SinhSphericalv2":
            xxmin = [sp.sympify(0), sp.sympify(0), -M_PI]
            xxmax = [sp.sympify(1),          M_PI,  M_PI]
            
            AMPL, SINHW = par.Cparameters("REAL",thismodule,["AMPL","SINHW"],[10.0,0.2])
            const_dr = par.Cparameters("REAL",thismodule,["const_dr"],0.0625)
            r = AMPL*( const_dr*xx[0] + (sp.exp(xx[0] / SINHW) - sp.exp(-xx[0] / SINHW)) /
                       (sp.exp(1 / SINHW) - sp.exp(-1 / SINHW)) )
            th = xx[1]
            ph = xx[2]

            # NO CLOSED-FORM EXPRESSION FOR RADIAL INVERSION.
            # Cart_to_xx[0] = "NewtonRaphson"
            # Cart_to_xx[1] = sp.acos(Cartz / sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2))
            # Cart_to_xx[2] = sp.atan2(Carty, Cartx)
                

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

        f0_of_xx0 = xxSph[0]
        f1_of_xx1 = sp.sin(xxSph[1])
        scalefactor_orthog_funcform[0] = sp.diff(f0_of_xx0_funcform,xx[0])
        scalefactor_orthog_funcform[1] = f0_of_xx0_funcform
        scalefactor_orthog_funcform[2] = f0_of_xx0_funcform*f1_of_xx1_funcform
        
        # Set the unit vectors
        UnitVectors = [[ sp.sin(xxSph[1])*sp.cos(xxSph[2]), sp.sin(xxSph[1])*sp.sin(xxSph[2]),  sp.cos(xxSph[1])],
                       [ sp.cos(xxSph[1])*sp.cos(xxSph[2]), sp.cos(xxSph[1])*sp.sin(xxSph[2]), -sp.sin(xxSph[1])],
                       [                 -sp.sin(xxSph[2]),                  sp.cos(xxSph[2]),  sp.sympify(0)   ]]
        
    ######################################################################
    # SPHERICAL-LIKE COORDINATE SYSTEMS WITH RADIAL AND THETA RESCALINGS #
    ######################################################################
    elif CoordSystem == "NobleSphericalThetaOptionOne" or CoordSystem == "NobleSphericalThetaOptionTwo":
        # WARNING: CANNOT BE USED FOR SENR RUNS; 
        #  THESE DO NOT DEFINE xxmin, xxmax, Cart_to_xx
        #  ALSO THE RADIAL RESCALINGS ARE NOT ODD FUNCTIONS OF xx0,
        #  MEANING THAT CURVI. BOUNDARY CONDITIONS WILL NOT WORK.
        Rin,R0 = par.Cparameters("REAL", thismodule, ["Rin","R0"],[1.08986052555408,0.0])
        x0beg = sp.log(Rin-R0)
        xx[0] = sp.symbols("xx0", real=True)
        r  = R0 + sp.exp(x0beg + xx[0])

        # 0.053407075111026485 == 0.017*pi
        th_c,xi,x1beg = par.Cparameters("REAL", thismodule, ["th_c","xi","x1beg"],[0.053407075111026485,0.25,0.0])
        xx[1] = sp.symbols("xx1", real=True)
        x1j = x1beg + xx[1]
        if CoordSystem == "NobleSphericalThetaOptionOne":
            th = th_c + (M_PI - 2*th_c)*x1j + xi*sp.sin(2*M_PI*x1j)
        elif CoordSystem == "NobleSphericalThetaOptionTwo":
            x1_n_exponent = par.Cparameters("REAL", thismodule, ["x1_n_exponent"],9.0)
            th = M_PI/2 * ( 1 + (1 - xi)*(2*x1j - 1) + (xi - 2*th_c/M_PI)*(2*x1j - 1)**x1_n_exponent )

        xx[2] = sp.symbols("xx2", real=True)
        ph = xx[2]

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

    ##########################################################################
    # CYLINDRICAL-LIKE COORDINATE SYSTEMS WITH & WITHOUT RADIAL/Z RESCALINGS #
    ##########################################################################
    elif CoordSystem == "Cylindrical" or CoordSystem == "SinhCylindrical" or CoordSystem == "SinhCylindricalv2":
        # Assuming the cylindrical radial coordinate
        #   is positive makes nice simplifications of
        #   unit vectors possible.
        xx[0] = sp.symbols("xx0", real=True)

        if CoordSystem == "Cylindrical":
            RHOMAX,ZMIN,ZMAX = par.Cparameters("REAL",thismodule,["RHOMAX","ZMIN","ZMAX"],[10.0,-10.0,10.0])
            xxmin = [sp.sympify(0), -M_PI, ZMIN]
            xxmax = [       RHOMAX,  M_PI, ZMAX]

            RHOCYL = xx[0]
            PHICYL = xx[1]
            ZCYL   = xx[2]

            Cart_to_xx[0] = sp.sqrt(Cartx ** 2 + Carty ** 2)
            Cart_to_xx[1] = sp.atan2(Carty, Cartx)
            Cart_to_xx[2] = Cartz

        elif CoordSystem == "SinhCylindrical":
            xxmin = [sp.sympify(0), -M_PI, sp.sympify(-1)]
            xxmax = [sp.sympify(1),  M_PI, sp.sympify(+1)]

            AMPLRHO, SINHWRHO, AMPLZ, SINHWZ = par.Cparameters("REAL",thismodule,
                                                               ["AMPLRHO","SINHWRHO","AMPLZ","SINHWZ"],
                                                               [     10.0,       0.2,   10.0,    0.2])

            # Set SinhCylindrical radial & z coordinates by default; overwrite later if CoordSystem == "SinhCylindricalv2".
            RHOCYL = AMPLRHO * (sp.exp(xx[0] / SINHWRHO) - sp.exp(-xx[0] / SINHWRHO)) / (sp.exp(1 / SINHWRHO) - sp.exp(-1 / SINHWRHO))
            # phi coordinate remains unchanged.
            PHICYL = xx[1]
            ZCYL   = AMPLZ   * (sp.exp(xx[2] / SINHWZ)   - sp.exp(-xx[2] / SINHWZ))   / (sp.exp(1 / SINHWZ)   - sp.exp(-1 / SINHWZ))
            Cart_to_xx[0] = SINHWRHO*sp.asinh(sp.sqrt(Cartx ** 2 + Carty ** 2)*sp.sinh(1/SINHWRHO)/AMPLRHO)
            Cart_to_xx[1] = sp.atan2(Carty, Cartx)
            Cart_to_xx[2] = SINHWZ*sp.asinh(Cartz*sp.sinh(1/SINHWZ)/AMPLZ)
            
        # SinhCylindricalv2 adds the parameters "const_drho", "const_dz", which allows for regions near xx[0]=0
        # and xx[2]=0 to have constant rho and z resolution of const_drho and const_dz, provided the sinh() terms
        # do not dominate near xx[0]=0 and xx[2]=0.
        elif CoordSystem == "SinhCylindricalv2":
            xxmin = [sp.sympify(0), -M_PI, sp.sympify(-1)]
            xxmax = [sp.sympify(1),  M_PI, sp.sympify(+1)]
            AMPLRHO, SINHWRHO, AMPLZ, SINHWZ = par.Cparameters("REAL",thismodule,
                                                               ["AMPLRHO","SINHWRHO","AMPLZ","SINHWZ"],
                                                               [     10.0,       0.2,   10.0,    0.2])
            const_drho, const_dz = par.Cparameters("REAL",thismodule,["const_drho","const_dz"],[0.0625,0.0625])

            RHOCYL = AMPLRHO * ( const_drho*xx[0] + (sp.exp(xx[0] / SINHWRHO) - sp.exp(-xx[0] / SINHWRHO)) / (sp.exp(1 / SINHWRHO) - sp.exp(-1 / SINHWRHO)) )
            PHICYL = xx[1]
            ZCYL   = AMPLZ   * ( const_dz  *xx[2] + (sp.exp(xx[2] / SINHWZ  ) - sp.exp(-xx[2] / SINHWZ  )) / (sp.exp(1 / SINHWZ  ) - sp.exp(-1 / SINHWZ  )) )
    
            # NO CLOSED-FORM EXPRESSION FOR RADIAL OR Z INVERSION.
            # Cart_to_xx[0] = "NewtonRaphson"
            # Cart_to_xx[1] = sp.atan2(Carty, Cartx)
            # Cart_to_xx[2] = "NewtonRaphson"

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
        bScale, AW, AMAX, RHOMAX, ZMIN, ZMAX = par.Cparameters("REAL",thismodule,
                                                               ["bScale","AW","AMAX","RHOMAX","ZMIN","ZMAX"],
                                                               [0.5,     0.2,   10.0,    10.0, -10.0,  10.0])

        # Assuming xx0, xx1, and bScale
        #   are positive makes nice simplifications of
        #   unit vectors possible.
        xx[0],xx[1] = sp.symbols("xx0 xx1", real=True)

        xxmin = [sp.sympify(0), sp.sympify(0),-M_PI]
        xxmax = [         AMAX,          M_PI, M_PI]
    
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
        xxCart[2] = ZSYMTP
    
        xxSph[0] = sp.sqrt(RHOSYMTP**2 + ZSYMTP**2)
        xxSph[1] = sp.acos(ZSYMTP / xxSph[0])
        xxSph[2] = PHSYMTP

        if CoordSystem == "SymTP":
            rSph  = sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2)
            thSph = sp.acos(Cartz / rSph)
            phSph = sp.atan2(Carty, Cartx)

            # Mathematica script to compute Cart_to_xx[]
#             AA = x1;
#             var2 = Sqrt[AA^2 + bScale^2];
#             RHOSYMTP = AA*Sin[x2];
#             ZSYMTP = var2*Cos[x2];
#             Solve[{rSph == Sqrt[RHOSYMTP^2 + ZSYMTP^2],
#                    thSph == ArcCos[ZSYMTP/Sqrt[RHOSYMTP^2 + ZSYMTP^2]],
#                    phSph == x3}, 
#                   {x1, x2, x3}]
            Cart_to_xx[0] = sp.sqrt(-bScale**2 + rSph**2 + 
                                    sp.sqrt(bScale**4 + 2*bScale**2*rSph**2 + rSph**4 - 
                                            4*bScale**2*rSph**2*sp.cos(thSph)**2))*M_SQRT1_2 # M_SQRT1_2 = 1/sqrt(2); define this way for UnitTesting

            # The sign() function in the following expression ensures the correct root is taken.
            Cart_to_xx[1] = sp.acos(sp.sign(Cartz)*(
                                      sp.sqrt(1 + rSph**2/bScale**2 - 
                                              sp.sqrt(bScale**4 + 2*bScale**2*rSph**2 + rSph**4 - 
                                                      4*bScale**2*rSph**2*sp.cos(thSph)**2)/bScale**2)*M_SQRT1_2)) # M_SQRT1_2 = 1/sqrt(2); define this way for UnitTesting

            Cart_to_xx[2] = phSph

        elif CoordSystem == "SinhSymTP":
            pass
            # Closed form expression for Cart_to_xx in SinhSymTP may exist, but has not yet been found

        scalefactor_orthog[0] = sp.diff(AA,xx[0]) * var1 / var2
        scalefactor_orthog[1] = var1
        scalefactor_orthog[2] = AA * sp.sin(xx[1])

        f0_of_xx0              = AA
        f1_of_xx1              = sp.sin(xxSph[1])
        f2_of_xx0_xx1_funcform = var1
        f3_of_xx0              = var2

        scalefactor_orthog_funcform[0] = sp.diff(f0_of_xx0_funcform,xx[0]) * f2_of_xx0_xx1_funcform/f3_of_xx0
        scalefactor_orthog_funcform[1] = f2_of_xx0_xx1_funcform
        scalefactor_orthog_funcform[2] = f0_of_xx0_funcform*f1_of_xx1_funcform

        # Set the transpose of the matrix of unit vectors
        UnitVectors = [[sp.sin(xx[1]) * sp.cos(xx[2]) * var2 / var1,
                        sp.sin(xx[1]) * sp.sin(xx[2]) * var2 / var1,
                        AA * sp.cos(xx[1]) / var1],
                       [AA * sp.cos(xx[1]) * sp.cos(xx[2]) / var1,
                        AA * sp.cos(xx[1]) * sp.sin(xx[2]) / var1,
                            -sp.sin(xx[1]) * var2 / var1],
                       [-sp.sin(xx[2]), sp.cos(xx[2]), sp.sympify(0)]]

    elif CoordSystem == "Cartesian":
        xmin, xmax, ymin, ymax, zmin, zmax = par.Cparameters("REAL",thismodule,
                                                             ["xmin","xmax","ymin","ymax","zmin","zmax"],
                                                             [ -10.0,  10.0, -10.0,  10.0, -10.0,  10.0])
        xxmin = ["xmin", "ymin", "zmin"]
        xxmax = ["xmax", "ymax", "zmax"]
    
        xxCart[0] = xx[0]
        xxCart[1] = xx[1]
        xxCart[2] = xx[2]

        xxSph[0] = sp.sqrt(xx[0] ** 2 + xx[1] ** 2 + xx[2] ** 2)
        xxSph[1] = sp.acos(xx[2] / xxSph[0])
        xxSph[2] = sp.atan2(xx[1], xx[0])
        
        Cart_to_xx[0] = Cartx
        Cart_to_xx[1] = Carty
        Cart_to_xx[2] = Cartz

        scalefactor_orthog[0] = sp.sympify(1)
        scalefactor_orthog[1] = sp.sympify(1)
        scalefactor_orthog[2] = sp.sympify(1)

        # Set the transpose of the matrix of unit vectors
        UnitVectors = [[sp.sympify(1), sp.sympify(0), sp.sympify(0)],
                       [sp.sympify(0), sp.sympify(1), sp.sympify(0)],
                       [sp.sympify(0), sp.sympify(0), sp.sympify(1)]]

    else:
        print("CoordSystem == " + CoordSystem + " is not supported.")
        sys.exit(1)

    # Finally, call ref_metric__hatted_quantities()
    #  to construct hatted metric, derivs of hatted
    #  metric, and Christoffel symbols
    ref_metric__hatted_quantities(SymPySimplifyExpressions)
    # ref_metric__hatted_quantities(scalefactor_orthog_funcform,SymPySimplifyExpressions)
    # ref_metric__hatted_quantities(scalefactor_orthog,SymPySimplifyExpressions)

def ref_metric__hatted_quantities(SymPySimplifyExpressions=True):

    enable_rfm_precompute = False
    if par.parval_from_str(thismodule+"::enable_rfm_precompute") == "True":
        enable_rfm_precompute = True

    # Step 0: Set dimension DIM
    DIM = par.parval_from_str("grid::DIM")

    global ReU,ReDD,ghatDD,ghatUU,detgammahat
    ReU    = ixp.zerorank1()
    ReDD   = ixp.zerorank2()
    ghatDD = ixp.zerorank2()

    # Step 1: Compute ghatDD (reference metric), ghatUU
    #         (inverse reference metric), as well as 
    #         rescaling vector ReU & rescaling matrix ReDD
    if enable_rfm_precompute == False:
        for i in range(DIM):
            scalefactor_orthog[i] = sp.sympify(scalefactor_orthog[i])
            ghatDD[i][i] = scalefactor_orthog[i]**2
            ReU[i] = 1/scalefactor_orthog[i]
            for j in range(DIM):
                ReDD[i][j] = scalefactor_orthog[i]*scalefactor_orthog[j]
    else:
        for i in range(DIM):
            scalefactor_orthog_funcform[i] = sp.sympify(scalefactor_orthog_funcform[i])
            ghatDD[i][i] = scalefactor_orthog_funcform[i]**2
            ReU[i] = 1/scalefactor_orthog_funcform[i]
            for j in range(DIM):
                ReDD[i][j] = scalefactor_orthog_funcform[i]*scalefactor_orthog_funcform[j]

    # Step 1b: Compute ghatUU
    ghatUU, detgammahat = ixp.symm_matrix_inverter3x3(ghatDD)

    # Step 1c: Sanity check: verify that ReDD, ghatDD, 
    #          and ghatUU are all symmetric rank-2: 
    for i in range(DIM):
        for j in range(DIM):
            if ReDD[i][j] != ReDD[j][i]:
                print("Error: ReDD["+ str(i) + "][" + str(j) + "] != ReDD["+ str(j) + "][" + str(i) + ": " + str(ReDD[i][j]) + "!=" + str(ReDD[j][i]))
                sys.exit(1)
            if ghatDD[i][j] != ghatDD[j][i]:
                print("Error: ghatDD["+ str(i) + "][" + str(j) + "] != ghatDD["+ str(j) + "][" + str(i) + ": " + str(ghatDD[i][j]) + "!=" + str(ghatDD[j][i]))
                sys.exit(1)
            if ghatUU[i][j] != ghatUU[j][i]:
                print("Error: ghatUU["+ str(i) + "][" + str(j) + "] != ghatUU["+ str(j) + "][" + str(i) + ": " + str(ghatUU[i][j]) + "!=" + str(ghatUU[j][i]))
                sys.exit(1)

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
#                    ghatDDdD[i][j][k] = sp.trigsimp(sp.diff(ghatDD[i][j],xx[k])) # FIXME: BAD: MUST BE SIMPLIFIED OR ANSWER IS INCORRECT! Must be some bug in sympy...
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
#                    GammahatUDD[i][k][l] += sp.trigsimp((sp.Rational(1,2))*ghatUU[i][m]*\
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

    # Step 4c: If rfm_precompute is disabled, then we are finished with this function.
    #          Otherwise continue to Step 5.
    if enable_rfm_precompute == False:
        return
    else:
        CoordSystem = par.parval_from_str("reference_metric::CoordSystem")
        # if not (("Spherical" in CoordSystem) or ("SymTP" in CoordSystem)):
        if not (( "Spherical" in CoordSystem)):
            print("Error: CoordSystem == "+CoordSystem+" does not yet support rfm precompute infrastructure.")
            sys.exit(1)

        # enable_rfm_precompute: precompute and store in memory complicated
        #     expressions related to the reference metric (a.k.a., "hatted
        #     quantities")

        # The precomputed "hatted quantity" expressions will be stored in
        #    a C struct called rfmstruct. As these expressions generally
        #    involve computationally expensive transcendental functions
        #    of xx0,xx1,or xx2, and xx0,xx1, and xx2 remain fixed across
        #    most (if not all) of a given simulation, setting up the
        #    rfmstruct can greatly improve performance.

        # The core challenge in setting up the rfmstruct is collecting
        #    all the information needed to automatically generate it.
        # Step 5 and onwards implements this algorithm, using the
        #    *generic functional form* of the hatted quantities (as
        #    opposed to the exact closed-form expressions of the
        #    hatted quantities) computed above.

        # Step 5: Now that all hatted quantities are written in terms of generic SymPy functions,
        #         we will now replace SymPy functions with simple variables using rigid NRPy+ syntax,
        #         and store these variables to globals defined above.
        def make_replacements(input):
            for i in ["0", "1", "2"]:
                inputnew = input.replace(
                    "Derivative(f" + i + "_of_xx" + i + "_funcform(xx" + i + "), (xx" + i + ", 3))",
                    "Derivative(f" + i + "_of_xx" + i + "_funcform(xx" + i + "), xx" + i + ", xx" + i + ", xx" + i + ")"). \
                    replace("Derivative(f" + i + "_of_xx" + i + "_funcform(xx" + i + "), (xx" + i + ", 2))",
                            "Derivative(f" + i + "_of_xx" + i + "_funcform(xx" + i + "), xx" + i + ", xx" + i + ")"). \
                    replace(", xx" + i + ", xx" + i + ", xx" + i + ")", "__DDD" + i + i + i). \
                    replace(", xx" + i + ", xx" + i + ")", "__DD" + i + i). \
                    replace(", xx" + i + ")", "__D" + i). \
                    replace("f" + i + "_of_xx" + i + "_funcform(xx" + i + ")", "f" + i + "_of_xx" + i)
                input = inputnew
            inputnew = input.replace("Derivative(", "")
            if "Derivative" in inputnew:
                print("Error: ", inputnew)
                sys.exit(1)
            input = inputnew
            return input
        #
        detgammahat = sp.sympify(make_replacements(str(detgammahat)))
        for i in range(DIM):
            ReU[i] = sp.sympify(make_replacements(str(ReU[i])))
            detgammahatdD[i] = sp.sympify(make_replacements(str(detgammahatdD[i])))
            for j in range(DIM):
                ReDD[i][j] = sp.sympify(make_replacements(str(ReDD[i][j])))
                ReUdD[i][j] = sp.sympify(make_replacements(str(ReUdD[i][j])))
                ghatDD[i][j] = sp.sympify(make_replacements(str(ghatDD[i][j])))
                ghatUU[i][j] = sp.sympify(make_replacements(str(ghatUU[i][j])))
                detgammahatdDD[i][j] = sp.sympify(make_replacements(str(detgammahatdDD[i][j])))
                for k in range(DIM):
                    ReDDdD[i][j][k] = sp.sympify(make_replacements(str(ReDDdD[i][j][k])))
                    ReUdDD[i][j][k] = sp.sympify(make_replacements(str(ReUdDD[i][j][k])))
                    ghatDDdD[i][j][k] = sp.sympify(make_replacements(str(ghatDDdD[i][j][k])))
                    GammahatUDD[i][j][k] = sp.sympify(make_replacements(str(GammahatUDD[i][j][k])))
                    for l in range(DIM):
                        ReDDdDD[i][j][k][l] = sp.sympify(make_replacements(str(ReDDdDD[i][j][k][l])))
                        ghatDDdDD[i][j][k][l] = sp.sympify(make_replacements(str(ghatDDdDD[i][j][k][l])))
                        GammahatUDDdD[i][j][k][l] = sp.sympify(make_replacements(str(GammahatUDDdD[i][j][k][l])))

        # Step 6: At this point, each expression is written in terms of the generic functions
        #         of xx0, xx1, and/or xx2 and their derivatives. Depending on the functions, some
        #         of these derivatives may be zero. In Step 5 we'll evaluate the function
        #         derivatives exactly and set the expressions to zero. Otherwise in the C code
        #         we'd be storing performing arithmetic with zeros -- wasteful!

        # Step 6.a: Construct the full list of *unique* NRPy+ variables representing the
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

        # Step 6.b: Using the expressions f?_of_xx? set in reference_metric(),
        #           evaluate each needed derivative and, in the case it is zero,
        #           set the corresponding "freevar" variable to zero.
        freevars_uniq_vals = []
        for i in range(len(freevars_uniq)):
            var = freevars_uniq[i]
            basename = str(var).split("__")[0].replace("_funcform", "")
            derivatv = ""
            if "__" in str(var):
                derivatv = str(var).split("__")[1].replace("_funcform", "")
            if basename == "f0_of_xx0":
                basefunc = f0_of_xx0
            elif basename == "f1_of_xx1":
                basefunc = f1_of_xx1
            else:
                print("Error: function inside " + str(var) + " undefined.")
                sys.exit(1)
            diff_result = basefunc
            if derivatv == "":
                pass
            else:
                derivorder = derivatv.replace("d", "").replace("D", "").replace("0", "0 ").replace("1", "1 ").replace(
                    "2", "2 ").split(" ")
                for derivdirn in derivorder:
                    if derivdirn != "":
                        derivwrt = xx[int(derivdirn)]
                        diff_result = sp.diff(diff_result, derivwrt)
            freevars_uniq_vals.append(diff_result)
            if diff_result == sp.sympify(0):
                freevars_uniq_zeroed[i] = 0

        # Step 6.c: Finally, substitute zero for all functions & derivatives that evaluate to zero.
        for varidx in range(len(freevars_uniq)):
            detgammahat = detgammahat.subs(freevars_uniq[varidx], freevars_uniq_zeroed[varidx])
            for i in range(DIM):
                ReU[i] = ReU[i].subs(freevars_uniq[varidx], freevars_uniq_zeroed[varidx])
                detgammahatdD[i] = detgammahatdD[i].subs(freevars_uniq[varidx], freevars_uniq_zeroed[varidx])
                for j in range(DIM):
                    ReDD[i][j] = ReDD[i][j].subs(freevars_uniq[varidx], freevars_uniq_zeroed[varidx])
                    ReUdD[i][j] = ReUdD[i][j].subs(freevars_uniq[varidx], freevars_uniq_zeroed[varidx])
                    ghatDD[i][j] = ghatDD[i][j].subs(freevars_uniq[varidx], freevars_uniq_zeroed[varidx])
                    ghatUU[i][j] = ghatUU[i][j].subs(freevars_uniq[varidx], freevars_uniq_zeroed[varidx])
                    # print(varidx,i,j,freevars_uniq[varidx],freevars_uniq_zeroed[varidx],detgammahatdDD[i][j])
                    detgammahatdDD[i][j] = detgammahatdDD[i][j].subs(freevars_uniq[varidx],
                                                                     freevars_uniq_zeroed[varidx])
                    # print(varidx,i,j,freevars_uniq[varidx],freevars_uniq_zeroed[varidx],detgammahatdDD[i][j])
                    for k in range(DIM):
                        ReDDdD[i][j][k] = ReDDdD[i][j][k].subs(freevars_uniq[varidx], freevars_uniq_zeroed[varidx])
                        ReUdDD[i][j][k] = ReUdDD[i][j][k].subs(freevars_uniq[varidx], freevars_uniq_zeroed[varidx])
                        ghatDDdD[i][j][k] = ghatDDdD[i][j][k].subs(freevars_uniq[varidx], freevars_uniq_zeroed[varidx])
                        GammahatUDD[i][j][k] = GammahatUDD[i][j][k].subs(freevars_uniq[varidx],
                                                                         freevars_uniq_zeroed[varidx])
                        for l in range(DIM):
                            ReDDdDD[i][j][k][l] = ReDDdDD[i][j][k][l].subs(freevars_uniq[varidx],
                                                                           freevars_uniq_zeroed[varidx])
                            ghatDDdDD[i][j][k][l] = ghatDDdDD[i][j][k][l].subs(freevars_uniq[varidx],
                                                                               freevars_uniq_zeroed[varidx])
                            GammahatUDDdD[i][j][k][l] = GammahatUDDdD[i][j][k][l].subs(freevars_uniq[varidx],
                                                                                       freevars_uniq_zeroed[varidx])

        # Step 7: Construct needed C code for declaring rfmstruct, allocating storage for
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
        readvr_str = ["", "", ""]
        readvr_SIMD_outer_str = ["", "", ""]
        readvr_SIMD_inner_str = ["", "", ""]
        freemm_str = ""
        for dirn in [0, 1, 2]:
            malloc_size = gri.Nxx_plus_2NGHOSTS[dirn]
            #        malloc_size = "Nxx_plus_2NGHOSTS["+str(dirn)+"]"

            numvars = 0
            for varidx in range(len(freevars_uniq)):
                if "xx" + str(dirn) in str(freevars_uniq_zeroed[varidx]):
                    numvars = numvars + 1
                    struct_str += "\tREAL *restrict " + str(freevars_uniq_zeroed[varidx]) + ";\n"
                    malloc_str += "rfmstruct." + str(
                        freevars_uniq_zeroed[varidx]) + " = (REAL *)malloc(sizeof(REAL)*" + str(malloc_size) + ");\n"
                    freemm_str += "free(rfmstruct." + str(freevars_uniq_zeroed[varidx]) + ");\n"

                    define_str += """
        for(int ii=0;ii<""" + str(malloc_size) + """;ii++) {
            const REAL xx""" + str(dirn) + """ = xx[""" + str(dirn) + """][ii];
            rfmstruct.""" + str(freevars_uniq_zeroed[varidx]) + """[ii] = """ + str(
                        sp.ccode(freevars_uniq_vals[varidx])) + """;
        }"""
                    readvr_str[dirn] += "const REAL " + str(freevars_uniq_zeroed[varidx]) + " = rfmstruct->" + str(
                        freevars_uniq_zeroed[varidx]) + "[i" + str(dirn) + "];\n"
                    readvr_SIMD_outer_str[dirn] += "const double NOSIMD" + str(
                        freevars_uniq_zeroed[varidx]) + " = rfmstruct->" + str(
                        freevars_uniq_zeroed[varidx]) + "[i" + str(dirn) + "]; "
                    readvr_SIMD_outer_str[dirn] += "const REAL_SIMD_ARRAY " + str(
                        freevars_uniq_zeroed[varidx]) + " = ConstSIMD(NOSIMD" + str(
                        freevars_uniq_zeroed[varidx]) + ");\n"
                    readvr_SIMD_inner_str[dirn] += "const REAL_SIMD_ARRAY " + str(
                        freevars_uniq_zeroed[varidx]) + " = ReadSIMD(&rfmstruct->" + str(
                        freevars_uniq_zeroed[varidx]) + "[i" + str(dirn) + "]);\n"

        struct_str += "} rfm_struct;\n\n"

        # Step 8: Output needed C code to files
        outdir = par.parval_from_str(thismodule+"::rfm_precompute_Ccode_outdir")
        with open(outdir + "/rfm_struct__declare.h", "w") as file:
            file.write(struct_str)
        with open(outdir + "/rfm_struct__malloc.h", "w") as file:
            file.write(malloc_str)
        with open(outdir + "/rfm_struct__define.h", "w") as file:
            file.write(define_str)
        for i in range(3):
            with open(outdir + "/rfm_struct__read" + str(i) + ".h", "w") as file:
                file.write(readvr_str[i])
            with open(outdir + "/rfm_struct__SIMD_outer_read" + str(i) + ".h", "w") as file:
                file.write(readvr_SIMD_outer_str[i])
            with open(outdir + "/rfm_struct__SIMD_inner_read" + str(i) + ".h", "w") as file:
                file.write(readvr_SIMD_inner_str[i])
        with open(outdir + "/rfm_struct__freemem.h", "w") as file:
            file.write(freemm_str)


# Compute proper distance in all 3 directions. Used to find the appropriate timestep for the CFL condition.
def ds_dirn(delxx):
    ds_dirn = ixp.zerorank1(3)
    for i in range(3):
        ds_dirn[i] = delxx[i]*scalefactor_orthog[i]
    return ds_dirn
