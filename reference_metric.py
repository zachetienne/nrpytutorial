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

import sympy as sp                         # Import SymPy
from outputC import outputC,superfast_uniq,outC_function_dict,add_to_Cfunction_dict # NRPy+: Core C code output module
import NRPy_param_funcs as par             # NRPy+: Parameter interface
import grid as gri                         # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp                   # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import sys                                 # Standard Python modules for multiplatform OS-level functions

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
    global f0_of_xx0_funcform, f1_of_xx1_funcform, f2_of_xx0_xx1_funcform, f3_of_xx0_funcform, f4_of_xx2_funcform
    global f0_of_xx0, f1_of_xx1, f2_of_xx0_xx1, f3_of_xx0, f4_of_xx2
    f0_of_xx0_funcform     = sp.Function('f0_of_xx0_funcform')(xx[0])
    f1_of_xx1_funcform     = sp.Function('f1_of_xx1_funcform')(xx[1])
    f2_of_xx0_xx1_funcform = sp.Function('f2_of_xx0_xx1_funcform')(xx[0], xx[1])
    f3_of_xx0_funcform     = sp.Function('f3_of_xx0_funcform')(xx[0])
    f4_of_xx2_funcform     = sp.Function('f4_of_xx2_funcform')(xx[2])
    f0_of_xx0, f1_of_xx1, f2_of_xx0_xx1, f3_of_xx0, f4_of_xx2 = par.Cparameters("REAL", thismodule,
                                  ["f0_of_xx0", "f1_of_xx1", "f2_of_xx0_xx1", "f3_of_xx0", "f4_of_xx2"], 1e300)
    # FIXME: Hack
    f0_of_xx0__D0, f0_of_xx0__DD00, f0_of_xx0__DDD000 = par.Cparameters("REAL", thismodule,
                                                                        ["f0_of_xx0__D0", "f0_of_xx0__DD00",
                                                                         "f0_of_xx0__DDD000"], 1e300)
    f1_of_xx1__D1, f1_of_xx1__DD11, f1_of_xx1__DDD111 = par.Cparameters("REAL", thismodule,
                                                                        ["f1_of_xx1__D1", "f1_of_xx1__DD11",
                                                                         "f1_of_xx1__DDD111"], 1e300)
    f2_of_xx0_xx1__D0,f2_of_xx0_xx1__D1,f2_of_xx0_xx1__DD00,f2_of_xx0_xx1__DD11 = \
        par.Cparameters("REAL", thismodule,
                        ["f2_of_xx0_xx1__D0","f2_of_xx0_xx1__D1","f2_of_xx0_xx1__DD00","f2_of_xx0_xx1__DD11"],
                        1e300)
    f3_of_xx0__D0,f3_of_xx0__DD00     = par.Cparameters("REAL", thismodule,["f3_of_xx0__D0","f3_of_xx0__DD00"], 1e300)
    f4_of_xx2__D2,f4_of_xx2__DD22     = par.Cparameters("REAL", thismodule,["f4_of_xx2__D2","f4_of_xx2__DD22"], 1e300)

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

        f0_of_xx0              = RHOCYL
        f4_of_xx2              = sp.diff(ZCYL,xx[2])

        scalefactor_orthog_funcform[0] = sp.diff(f0_of_xx0_funcform,xx[0])
        scalefactor_orthog_funcform[1] = f0_of_xx0_funcform
        scalefactor_orthog_funcform[2] = f4_of_xx2_funcform

        # Set the unit vectors
        UnitVectors = [[ sp.cos(PHICYL), sp.sin(PHICYL), sp.sympify(0)],
                       [-sp.sin(PHICYL), sp.cos(PHICYL), sp.sympify(0)],
                       [ sp.sympify(0),  sp.sympify(0),  sp.sympify(1)]]

    elif CoordSystem == "SymTP" or CoordSystem == "SinhSymTP":
        # var1, var2= sp.symbols('var1 var2',real=True)
        bScale, SINHWAA, AMAX = par.Cparameters("REAL",thismodule,
                                                ["bScale","SINHWAA","AMAX"],
                                                [0.5,     0.2,      10.0  ])

        # Assuming xx0, xx1, and bScale
        #   are positive makes nice simplifications of
        #   unit vectors possible.
        xx[0],xx[1] = sp.symbols("xx0 xx1", real=True)

        xxmin = [sp.sympify(0), sp.sympify(0),-M_PI]
        xxmax = [         AMAX,          M_PI, M_PI]

        AA = xx[0]

        if CoordSystem == "SinhSymTP":
            # With xxmax[0] == AMAX, sinh(xx0/AMAX) will evaluate to a number between 0 and 1.
            #   Similarly, sinh(xx0/(AMAX*SINHWAA)) / sinh(1/SINHWAA) will also evaluate to a number between 0 and 1.
            #   Then AA = AMAX*sinh(xx0/(AMAX*SINHWAA)) / sinh(1/SINHWAA) will evaluate to a number between 0 and AMAX.
            AA = AMAX * (sp.exp(xx[0] / (AMAX*SINHWAA)) - sp.exp(-xx[0] / (AMAX*SINHWAA))) / (sp.exp(1 / SINHWAA) - sp.exp(-1 / AMAX))

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
        f1_of_xx1              = sp.sin(xx[1])
        f2_of_xx0_xx1          = var1
        f3_of_xx0              = var2

        scalefactor_orthog_funcform[0] = sp.diff(f0_of_xx0_funcform,xx[0]) * f2_of_xx0_xx1_funcform/f3_of_xx0_funcform
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

        scalefactor_orthog_funcform[0] = sp.sympify(1)
        scalefactor_orthog_funcform[1] = sp.sympify(1)
        scalefactor_orthog_funcform[2] = sp.sympify(1)

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
        # if not (( "Spherical" in CoordSystem)):
        #     print("Error: CoordSystem == "+CoordSystem+" does not yet support rfm precompute infrastructure.")
        #     sys.exit(1)

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
        def make_replacements(expr):
            sympy_version = sp.__version__.replace("rc","...").replace("b","...") # Ignore the rc's and b's for release candidates & betas.
            sympy_major_version = int(sympy_version.split(".")[0])
            sympy_minor_version = int(sympy_version.split(".")[1])
            if sympy_major_version < 1 or (sympy_major_version >= 1 and sympy_minor_version < 2):
                print("Detected SymPy version "+sympy_version)
                print("Sorry, reference metric precomputation unsupported in SymPy < 1.2!")
                sys.exit(1)

            for item in sp.preorder_traversal(expr):
                if item.func == sp.Derivative:
                    stringfunc = str(item.args[0]).split("_funcform(", 1)[0]  # store everything before _funcform(...
                    stringderv = str(item.args[1]).replace(" ", "")  # Ignore whitespace
                    deriv_wrt = stringderv.split(",")[0].replace("(xx", "")
                    derivorder = int(stringderv.split(",")[1].replace(")", ""))

                    derivop = "__D"
                    for i in range(derivorder - 1):
                        derivop += "D"
                    derivop += deriv_wrt
                    for i in range(derivorder - 1):
                        derivop += deriv_wrt
                    expr = expr.xreplace(
                        {item: sp.sympify(stringfunc + derivop)})

            for item in sp.preorder_traversal(expr):
                if "_funcform" in str(item.func):
                    stringfunc = str(item.func).split("_funcform", 1)[0]  # store everything before _funcform(...
                    expr = expr.xreplace({item: sp.sympify(stringfunc)})
            return expr

        detgammahat = make_replacements(detgammahat)
        for i in range(DIM):
            ReU[i] = make_replacements(ReU[i])
            detgammahatdD[i] = make_replacements(detgammahatdD[i])
            for j in range(DIM):
                ReDD[i][j] = make_replacements(ReDD[i][j])
                ReUdD[i][j] = make_replacements(ReUdD[i][j])
                ghatDD[i][j] = make_replacements(ghatDD[i][j])
                ghatUU[i][j] = make_replacements(ghatUU[i][j])
                detgammahatdDD[i][j] = make_replacements(detgammahatdDD[i][j])
                for k in range(DIM):
                    ReDDdD[i][j][k] = make_replacements(ReDDdD[i][j][k])
                    ReUdDD[i][j][k] = make_replacements(ReUdDD[i][j][k])
                    ghatDDdD[i][j][k] = make_replacements(ghatDDdD[i][j][k])
                    GammahatUDD[i][j][k] = make_replacements(GammahatUDD[i][j][k])
                    for l in range(DIM):
                        ReDDdDD[i][j][k][l] = make_replacements(ReDDdDD[i][j][k][l])
                        ghatDDdDD[i][j][k][l] = make_replacements(ghatDDdDD[i][j][k][l])
                        GammahatUDDdD[i][j][k][l] = make_replacements(GammahatUDDdD[i][j][k][l])

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

        freevars_uniq_xx_indep = []
        for i in range(len(freevars_uniq)):
            freevars_uniq_xx_indep.append(freevars_uniq[i])

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
            elif basename == "f2_of_xx0_xx1":
                basefunc = f2_of_xx0_xx1
            elif basename == "f3_of_xx0":
                basefunc = f3_of_xx0
            elif basename == "f4_of_xx2":
                basefunc = f4_of_xx2
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

            frees_uniq = superfast_uniq(diff_result.free_symbols)
            xx_dep = False
            for dirn in range(3):
                if gri.xx[dirn] in frees_uniq:
                    xx_dep = True
            if xx_dep == False:
                freevars_uniq_xx_indep[i] = diff_result

        # Step 6.c: Finally, substitute integers for all functions & derivatives that evaluate to integers
        for varidx in range(len(freevars_uniq)):
            detgammahat = detgammahat.subs(freevars_uniq[varidx], freevars_uniq_xx_indep[varidx])
            for i in range(DIM):
                ReU[i] = ReU[i].subs(freevars_uniq[varidx], freevars_uniq_xx_indep[varidx])
                detgammahatdD[i] = detgammahatdD[i].subs(freevars_uniq[varidx], freevars_uniq_xx_indep[varidx])
                for j in range(DIM):
                    ReDD[i][j] = ReDD[i][j].subs(freevars_uniq[varidx], freevars_uniq_xx_indep[varidx])
                    ReUdD[i][j] = ReUdD[i][j].subs(freevars_uniq[varidx], freevars_uniq_xx_indep[varidx])
                    ghatDD[i][j] = ghatDD[i][j].subs(freevars_uniq[varidx], freevars_uniq_xx_indep[varidx])
                    ghatUU[i][j] = ghatUU[i][j].subs(freevars_uniq[varidx], freevars_uniq_xx_indep[varidx])
                    detgammahatdDD[i][j] = detgammahatdDD[i][j].subs(freevars_uniq[varidx],
                                                                     freevars_uniq_xx_indep[varidx])
                    for k in range(DIM):
                        ReDDdD[i][j][k] = ReDDdD[i][j][k].subs(freevars_uniq[varidx], freevars_uniq_xx_indep[varidx])
                        ReUdDD[i][j][k] = ReUdDD[i][j][k].subs(freevars_uniq[varidx], freevars_uniq_xx_indep[varidx])
                        ghatDDdD[i][j][k] = ghatDDdD[i][j][k].subs(freevars_uniq[varidx], freevars_uniq_xx_indep[varidx])
                        GammahatUDD[i][j][k] = GammahatUDD[i][j][k].subs(freevars_uniq[varidx],
                                                                         freevars_uniq_xx_indep[varidx])
                        for l in range(DIM):
                            ReDDdDD[i][j][k][l] = ReDDdDD[i][j][k][l].subs(freevars_uniq[varidx],
                                                                           freevars_uniq_xx_indep[varidx])
                            ghatDDdDD[i][j][k][l] = ghatDDdDD[i][j][k][l].subs(freevars_uniq[varidx],
                                                                               freevars_uniq_xx_indep[varidx])
                            GammahatUDDdD[i][j][k][l] = GammahatUDDdD[i][j][k][l].subs(freevars_uniq[varidx],
                                                                                       freevars_uniq_xx_indep[varidx])

        # Step 7: Construct needed C code for declaring rfmstruct, allocating storage for
        #         rfmstruct arrays, defining each element in each array, reading the
        #         rfmstruct data from memory (both with and without SIMD enabled), and
        #         freeing allocated memory for the rfmstrcut arrays.
        # struct_str: String that declares the rfmstruct struct.
        struct_str = "typedef struct __rfmstruct__ {\n"
        define_str = ""
        # rfmstruct stores pointers to (so far) 1D arrays. The malloc_str string allocates space for the arrays.
        malloc_str = "rfm_struct rfmstruct;\n"
        freemm_str = ""

        # readvr_str reads the arrays from memory as needed
        readvr_str = ["", "", ""]
        readvr_SIMD_outer_str = ["", "", ""]
        readvr_SIMD_inner_str = ["", "", ""]

        # Sort freevars_uniq_vals and freevars_uniq_xx_indep, according to alphabetized freevars_uniq_xx_indep.
        #    Without this step, the ordering of elements in rfmstruct would be random, and would change each time
        #    this function was called.
        if len(freevars_uniq_xx_indep) > 0:
            freevars_uniq_xx_indep, freevars_uniq_vals = (list(x) for x in zip(*sorted(zip(freevars_uniq_xx_indep, freevars_uniq_vals),key=str)))

        # Tease out how many variables each function in freevars_uniq_vals
        which_freevar = 0
        for expr in freevars_uniq_vals:
            if "_of_xx" in str(freevars_uniq_xx_indep[which_freevar]):
                frees = expr.free_symbols
                frees_uniq = superfast_uniq(frees)
                xx_list = []
                malloc_size = 1
                for i in range(3):
                    if gri.xx[i] in frees_uniq:
                        xx_list.append(gri.xx[i])
                        malloc_size *= gri.Nxx_plus_2NGHOSTS[i]

                struct_str += "\tREAL *restrict " + str(freevars_uniq_xx_indep[which_freevar]) + ";\n"
                malloc_str += "rfmstruct." + str(
                    freevars_uniq_xx_indep[which_freevar]) + " = (REAL *)malloc(sizeof(REAL)*" + str(malloc_size) + ");\n"
                freemm_str += "free(rfmstruct." + str(freevars_uniq_xx_indep[which_freevar]) + ");\n"
                output_define_and_readvr = False
                for dirn in range(3):
                    if (gri.xx[dirn] in frees_uniq) and not (gri.xx[(dirn+1)%3] in frees_uniq) and not (gri.xx[(dirn+2)%3] in frees_uniq):
                        define_str += "for(int i"+str(dirn)+"=0;i"+str(dirn)+"<Nxx_plus_2NGHOSTS"+str(dirn)+";i"+str(dirn)+"++) {\n"
                        define_str += "    const REAL xx"+str(dirn)+" = xx["+str(dirn)+"][i"+str(dirn)+"];\n"
                        define_str += "    rfmstruct." + str(freevars_uniq_xx_indep[which_freevar]) + "[i"+str(dirn)+"] = " + str(sp.ccode(freevars_uniq_vals[which_freevar])) + ";\n"
                        define_str += "}\n\n"
                        readvr_str[dirn] += "const REAL " + str(freevars_uniq_xx_indep[which_freevar]) + " = rfmstruct->" + \
                                         str(freevars_uniq_xx_indep[which_freevar]) + "[i"+str(dirn)+"];\n"
                        readvr_SIMD_outer_str[dirn] += "const double NOSIMD" + str(
                            freevars_uniq_xx_indep[which_freevar]) + " = rfmstruct->" + str(freevars_uniq_xx_indep[which_freevar]) + "[i"+str(dirn)+"]; "
                        readvr_SIMD_outer_str[dirn] += "const REAL_SIMD_ARRAY " + str(freevars_uniq_xx_indep[which_freevar]) + \
                                                    " = ConstSIMD(NOSIMD" + str(freevars_uniq_xx_indep[which_freevar]) + ");\n"
                        readvr_SIMD_inner_str[dirn] += "const REAL_SIMD_ARRAY " + str(freevars_uniq_xx_indep[which_freevar]) + \
                                                    " = ReadSIMD(&rfmstruct->" + str(freevars_uniq_xx_indep[which_freevar]) + "[i"+str(dirn)+"]);\n"
                        output_define_and_readvr = True

                if (output_define_and_readvr == False) and (gri.xx[0] in frees_uniq) and (gri.xx[1] in frees_uniq):
                    define_str += """
for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {
    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    rfmstruct.""" + str(freevars_uniq_xx_indep[which_freevar]) + """[i0 + Nxx_plus_2NGHOSTS0*i1] = """ + str(sp.ccode(freevars_uniq_vals[which_freevar])) + """;
}\n\n"""
                    readvr_str[0] += "const REAL " + str(freevars_uniq_xx_indep[which_freevar]) + " = rfmstruct->" + \
                                     str(freevars_uniq_xx_indep[which_freevar]) + "[i0 + Nxx_plus_2NGHOSTS0*i1];\n"
                    readvr_SIMD_outer_str[0] += "const double NOSIMD" + str(freevars_uniq_xx_indep[which_freevar]) + \
                                                " = rfmstruct->" + str(freevars_uniq_xx_indep[which_freevar]) + "[i0 + Nxx_plus_2NGHOSTS0*i1]; "
                    readvr_SIMD_outer_str[0] += "const REAL_SIMD_ARRAY " + str(freevars_uniq_xx_indep[which_freevar]) + \
                                                " = ConstSIMD(NOSIMD" + str(freevars_uniq_xx_indep[which_freevar]) + ");\n"
                    readvr_SIMD_inner_str[0] += "const REAL_SIMD_ARRAY " + str(freevars_uniq_xx_indep[which_freevar]) + \
                                                " = ReadSIMD(&rfmstruct->" + str(freevars_uniq_xx_indep[which_freevar]) + "[i0 + Nxx_plus_2NGHOSTS0*i1]);\n"
                    output_define_and_readvr = True

                if output_define_and_readvr == False:
                    print("ERROR: Could not figure out the (xx0,xx1,xx2) dependency within the expression for "+str(freevars_uniq_xx_indep[which_freevar])+":")
                    print(str(freevars_uniq_vals[which_freevar]))
                    sys.exit(1)

            which_freevar += 1

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

def get_EigenCoord():
    CoordSystem_orig = par.parval_from_str("reference_metric::CoordSystem")
    for EigenCoordstr in ["Spherical","Cylindrical","SymTP","Cartesian"]:
        if EigenCoordstr in CoordSystem_orig:
            return EigenCoordstr
    print("Error: Could not find EigenCoord for reference_metric::CoordSystem == "+CoordSystem_orig)
    sys.exit(1)

def set_Nxx_dxx_invdx_params__and__xx_h(outdir="."):
    import os
    with open(os.path.join(outdir,"set_Nxx_dxx_invdx_params__and__xx.h"),"w") as file:
        file.write("""
void set_Nxx_dxx_invdx_params__and__xx(const int EigenCoord, const int Nxx[3], 
                                       paramstruct *restrict params, REAL *restrict xx[3]) {
    // Override parameter defaults with values based on command line arguments and NGHOSTS.
    params->Nxx0 = Nxx[0];
    params->Nxx1 = Nxx[1];
    params->Nxx2 = Nxx[2];
    params->Nxx_plus_2NGHOSTS0 = Nxx[0] + 2*NGHOSTS;
    params->Nxx_plus_2NGHOSTS1 = Nxx[1] + 2*NGHOSTS;
    params->Nxx_plus_2NGHOSTS2 = Nxx[2] + 2*NGHOSTS;
    // Step 0d: Set up space and time coordinates
    // Step 0d.i: Declare \Delta x^i=dxx{0,1,2} and invdxx{0,1,2}, as well as xxmin[3] and xxmax[3]:
#include "set_Cparameters.h"
    REAL xxmin[3],xxmax[3];
    if(EigenCoord == 0) {
""")
        for i in range(3):
            file.write("        xxmin["+str(i)+"] = "+str(xxmin[i])+";\n")
            file.write("        xxmax["+str(i)+"] = "+str(xxmax[i])+";\n")
        file.write("""
    } else if (EigenCoord == 1) {
""")
        CoordSystem_orig = par.parval_from_str("reference_metric::CoordSystem")
        par.set_parval_from_str("reference_metric::CoordSystem",get_EigenCoord())
        reference_metric()
        for i in range(3):
            file.write("        xxmin["+str(i)+"] = "+str(xxmin[i])+";\n")
            file.write("        xxmax["+str(i)+"] = "+str(xxmax[i])+";\n")
        par.set_parval_from_str("reference_metric::CoordSystem",CoordSystem_orig)
        reference_metric()
        file.write("""
    }

    params->dxx0 = (xxmax[0] - xxmin[0]) / ((REAL)Nxx[0]);
    params->dxx1 = (xxmax[1] - xxmin[1]) / ((REAL)Nxx[1]);
    params->dxx2 = (xxmax[2] - xxmin[2]) / ((REAL)Nxx[2]);
    params->invdx0 = 1.0/params->dxx0;
    params->invdx1 = 1.0/params->dxx1;
    params->invdx2 = 1.0/params->dxx2;

    // Now that params.dxx{0,1,2} and params.invdxx{0,1,2} have been set,
    // Step 0d.iii: Set up uniform coordinate grids
    xx[0] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);
    for(int j=0;j<Nxx_plus_2NGHOSTS0;j++) 
        xx[0][j] = xxmin[0] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx0; // Cell-centered grid.
    xx[1] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS1);
    for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) 
        xx[1][j] = xxmin[1] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx1; // Cell-centered grid.
    xx[2] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS2);
    for(int j=0;j<Nxx_plus_2NGHOSTS2;j++) 
        xx[2][j] = xxmin[2] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx2; // Cell-centered grid.
    //fprintf(stderr,"hey inside setxx: %e %e %e | %e %e\\n",xxmin[0],xxmin[1],xxmin[2],xx[0][0],params->dxx0);
}
""")

def xxCart_h(funcname,cparamsloc,outfile):
    import outputC
    # Arbitrary-coordinate NRPy+ file output, Part 1: output the conversion from (x0,x1,x2) to Cartesian (x,y,z)
    Cout = outputC.outputC([xxCart[0],xxCart[1],xxCart[2]],
                           ["xCart[0]","xCart[1]","xCart[2]"],
                           "returnstring",params="preindent=1")

    with open(outfile, "w") as file:
        file.write("""
inline void """+funcname+"""(const paramstruct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]) {
#include """+"\""+cparamsloc+"\""+"""
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];\n"""+Cout+"}\n")

# Compute proper distance in all 3 directions. Used to find the appropriate timestep for the CFL condition.
def ds_dirn(delxx):
    ds_dirn = ixp.zerorank1(3)
    for i in range(3):
        ds_dirn[i] = delxx[i]*scalefactor_orthog[i]
    return ds_dirn

# Find the appropriate timestep for the CFL condition.
def add_find_timestep_func_to_dict():
    # Compute proper distance in all 3 directions.
    delxx = ixp.declarerank1("dxx", DIM=3)
    ds_drn = ds_dirn(delxx)

    ds_dirn_h = outputC([ds_drn[0], ds_drn[1], ds_drn[2]], ["ds_dirn0", "ds_dirn1", "ds_dirn2"],"returnstring")

    desc="Find the CFL-constrained timestep"
    add_to_Cfunction_dict(
        desc     =desc,
        type     ="REAL",
        name     ="find_timestep",
        params   ="const paramstruct *restrict params, REAL *restrict xx[3]",
        preloop  ="REAL dsmin = 1e38; // Start with a crazy high value... close to the largest number in single precision.",
        body     ="REAL ds_dirn0, ds_dirn1, ds_dirn2;\n"+ds_dirn_h+"""
#ifndef MIN
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
#endif
        // Set dsmin = MIN(dsmin, ds_dirn0, ds_dirn1, ds_dirn2);
        dsmin = MIN(dsmin,MIN(ds_dirn0,MIN(ds_dirn1,ds_dirn2)));
""",
        loopopts ="InteriorPoints,Read_xxs,DisableOpenMP",
        postloop ="return dsmin*CFL_FACTOR/wavespeed;\n")

def out_timestep_func_to_file(outfile):
    add_find_timestep_func_to_dict()
    with open(outfile, "w") as file:
        file.write(outC_function_dict["find_timestep"])

def out_default_free_parameters_for_rfm(free_parameters_file,
                                        domain_size=1.0,sinh_width=0.4,sinhv2_const_dr=0.05,SymTP_bScale=0.5):
    CoordSystem = par.parval_from_str("reference_metric::CoordSystem")

    with open(free_parameters_file, "a") as file:
        file.write("""
// Set free-parameter values.

const REAL domain_size    = """ + str(domain_size) + """;
const REAL sinh_width     = """ + str(sinh_width) + """;
const REAL sinhv2_const_dr= """ + str(sinhv2_const_dr) + """;
const REAL SymTP_bScale   = """ + str(SymTP_bScale) + ";\n")

        coordparams = ""
        if CoordSystem == "Spherical":
            coordparams += """
params.RMAX = domain_size;\n"""
        elif "SinhSpherical" in CoordSystem:
            coordparams += """
params.AMPL = domain_size;
params.SINHW=  sinh_width;\n"""
            if CoordSystem == "SinhSphericalv2":
                coordparams += "        params.const_dr = sinhv2_const_dr;\n"
        elif "SymTP" in CoordSystem:
            coordparams += """
params.bScale =  SymTP_bScale;
params.AMAX   =  domain_size;\n"""
            if CoordSystem == "SinhSymTP":
                coordparams += "        params.SINHWAA = sinh_width;\n"
        elif CoordSystem == "Cartesian":
            coordparams += """
params.xmin = -domain_size, params.xmax = domain_size;
params.ymin = -domain_size, params.ymax = domain_size;
params.zmin = -domain_size, params.zmax = domain_size;\n"""
        elif CoordSystem == "Cylindrical":
            coordparams += """
params.ZMIN   = -domain_size;
params.ZMAX   =  domain_size;
params.RHOMAX =  domain_size;\n"""
        elif "SinhCylindrical" in CoordSystem:
            coordparams += """
params.AMPLRHO = domain_size;
params.SINHWRHO= sinh_width;
params.AMPLZ   = domain_size;
params.SINHWZ  = sinh_width;\n"""
            if CoordSystem == "SinhCylindricalv2":
                coordparams += """
params.const_drho = sinhv2_const_dr;
params.const_dz   = sinhv2_const_dr;\n"""
        file.write(coordparams + "\n")
