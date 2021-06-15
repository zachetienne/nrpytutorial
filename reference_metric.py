# reference_metric.py: Define all needed quantities
#     for a reference metric.
# Given uniform (reference metric) coordinate
#    (xx[0],xx[1],xx[2]), you must define:
#     1) xxmin[3],xxmax[3]: Valid ranges for each
#       uniform coordinate xx0,xx1,xx2
#     2) xxSph[3]: Spherical coordinate (r,theta,phi),
#       in terms of uniform coordinate xx0,xx1,xx2
#     3) xx_to_Cart[3]: Cartesian coordinate (x,y,z),
#       in terms of uniform coordinate xx0,xx1,xx2
#     4) scalefactor_orthog:
#       orthogonal coordinate scale factor
#       (positive root of diagonal reference metric
#       components)
#     5) Cart_to_xx[3]: Inverse of xx_to_Cart:
#       xx0,xx1,xx2 as functions of (x,y,z).
#       In the case that there exists no closed-form
#       expression, then a root finder might be needed
#     6) UnitVectors[3][3]: Unit vectors of reference
#       metric.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

import sympy as sp                  # SymPy: The Python computer algebra package upon which NRPy+ depends
from outputC import outputC,superfast_uniq,add_to_Cfunction_dict # NRPy+: Core C code output module
# VVVVVVVVVVVVVVVVV
## TO BE DEPRECATED
from outputC import outC_function_dict
# ^^^^^^^^^^^^^^^^^
import NRPy_param_funcs as par      # NRPy+: Parameter interface
import grid as gri                  # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp            # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import os, sys                      # Standard Python modules for multiplatform OS-level functions

# Step 0a: Initialize parameters
thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "CoordSystem", "Spherical"))
par.initialize_param(par.glb_param("char", thismodule, "enable_rfm_precompute", "False"))
par.initialize_param(par.glb_param("char", thismodule, "rfm_precompute_Ccode_outdir", "Ccode"))

# Step 0b: Declare global variables
xx = gri.xx
xx_to_Cart = ixp.zerorank1(DIM=4)          # Must be set as a function of (xx[0],xx[1],xx[2])
Cart_to_xx = ixp.zerorank1(DIM=4)          # Must be set as a function of (xx[0],xx[1],xx[2])
xxSph      = ixp.zerorank1(DIM=4)          # Must be set as a function of (xx[0],xx[1],xx[2])
scalefactor_orthog = ixp.zerorank1(DIM=4)  # Must be set as a function of (xx[0],xx[1],xx[2])

Cartx, Carty, Cartz = sp.symbols("Cartx Carty Cartz", real=True)
Cart = [Cartx, Carty, Cartz]

scalefactor_orthog_funcform = ixp.zerorank1(DIM=4) # Must be set in terms of generic functions of xx[]s

# The following are necessary since SymPy has trouble with its native sinh and cosh functions.
def nrpysinh(x):
    return (sp.exp(x) - sp.exp(-x)) * sp.Rational(1, 2)
def nrpycosh(x):
    return (sp.exp(x) + sp.exp(-x)) * sp.Rational(1, 2)

have_already_called_reference_metric_function = False

def reference_metric(SymPySimplifyExpressions=True, enable_compute_hatted_quantities=True):
    global f0_of_xx0_funcform, f1_of_xx1_funcform, f2_of_xx0_xx1_funcform, f3_of_xx0_funcform, f4_of_xx2_funcform
    global f0_of_xx0, f1_of_xx1, f2_of_xx1, f2_of_xx0_xx1, f3_of_xx0, f4_of_xx2
    f0_of_xx0_funcform     = sp.Function('f0_of_xx0_funcform')(xx[0])
    f1_of_xx1_funcform     = sp.Function('f1_of_xx1_funcform')(xx[1])
    f2_of_xx1_funcform     = sp.Function('f2_of_xx1_funcform')(xx[1])
    f2_of_xx0_xx1_funcform = sp.Function('f2_of_xx0_xx1_funcform')(xx[0], xx[1])
    f3_of_xx0_funcform     = sp.Function('f3_of_xx0_funcform')(xx[0])
    f4_of_xx2_funcform     = sp.Function('f4_of_xx2_funcform')(xx[2])
    f0_of_xx0, f1_of_xx1, f2_of_xx1, f2_of_xx0_xx1, f3_of_xx0, f4_of_xx2 = par.Cparameters("REAL", thismodule,
                                  ["f0_of_xx0", "f1_of_xx1", "f2_of_xx1", "f2_of_xx0_xx1", "f3_of_xx0", "f4_of_xx2"], 1e300)
    # FIXME: Hack
    # return values of par.Cparameters() in the following code block are unused, so we ignore them.
    par.Cparameters("REAL", thismodule, ["f0_of_xx0__D0", "f0_of_xx0__DD00","f0_of_xx0__DDD000"], 1e300)
    par.Cparameters("REAL", thismodule, ["f1_of_xx1__D1", "f1_of_xx1__DD11","f1_of_xx1__DDD111"], 1e300)
    par.Cparameters("REAL", thismodule, ["f2_of_xx1__D1", "f2_of_xx1__DD11","f2_of_xx1__DDD111"], 1e300)
    par.Cparameters("REAL", thismodule,
                    ["f2_of_xx0_xx1__D0", "f2_of_xx0_xx1__D1", "f2_of_xx0_xx1__DD00", "f2_of_xx0_xx1__DD11"], 1e300)
    par.Cparameters("REAL", thismodule, ["f3_of_xx0__D0", "f3_of_xx0__DD00"], 1e300)
    par.Cparameters("REAL", thismodule, ["f4_of_xx2__D2", "f4_of_xx2__DD22"], 1e300)

    global have_already_called_reference_metric_function # setting to global enables other modules to see updated value.
    have_already_called_reference_metric_function = True

    CoordSystem = par.parval_from_str("reference_metric::CoordSystem")
    M_PI, M_SQRT1_2 = par.Cparameters("#define", thismodule, ["M_PI", "M_SQRT1_2"], "")

    global xxmin
    global xxmax

    global UnitVectors
    UnitVectors = ixp.zerorank2(DIM=3)

    # Set up hatted metric tensor, rescaling matrix, and rescaling vector

    #####################################################################
    # SPHERICAL-LIKE COORDINATE SYSTEMS WITH & WITHOUT RADIAL RESCALING #
    #####################################################################
    if CoordSystem in ('Spherical', 'SinhSpherical', 'SinhSphericalv2'):

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
        xx_to_Cart[0] = xxSph[0]*sp.sin(xxSph[1])*sp.cos(xxSph[2])
        xx_to_Cart[1] = xxSph[0]*sp.sin(xxSph[1])*sp.sin(xxSph[2])
        xx_to_Cart[2] = xxSph[0]*sp.cos(xxSph[1])

        scalefactor_orthog[0] = sp.diff(xxSph[0],xx[0])
        scalefactor_orthog[1] = xxSph[0]
        scalefactor_orthog[2] = xxSph[0]*sp.sin(xxSph[1])

        f0_of_xx0 = xxSph[0]
        f1_of_xx1 = sp.sin(xxSph[1])
        scalefactor_orthog_funcform[0] = sp.diff(f0_of_xx0_funcform,xx[0])
        scalefactor_orthog_funcform[1] = f0_of_xx0_funcform
        scalefactor_orthog_funcform[2] = f0_of_xx0_funcform*f1_of_xx1_funcform

        # Set the unit vectors
        UnitVectors = [[sp.sin(xxSph[1])*sp.cos(xxSph[2]), sp.sin(xxSph[1])*sp.sin(xxSph[2]),  sp.cos(xxSph[1])],
                       [sp.cos(xxSph[1])*sp.cos(xxSph[2]), sp.cos(xxSph[1])*sp.sin(xxSph[2]), -sp.sin(xxSph[1])],
                       [                -sp.sin(xxSph[2]),                  sp.cos(xxSph[2]),  sp.sympify(0)   ]]

    ######################################################################
    # SPHERICAL-LIKE COORDINATE SYSTEMS WITH RADIAL AND THETA RESCALINGS #
    ######################################################################
    elif CoordSystem in ('NobleSphericalThetaOptionOne', 'NobleSphericalThetaOptionTwo'):
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
        xx_to_Cart[0] = xxSph[0]*sp.sin(xxSph[1])*sp.cos(xxSph[2])
        xx_to_Cart[1] = xxSph[0]*sp.sin(xxSph[1])*sp.sin(xxSph[2])
        xx_to_Cart[2] = xxSph[0]*sp.cos(xxSph[1])

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
    elif CoordSystem in ('Cylindrical', 'SinhCylindrical', 'SinhCylindricalv2'):
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

        xx_to_Cart[0] = RHOCYL*sp.cos(PHICYL)
        xx_to_Cart[1] = RHOCYL*sp.sin(PHICYL)
        xx_to_Cart[2] = ZCYL

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

    elif CoordSystem in ('SymTP', 'SinhSymTP'):
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
            xxmax[0] = sp.sympify(1)
            # With xxmax[0] = 1, sinh(xx0/SINHWAA) / sinh(1/SINHWAA) will evaluate to a number between 0 and 1.
            #   Then AA = AMAX * sinh(xx0/SINHWAA) / sinh(1/SINHWAA) will evaluate to a number between 0 and AMAX.
            AA = AMAX * (sp.exp(xx[0] / SINHWAA) - sp.exp(-xx[0] / SINHWAA)) / (sp.exp(1 / SINHWAA) - sp.exp(-1 / SINHWAA))

        var1 = sp.sqrt(AA**2 + (bScale * sp.sin(xx[1]))**2)
        var2 = sp.sqrt(AA**2 + bScale**2)

        RHOSYMTP = AA*sp.sin(xx[1])
        PHSYMTP = xx[2]
        ZSYMTP = var2*sp.cos(xx[1])

        xx_to_Cart[0] = AA  *sp.sin(xx[1])*sp.cos(xx[2])
        xx_to_Cart[1] = AA  *sp.sin(xx[1])*sp.sin(xx[2])
        xx_to_Cart[2] = ZSYMTP

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

    #####################################
    # CARTESIAN-LIKE COORDINATE SYSTEMS #
    #####################################
    elif CoordSystem == "Cartesian":
        # return values of par.Cparameters() in the following line of code are unused, so we ignore them.
        par.Cparameters("REAL",thismodule, ["xmin","xmax","ymin","ymax","zmin","zmax"],
                                           [ -10.0,  10.0, -10.0,  10.0, -10.0,  10.0])
        xxmin = ["xmin", "ymin", "zmin"]
        xxmax = ["xmax", "ymax", "zmax"]

        xx_to_Cart[0] = xx[0]
        xx_to_Cart[1] = xx[1]
        xx_to_Cart[2] = xx[2]

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

    elif CoordSystem == "SinhCartesian":
        # SinhCartesian coordinates allows us to push the outer boundary of the
        # computational domain a lot further away, while keeping reasonably high
        # resolution towards the center of the computational grid.

        # Set default values for min and max (x,y,z)
        xxmin = [sp.sympify(-1), sp.sympify(-1), sp.sympify(-1)]
        xxmax = [sp.sympify(+1), sp.sympify(+1), sp.sympify(+1)]

        # Declare basic parameters of the coordinate system and their default values
        AMPLXYZ, SINHWXYZ = par.Cparameters("REAL", thismodule,
                                            ["AMPLXYZ", "SINHWXYZ"],
                                            [     10.0,        0.2])

        # Compute (xx_to_Cart0,xx_to_Cart1,xx_to_Cart2) from (xx0,xx1,xx2)
        for ii in [0, 1, 2]:
            xx_to_Cart[ii] = AMPLXYZ*(sp.exp(xx[ii]/SINHWXYZ) - sp.exp(-xx[ii]/SINHWXYZ))/(sp.exp(1/SINHWXYZ) - sp.exp(-1/SINHWXYZ))

        # Compute (r,th,ph) from (xx_to_Cart2,xx_to_Cart1,xx_to_Cart2)
        xxSph[0] = sp.sqrt(xx_to_Cart[0] ** 2 + xx_to_Cart[1] ** 2 + xx_to_Cart[2] ** 2)
        xxSph[1] = sp.acos(xx_to_Cart[2] / xxSph[0])
        xxSph[2] = sp.atan2(xx_to_Cart[1], xx_to_Cart[0])

        # Compute (xx0,xx1,xx2) from (Cartx,Carty,Cartz)
        Cart_to_xx[0] = SINHWXYZ*sp.asinh(Cartx*sp.sinh(1/SINHWXYZ)/AMPLXYZ)
        Cart_to_xx[1] = SINHWXYZ*sp.asinh(Carty*sp.sinh(1/SINHWXYZ)/AMPLXYZ)
        Cart_to_xx[2] = SINHWXYZ*sp.asinh(Cartz*sp.sinh(1/SINHWXYZ)/AMPLXYZ)

        # Compute scale factors
        scalefactor_orthog[0] = sp.diff(xx_to_Cart[0], xx[0])
        scalefactor_orthog[1] = sp.diff(xx_to_Cart[1], xx[1])
        scalefactor_orthog[2] = sp.diff(xx_to_Cart[2], xx[2])

        f0_of_xx0             = sp.diff(xx_to_Cart[0], xx[0])
        f1_of_xx1             = sp.diff(xx_to_Cart[1], xx[1])
        f4_of_xx2             = sp.diff(xx_to_Cart[2], xx[2])

        scalefactor_orthog_funcform[0] = f0_of_xx0_funcform
        scalefactor_orthog_funcform[1] = f1_of_xx1_funcform
        scalefactor_orthog_funcform[2] = f4_of_xx2_funcform

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

    # enable_rfm_precompute: precompute and store in memory possibly
    #     complex expressions related to the reference metric (a.k.a.,
    #      "hatted quantities")

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
        sympy_version_decimal = float(int(sympy_version.split(".")[0]) + int(sympy_version.split(".")[1])/10.0)
        is_old_sympy_version = sympy_version_decimal < 1.2
        # The derivative representation changed with SymPy 1.2, forcing version-dependent behavior.

        # Example: Derivative(f0_of_xx0_funcform(xx0)(xx0), (xx0, 2)) >> f0_of_xx0__DD00
        rule = {} # replacement dictionary
        for item in sp.preorder_traversal(expr):
            if item.func == sp.Derivative:
                # extract function name before '_funcform'
                strfunc = str(item.args[0]).split('_funcform(', 1)[0]
                if is_old_sympy_version:
                    # extract differentiation variable and derivative order (SymPy <= 1.1)
                    var, order = str(item.args[1])[2:], len(item.args) - 1
                else:
                    # extract differentiation variable and derivative order (SymPy >= 1.2)
                    var, order = str(item.args[1][0])[2:], item.args[1][1]
                # build derivative operator with format: __DD...D(var)(var)...(var) where
                # D and (var) are repeated for every derivative order
                oper = '__D' + 'D'*(order - 1) + var*order
                # add replacement rule to dictionary
                rule[item] = sp.sympify(strfunc + oper)
        expr = expr.xreplace(rule); rule = {}

        # Example: f0_of_xx0_funcform(xx0)(xx0) >> f0_of_xx0
        for item in sp.preorder_traversal(expr):
            if "_funcform" in str(item.func):
                # extract function name before '_funcform'
                strfunc = str(item.func).split("_funcform", 1)[0]
                # add replacement rule to dictionary
                rule[item] = sp.sympify(strfunc)
        return expr.xreplace(rule)

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
    for i, var in enumerate(freevars_uniq):
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
        has_xx_dependence = False
        for dirn in range(3):
            if gri.xx[dirn] in frees_uniq:
                has_xx_dependence = True
        if not has_xx_dependence:
            freevars_uniq_xx_indep[i] = diff_result

    # Step 6.c: Finally, substitute integers for all functions & derivatives that evaluate to integers
    for varidx, freevar in enumerate(freevars_uniq):
        detgammahat = detgammahat.subs(freevar, freevars_uniq_xx_indep[varidx])
        for i in range(DIM):
            ReU[i] = ReU[i].subs(freevar, freevars_uniq_xx_indep[varidx])
            detgammahatdD[i] = detgammahatdD[i].subs(freevar, freevars_uniq_xx_indep[varidx])
            for j in range(DIM):
                ReDD[i][j] = ReDD[i][j].subs(freevar, freevars_uniq_xx_indep[varidx])
                ReUdD[i][j] = ReUdD[i][j].subs(freevar, freevars_uniq_xx_indep[varidx])
                ghatDD[i][j] = ghatDD[i][j].subs(freevar, freevars_uniq_xx_indep[varidx])
                ghatUU[i][j] = ghatUU[i][j].subs(freevar, freevars_uniq_xx_indep[varidx])
                detgammahatdDD[i][j] = detgammahatdDD[i][j].subs(freevar,
                                                                 freevars_uniq_xx_indep[varidx])
                for k in range(DIM):
                    ReDDdD[i][j][k] = ReDDdD[i][j][k].subs(freevar, freevars_uniq_xx_indep[varidx])
                    ReUdDD[i][j][k] = ReUdDD[i][j][k].subs(freevar, freevars_uniq_xx_indep[varidx])
                    ghatDDdD[i][j][k] = ghatDDdD[i][j][k].subs(freevar, freevars_uniq_xx_indep[varidx])
                    GammahatUDD[i][j][k] = GammahatUDD[i][j][k].subs(freevar,
                                                                     freevars_uniq_xx_indep[varidx])
                    for l in range(DIM):
                        ReDDdDD[i][j][k][l] = ReDDdDD[i][j][k][l].subs(freevar,
                                                                       freevars_uniq_xx_indep[varidx])
                        ghatDDdDD[i][j][k][l] = ghatDDdDD[i][j][k][l].subs(freevar,
                                                                           freevars_uniq_xx_indep[varidx])
                        GammahatUDDdD[i][j][k][l] = GammahatUDDdD[i][j][k][l].subs(freevar,
                                                                                   freevars_uniq_xx_indep[varidx])

    # Step 7: Construct needed C code for declaring rfmstruct, allocating storage for
    #         rfmstruct arrays, defining each element in each array, reading the
    #         rfmstruct data from memory (both with and without SIMD enabled), and
    #         freeing allocated memory for the rfmstruct arrays.
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
    with open(os.path.join(outdir, "rfm_struct__declare.h"), "w") as file:
        file.write(struct_str)
    with open(os.path.join(outdir, "rfm_struct__malloc.h"), "w") as file:
        file.write(malloc_str)
    with open(os.path.join(outdir, "rfm_struct__define.h"), "w") as file:
        file.write(define_str)
    for i in range(3):
        with open(os.path.join(outdir, "rfm_struct__read" + str(i) + ".h"), "w") as file:
            file.write(readvr_str[i])
        with open(os.path.join(outdir, "rfm_struct__SIMD_outer_read" + str(i) + ".h"), "w") as file:
            file.write(readvr_SIMD_outer_str[i])
        with open(os.path.join(outdir, "rfm_struct__SIMD_inner_read" + str(i) + ".h"), "w") as file:
            file.write(readvr_SIMD_inner_str[i])
    with open(os.path.join(outdir, "rfm_struct__freemem.h"), "w") as file:
        file.write(freemm_str)

####################################################
# Core Jacobian (basis) transformation functions,
#      for reference metric basis to/from the
#      Cartesian basis.

# We define Jacobians relative to the reference metric
#   basis at a point x^j_rfm=(xx0,xx1,xx2)_rfm on the source grid:
#
#  Jac_dUCart_dDrfmUD[i][j] = dx^i_Cart / dx^j_rfm
#
# via exact differentiation (courtesy SymPy), and the inverse Jacobian
#
#  Jac_dUrfm_dDCartUD[i][j] = dx^i_rfm / dx^j_Cart
#
# using NRPy+'s generic_matrix_inverter3x3() function

def compute_Jacobian_and_inverseJacobian_tofrom_Cartesian():
    # Step 2.a: First construct Jacobian matrix:
    Jac_dUCart_dDrfmUD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            Jac_dUCart_dDrfmUD[i][j] = sp.diff(xx_to_Cart[i], xx[j])
    Jac_dUrfm_dDCartUD, dummyDET = ixp.generic_matrix_inverter3x3(Jac_dUCart_dDrfmUD)
    return Jac_dUCart_dDrfmUD, Jac_dUrfm_dDCartUD

def basis_transform_vectorU_from_rfmbasis_to_Cartesian(Jac_dUCart_dDrfmUD, src_vectorU):
    Cart_dst_vectorU = ixp.zerorank1()
    for i in range(3):
        for l in range(3):
            Cart_dst_vectorU[i] += Jac_dUCart_dDrfmUD[i][l] * src_vectorU[l]
    return Cart_dst_vectorU

def basis_transform_tensorDD_from_rfmbasis_to_Cartesian(Jac_dUrfm_dDCartUD, src_tensorDD):
    Cart_dst_tensorDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            for l in range(3):
                for m in range(3):
                    Cart_dst_tensorDD[i][j] += Jac_dUrfm_dDCartUD[l][i]*Jac_dUrfm_dDCartUD[m][j]*src_tensorDD[l][m]
    return Cart_dst_tensorDD

def basis_transform_vectorU_from_Cartesian_to_rfmbasis(Jac_dUrfm_dDCartUD, Cart_src_vectorU):
    rfm_dst_vectorU = ixp.zerorank1()
    for i in range(3):
        for l in range(3):
            rfm_dst_vectorU[i] += Jac_dUrfm_dDCartUD[i][l] * Cart_src_vectorU[l]
    return rfm_dst_vectorU

def basis_transform_tensorDD_from_Cartesian_to_rfmbasis(Jac_dUCart_dDrfmUD, Cart_src_tensorDD):
    rfm_dst_tensorDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            for l in range(3):
                for m in range(3):
                    rfm_dst_tensorDD[i][j] += Jac_dUCart_dDrfmUD[l][i]*Jac_dUCart_dDrfmUD[m][j]*Cart_src_tensorDD[l][m]
    return rfm_dst_tensorDD
##################################################

def get_EigenCoord():
    CoordSystem_orig = par.parval_from_str("reference_metric::CoordSystem")
    for EigenCoordstr in ["Spherical", "Cylindrical", "SymTP", "Cartesian"]:
        if EigenCoordstr in CoordSystem_orig:
            return EigenCoordstr
    print("Error: Could not find EigenCoord for reference_metric::CoordSystem == "+CoordSystem_orig)
    sys.exit(1)

def add_to_Cfunc_dict__set_Nxx_dxx_invdx_params__and__xx(rel_path_to_Cparams=os.path.join("./"), NGHOSTS_is_a_param=False):
    gridsuffix = ""  # Disable for now.

    def set_xxmin_xxmax():
        outstr = ""
        for dirn in range(3):
            outstr += "        xxmin[" + str(dirn) + "] = " + str(xxmin[dirn]) + ";\n"
            outstr += "        xxmax[" + str(dirn) + "] = " + str(xxmax[dirn]) + ";\n"
        return outstr
    body = """
    // Override parameter defaults with values based on command line arguments and NGHOSTS.
    params->Nxx0""" + gridsuffix + r""" = Nxx[0];
    params->Nxx1""" + gridsuffix + r""" = Nxx[1];
    params->Nxx2""" + gridsuffix + r""" = Nxx[2];
"""
    NGHOSTS_prefix=""
    if NGHOSTS_is_a_param:
        NGHOSTS_prefix="params->"
    body += """
    params->Nxx_plus_2NGHOSTS0""" + gridsuffix + """ = Nxx[0] + 2*"""+NGHOSTS_prefix+"""NGHOSTS;
    params->Nxx_plus_2NGHOSTS1""" + gridsuffix + """ = Nxx[1] + 2*"""+NGHOSTS_prefix+"""NGHOSTS;
    params->Nxx_plus_2NGHOSTS2""" + gridsuffix + """ = Nxx[2] + 2*"""+NGHOSTS_prefix+"""NGHOSTS;
    // Now that params->Nxx_plus_2NGHOSTS* has been set, and below we need e.g., Nxx_plus_2NGHOSTS*, we include set_Cparameters.h here:
#include \"""" + os.path.join(rel_path_to_Cparams, "set_Cparameters.h") + """\"
    // Step 0d: Set up space and time coordinates
    // Step 0d.i: Declare Delta x^i=dxx{0,1,2} and invdxx{0,1,2}, as well as xxmin[3] and xxmax[3]:
    REAL xxmin[3],xxmax[3];
    if(EigenCoord == 0) {
"""
    body += set_xxmin_xxmax() + """    } else if (EigenCoord == 1) {
"""
    CoordSystem_orig = par.parval_from_str("reference_metric::CoordSystem")
    # If we are using a "holey" Spherical-like coordinate, for certain grids xx0min>0 is
    #    such that xx[0][0] is negative, which causes "Cartesian disagreement" errors.
    if "Spherical" not in CoordSystem_orig:
        par.set_parval_from_str("reference_metric::CoordSystem", get_EigenCoord())
        reference_metric()
        body += set_xxmin_xxmax()
        par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem_orig)
        reference_metric()
    else:
        body += set_xxmin_xxmax()

    # Now set grid spacing dxx, invdx = 1/dxx, and xx[]
    body += """    }
    // Step 0d.iii: Set params.dxx{0,1,2}, params.invdx{0,1,2}, and uniform coordinate grids xx[3][]
"""
    for dirn in ["0", "1", "2"]:
        body += "    params->dxx"+dirn+gridsuffix+" = (xxmax["+dirn+"] - xxmin["+dirn+"]) / ((REAL)Nxx["+dirn+"]);\n"
        body += "    params->invdx"+dirn+gridsuffix+" = 1.0/params->dxx"+dirn+gridsuffix+";\n"
        body += """    xx["""+dirn+"""] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS"""+dirn+gridsuffix + """);
    for(int j=0;j<Nxx_plus_2NGHOSTS"""+dirn+gridsuffix+""";j++)
        xx["""+dirn+"""][j] = xxmin["""+dirn+"""] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx"""+dirn+gridsuffix+"""; // Cell-centered grid.\n"""
        if dirn != "2":
            body += "\n"

    add_to_Cfunction_dict(
        includes=["stdio.h", "math.h", "stdlib.h",
                  os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h"),
                  os.path.join(rel_path_to_Cparams, "declare_Cparameters_struct.h")],
        desc  ="Override default values for Nxx{0,1,2}, Nxx_plus_2NGHOSTS{0,1,2}, dxx{0,1,2}, and invdx{0,1,2}; and set xx[3][]",
        type  ="void",
        name  ="set_Nxx_dxx_invdx_params__and__xx"+gridsuffix,
        params="const int EigenCoord, const int Nxx[3],paramstruct *restrict params, REAL *restrict xx[3]",
        body  =body,
        opts  ="DisableCparameters")  # Cparameters here must be #include'd in body, not at top of function as usual.

def add_to_Cfunc_dict__xx_to_Cart(rel_path_to_Cparams=os.path.join("./")):
    gridsuffix = ""  # Disable for now
    # Arbitrary-coordinate NRPy+ file output, Part 1: output the conversion from (x0,x1,x2) to Cartesian (x,y,z)
    #    Suppose grid origin is at 1,1,1. Then the Cartesian gridpoint at 1,2,3 will be 2,3,4; hence
    #    the xx_to_Cart[i]+gri.Cart_origin[i] below:
    body = """
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
    """ + outputC([xx_to_Cart[0]+gri.Cart_origin[0],
                   xx_to_Cart[1]+gri.Cart_origin[1],
                   xx_to_Cart[2]+gri.Cart_origin[2]],
                  ["xCart[0]", "xCart[1]", "xCart[2]"],
                  "returnstring", params="preindent=1"). \
        replace("Cart_originx", "Cart_originx" + gridsuffix).\
        replace("Cart_originy", "Cart_originy" + gridsuffix).\
        replace("Cart_originz", "Cart_originz" + gridsuffix)

    add_to_Cfunction_dict(
        includes=["stdio.h", "math.h",
                  os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h"),
                  os.path.join(rel_path_to_Cparams, "declare_Cparameters_struct.h")],
        desc    ="Compute Cartesian coordinates given local grid coordinate (xx0,xx1,xx2), "
                 "  accounting for the origin of this grid being possibly offcenter.",
        type    ="void",
        name    ="xx_to_Cart"+gridsuffix,
        params  ="const paramstruct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]",
        body    =body,
        rel_path_to_Cparams=rel_path_to_Cparams)


# Compute proper distance in all 3 directions. Used to find the appropriate timestep for the CFL condition.
def ds_dirn(delxx, append_gridsuffix_to_xx=False):
    gridsuffix = ""  # Disable for now
    scalefactor_orthog_inj = []
    for i in range(3):
        if append_gridsuffix_to_xx:
            scalefactor_orthog_inj.append(scalefactor_orthog[i].
                                          subs(xx[0], sp.sympify(str(xx[0]) + gridsuffix)).
                                          subs(xx[1], sp.sympify(str(xx[1]) + gridsuffix)).
                                          subs(xx[2], sp.sympify(str(xx[2]) + gridsuffix)))
        else:
            scalefactor_orthog_inj.append(scalefactor_orthog[i])

    ds_dirn = ixp.zerorank1(3)
    for i in range(3):
        ds_dirn[i] = delxx[i]*scalefactor_orthog_inj[i]
    return ds_dirn


# Find the appropriate timestep for the CFL condition.
def add_to_Cfunc_dict__find_timestep(rel_path_to_Cparams=os.path.join("./"), enable_mask=False,
                                     output_dt_local_h_only=False):
    gridsuffix = ""  # Disable for now
    ##############################
    # Step 1: Function description
    desc = "Find the CFL-constrained timestep"
    ##############################
    # Step 2: Function return type
    type = "REAL"
    ##############################
    # Step 3: Function name
    name = "find_timestep" + gridsuffix
    ##############################
    # Step 4: Prior to the main loop
    preloop = "    REAL dsmin = 1e38; // Start with a crazy high value... close to the largest number in single precision."
    ##############################
    # Step 5: Loop options
    loopopts = "AllPoints,Read_xxs,DisableOpenMP"
    if gridsuffix != "":
        loopopts += ","+gridsuffix
    ##############################
    # Step 6: function input parameters
    params = "const paramstruct *restrict params, REAL *restrict xx[3]"
    if enable_mask:
        params += ", int8_t *restrict mask"
    ##############################
    # Step 7: function body
    # Compute proper distance in all 3 directions.
    if output_dt_local_h_only:
        ds_drn = ds_dirn(gri.dxx, append_gridsuffix_to_xx=True)
    else:
        ds_drn = ds_dirn(gri.dxx, append_gridsuffix_to_xx=False)
    ds_dirn_h = outputC([ds_drn[0], ds_drn[1], ds_drn[2]], ["ds_dirn0", "ds_dirn1", "ds_dirn2"], "returnstring").\
        replace("dxx0", "dxx0"+gridsuffix).\
        replace("dxx1", "dxx1"+gridsuffix).\
        replace("dxx2", "dxx2"+gridsuffix)
    body = "REAL ds_dirn0, ds_dirn1, ds_dirn2;\n" + ds_dirn_h + """
#ifndef MIN
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
#endif\n"""
    indent = ""
    if enable_mask:
        body += "if(mask[IDX3S" + gridsuffix + "(i0,i1,i2)] >= 0) {\n"
        indent = "    "
    if not output_dt_local_h_only:
        # not output_dt_local_h_only -> seeking dsmin over the entire grid, over all directions
        body += indent + "// Set dsmin = MIN(dsmin, ds_dirn0, ds_dirn1, ds_dirn2):\n"
        body += indent + "dsmin = MIN(dsmin, MIN(ds_dirn0, MIN(ds_dirn1, ds_dirn2)));\n"
    else:
        # output_dt_local_h_only means we seek a minimum over all directions at given gridpoint only
        body += indent + "// Set dt_local["+gridsuffix.replace("_grid", "")+"] = MIN(ds_dirn0, ds_dirn1, ds_dirn2) * CFL_FACTOR/wavespeed :\n"
        body += indent + "dt_local["+gridsuffix.replace("_grid", "")+"] = MIN(ds_dirn0, MIN(ds_dirn1, ds_dirn2)) * CFL_FACTOR/wavespeed;\n"
    if enable_mask:
        body += "}\n"

    if output_dt_local_h_only:
        return body.\
            replace("ds_dirn0", "ds_dirn0"+gridsuffix).\
            replace("ds_dirn1", "ds_dirn1"+gridsuffix).\
            replace("ds_dirn2", "ds_dirn2"+gridsuffix)

    ##############################
    # Step 8: after the loop
    postloop = "    return dsmin*CFL_FACTOR/wavespeed;\n"
    ##############################
    # Step 9: add to Cfunction dictionary
    add_to_Cfunction_dict(
        includes=["stdio.h", "math.h",
                  os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h"),
                  os.path.join(rel_path_to_Cparams, "declare_Cparameters_struct.h")],
        desc    =desc,
        type    =type,
        name    =name,
        params  =params,
        preloop =preloop,
        body    =body,
        loopopts=loopopts,
        postloop=postloop,
        rel_path_to_Cparams=rel_path_to_Cparams)


# Find the appropriate timestep for the CFL condition.
def add_to_Cfunc_dict__find_dsmin(rel_path_to_Cparams=os.path.join("./")):
    gridsuffix = ""  # Disable for now
    desc = "Find dsmin = min_i sqrt(ghat_{ii} dx^i dx^i)"
    type = "REAL"
    name = "find_dsmin" + gridsuffix
    params = "const paramstruct *restrict params, const int i0i1i2[3], const REAL *restrict xx[3]"
    body = """
  REAL dsmin = 1e38; // Start with a crazy high value... close to the largest number in single precision.
  const REAL xx0 = xx[0][i0i1i2[0]];
  const REAL xx1 = xx[1][i0i1i2[1]];
  const REAL xx2 = xx[2][i0i1i2[2]];
"""
    # Compute proper distance in all 3 directions.
    ds_drn = ds_dirn(gri.dxx, append_gridsuffix_to_xx=False)
    body += "  REAL ds_dirn0, ds_dirn1, ds_dirn2;\n"
    body += outputC([ds_drn[0], ds_drn[1], ds_drn[2]], ["ds_dirn0", "ds_dirn1", "ds_dirn2"], "returnstring",
                        params="outCverbose=false,includebraces=false").\
                    replace("dxx0", "dxx0"+gridsuffix).\
                    replace("dxx1", "dxx1"+gridsuffix).\
                    replace("dxx2", "dxx2"+gridsuffix)
    body += """
#ifndef MIN
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
#endif\n"""
    indent = ""
    body += indent + "  // Set dsmin = MIN(dsmin, ds_dirn0, ds_dirn1, ds_dirn2):\n"
    body += indent + "  return MIN(dsmin, MIN(ds_dirn0, MIN(ds_dirn1, ds_dirn2)));\n"
    add_to_Cfunction_dict(
        includes=["stdio.h", "math.h",
                  os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h"),
                  os.path.join(rel_path_to_Cparams, "declare_Cparameters_struct.h")],
        desc    =desc,
        type    =type,
        name    =name,
        params  =params,
        body    =body,
        rel_path_to_Cparams=rel_path_to_Cparams)




# Step 8.a.iv: Generate Cart_to_xx.h, which contains Cart_to_xx_grid*()
#              for mapping from Cartesian->xx for the chosen CoordSystem.
def add_to_Cfunc_dict__Cart_to_xx_and_nearest_i0i1i2(rel_path_to_Cparams=os.path.join("./"), relative_to="local_grid_center"):
    gridsuffix = ""  # Disable for now
    CoordSystem = par.parval_from_str("reference_metric::CoordSystem")

    prefunc = ""
    desc = """Given Cartesian point (x,y,z), this function outputs the corresponding
  (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid"""
    namesuffix = ""
    if relative_to == "global_grid_center":
        namesuffix = "_" + relative_to
    name = "Cart_to_xx_and_nearest_i0i1i2" + namesuffix + gridsuffix
    params = "const paramstruct *restrict params, const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]"

    preloop = ""
    if relative_to == "local_grid_center":
        preloop = """
    // First compute the closest (xx0,xx1,xx2) to the given Cartesian gridpoint (x,y,z),
    //   *relative* to the center of the local grid.
    //   So for example,
    //   1) if global xCart[012] = (1,1,1), and the
    //      origin of the grid is at global xCart (x,y,z) = (1,1,1), then
    //      (Cartx,Carty,Cartz) = (0,0,0)
    //   2) if global xCart[012] = (0,0,0), and the
    //      origin of the grid is at global xCart (x,y,z) = (1,1,1), then
    //      (Cartx,Carty,Cartz) = (-1,-1,-1)
    // Therefore, (Cartx,Carty,Cartz) = (xCart[0]-originx, xCart[1]-originy, xCart[2]-originz)
    const REAL Cartx = xCart[0] - Cart_originx_GRIDSFX_;
    const REAL Carty = xCart[1] - Cart_originy_GRIDSFX_;
    const REAL Cartz = xCart[2] - Cart_originz_GRIDSFX_;
""".replace("_GRIDSFX_", gridsuffix)
    elif relative_to == "global_grid_center":
        preloop = """
    const REAL Cartx = xCart[0];
    const REAL Carty = xCart[1];
    const REAL Cartz = xCart[2];
"""
    else:
        print("Error: relative_to must be set to either local_grid_center or global_grid_center. " + relative_to + " was chosen.")
        sys.exit(1)

    if "theta_adj" in CoordSystem:
        body = outputC([Cart_to_xx[0], Cart_to_xx[1], Cart_to_xx[2]],
                       ["xx[0]", "const REAL target_th", "xx[2]"], "returnstring", params="includebraces=False,preindent=1")
        body += "       xx[1] = NewtonRaphson_get_xx1_from_th(params, target_th);\n"
    else:
        body = outputC([Cart_to_xx[0], Cart_to_xx[1], Cart_to_xx[2]],
                       ["xx[0]", "xx[1]", "xx[2]"], "returnstring", params="includebraces=False,preindent=1")

    body += """
    // Then find the nearest index (i0,i1,i2) on underlying grid to (x,y,z)
    // Recall that:
    // xx[0][j] = xxmin[0] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx0"""+gridsuffix+"""; // Cell-centered grid.
    //   --> j = (int) ( (xx[0][j] - xxmin[0]) / params->dxx0"""+gridsuffix+""" + (1.0/2.0) + NGHOSTS )
    Cart_to_i0i1i2[0] = (int)( ( xx[0] - ("""+str(xxmin[0])+""") ) / params->dxx0"""+gridsuffix+""" + (1.0/2.0) + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
    Cart_to_i0i1i2[1] = (int)( ( xx[1] - ("""+str(xxmin[1])+""") ) / params->dxx1"""+gridsuffix+""" + (1.0/2.0) + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
    Cart_to_i0i1i2[2] = (int)( ( xx[2] - ("""+str(xxmin[2])+""") ) / params->dxx2"""+gridsuffix+""" + (1.0/2.0) + NGHOSTS - 0.5 ); // Account for (int) typecast rounding down
"""
    add_to_Cfunction_dict(
        includes=["stdio.h", "math.h", "stdlib.h",
                  os.path.join(rel_path_to_Cparams, "NRPy_basic_defines.h"),
                  os.path.join(rel_path_to_Cparams, "declare_Cparameters_struct.h")],
        prefunc=prefunc,
        desc   =desc,
        type   ="void",
        name   =name,
        params =params,
        preloop=preloop,
        body   =body,
        rel_path_to_Cparams=rel_path_to_Cparams)

    
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
        elif CoordSystem =="SinhCartesian":
            coordparams += """
params.AMPLX  = domain_size;
params.SINHWX = sinh_width;
params.AMPLY  = domain_size;
params.SINHWY = sinh_width;
params.AMPLZ  = domain_size;
params.SINHWZ = sinh_width;\n"""
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


############################
## TO BE DEPRECATED:
def set_Nxx_dxx_invdx_params__and__xx_h(outdir=".",grid_centering="cell"):
    if grid_centering not in ('cell', 'vertex'):
        print("rfm.set_Nxx_dxx_invdx_params__and__xx_h(): grid_centering = \""+grid_centering+"\" not supported!")
        sys.exit(1)

    with open(os.path.join(outdir,"set_Nxx_dxx_invdx_params__and__xx.h"),"w") as file:
        file.write(r"""
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
    params->invdx2 = 1.0/params->dxx2;\n""")
        # The following capability was suggested by Terrence Pierre Jacques (Thanks!)
        cell_offset = "(1.0/2.0)" # default cell-centered offset
        cell_comment = "Cell-centered grid."
        if grid_centering == "vertex":
            cell_offset = "0.0"
            cell_comment = "Vertex-centered grid."
        file.write("""
    // Now that params.dxx{0,1,2} and params.invdxx{0,1,2} have been set,
    // Step 0d.iii: Set up uniform coordinate grids
    xx[0] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);
    for(int j=0;j<Nxx_plus_2NGHOSTS0;j++)
        xx[0][j] = xxmin[0] + ((REAL)(j-NGHOSTS) + """+cell_offset+""")*params->dxx0; // """+cell_comment+"""
    xx[1] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS1);
    for(int j=0;j<Nxx_plus_2NGHOSTS1;j++)
        xx[1][j] = xxmin[1] + ((REAL)(j-NGHOSTS) + """+cell_offset+""")*params->dxx1; // """+cell_comment+"""
    xx[2] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS2);
    for(int j=0;j<Nxx_plus_2NGHOSTS2;j++)
        xx[2][j] = xxmin[2] + ((REAL)(j-NGHOSTS) + """+cell_offset+""")*params->dxx2; // """+cell_comment+"""
    //fprintf(stderr,"hey inside setxx: %e %e %e | %e %e\\n",xxmin[0],xxmin[1],xxmin[2],xx[0][0],params->dxx0);
}
""")

def xx_to_Cart_h(funcname,cparamsloc,outfile):
    # Arbitrary-coordinate NRPy+ file output, Part 1: output the conversion from (x0,x1,x2) to Cartesian (x,y,z)
    Cout = outputC([xx_to_Cart[0],xx_to_Cart[1],xx_to_Cart[2]],
                   ["xCart[0]","xCart[1]","xCart[2]"],
                   "returnstring",params="preindent=1")

    with open(outfile, "w") as file:
        file.write("""
static inline void """+funcname+"""(const paramstruct *restrict params, REAL *restrict xx[3],const int i0,const int i1,const int i2, REAL xCart[3]) {
#include """+"\""+cparamsloc+"\""+"""
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];\n"""+Cout+"}\n")

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
