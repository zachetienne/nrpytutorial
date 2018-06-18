# reference_metric.py: Define reference metric ghatDD, rescaling "matrix" ReDD, and rescaling "vector" ReU
# reference_metric.py: Define reference metric ghatDD, rescaling "matrix" ReDD, and rescaling "vector" ReU

import time
import sympy as sp

import NRPy_param_funcs as par
# from outputC import *
import indexedexp as ixp

# Step 0a: Initialize parameters
thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "CoordSystem", "SymTP"))

# Step 0b: Declare global variables
#xx0,xx1,xx2,xx3 = par.Cparameters("REALARRAY",thismodule,['xx0','xx1','xx2','xx3'])
# Assuming xx's components are *positive* makes
#   certain simplifications possible later.
xx = ixp.declarerank1("xx",DIM=4)
xxCart = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s
xxSph  = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s
scalefactor_orthog = ixp.zerorank1(DIM=4) # Must be set in terms of xx[]s

xxmin = []
xxmax = []

def reference_metric():
    CoordSystem = par.parval_from_str("reference_metric::CoordSystem")

    # Set up hatted metric tensor, rescaling matrix, and rescaling vector
    if CoordSystem == "Spherical" or CoordSystem == "SinhSpherical" or CoordSystem == "SinhSphericalv2":
        # Assuming the spherical radial & theta coordinates
        #   are positive makes nice simplifications of
        #   unit vectors possible.
        xx[0] = sp.symbols("xx0", positive=True)
        xx[1] = sp.symbols("xx1", positive=True)

        r  = xx[0]
        th = xx[1]
        ph = xx[2]

        if CoordSystem == "Spherical":
            RMAX = par.Cparameters("REAL", thismodule, ["RMAX"])
            xxmin = ["0.0", "0.0", "0.0"]
            xxmax = ["params.RMAX", "M_PI", "2.0*M_PI"]

        elif CoordSystem == "SinhSpherical" or CoordSystem == "SinhSphericalv2":
            AMPL, SINHW = par.Cparameters("REAL",thismodule,["AMPL","SINHW"])
    
            xxmin = ["0.0", "0.0", "0.0"]
            xxmax = ["1.0", "M_PI", "2.0*M_PI"]
    
            # Set SinhSpherical radial coordinate by default; overwrite later if CoordSystem == "SinhSphericalv2".
            r = AMPL * (sp.exp(xx[0] / SINHW) - sp.exp(-xx[0] / SINHW)) / \
                       (sp.exp(1 / SINHW) - sp.exp(-1 / SINHW))
            # SinhSphericalv2 adds the parameter "const_dr", which allows for a region near xx[0]=0 to have
            # constant radial resolution of const_dr, provided the sinh() term does not dominate near xx[0]=0.
            if CoordSystem == "SinhSphericalv2":
                const_dr = par.Cparameters("REAL",thismodule,["const_dr"])
                r = AMPL*( const_dr*xx[0] + (sp.exp(xx[0] / SINHW) - sp.exp(-xx[0] / SINHW)) /
                           (sp.exp(1 / SINHW) - sp.exp(-1 / SINHW)) )

        # xxhat = sp.Matrix([[sp.sin(xxSph[1])*sp.cos(xxSph[2]), sp.sin(xxSph[1])*sp.sin(xxSph[2]), sp.cos(xxSph[1])],
        #                    [sp.cos(xxSph[1])*sp.cos(xxSph[2]), sp.cos(xxSph[1])*sp.sin(xxSph[2]),-sp.sin(xxSph[1])],
        #                    [-sp.sin(xxSph[2]),           sp.cos(xxSph[2]),            0         ]])

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

    elif CoordSystem == "Cylindrical" or CoordSystem == "SinhCylindrical" or CoordSystem == "SinhCylindricalv2":
        # Assuming the cylindrical radial coordinate
        #   is positive makes nice simplifications of
        #   unit vectors possible.
        xx[0] = sp.symbols("xx0", positive=True)

        RHOCYL = xx[0]
        PHICYL = xx[1]
        ZCYL   = xx[2]

        if CoordSystem == "Cylindrical":
            RHOMAX,ZMIN,ZMAX = par.Cparameters("REAL",thismodule,["RHOMAX","ZMIN","ZMAX"])
            xxmin = ["0.0", "0.0", "params.ZMIN"]
            xxmax = ["params.RHOMAX", "2.0*M_PI", "params.ZMAX"]
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
    
            xxmin = ["0.0","0.0","-1.0"]
            xxmax = ["1.0","2.0*M_PI","1.0"]
    
        # xxhat = sp.Matrix([[ sp.cos(PHICYL), sp.sin(PHICYL), 0],
        #                    [-sp.sin(PHICYL), sp.cos(PHICYL), 0],
        #                    [0,               0,              1]])
    
        xxCart[0] = RHOCYL*sp.cos(PHICYL)
        xxCart[1] = RHOCYL*sp.sin(PHICYL)
        xxCart[2] = ZCYL
    
        xxSph[0] = sp.sqrt(RHOCYL**2 + ZCYL**2)
        xxSph[1] = sp.cos(ZCYL / xxSph[0])
        xxSph[2] = PHICYL
    
        scalefactor_orthog[0] = sp.diff(RHOCYL,xx[0])
        scalefactor_orthog[1] = RHOCYL
        scalefactor_orthog[2] = sp.diff(ZCYL,xx[2])
    
    elif CoordSystem == "SymTP" or CoordSystem == "SinhSymTP":
        var1, var2= sp.symbols('var1 var2',positive=True)
        bScale, AW, AA, AMAX, RHOMAX, ZMIN, ZMAX = par.Cparameters("REAL",thismodule,["bScale","AW","AA","AMAX","RHOMAX","ZMIN","ZMAX"])

        # Assuming xx0, xx1, and bScale
        #   are positive makes nice simplifications of
        #   unit vectors possible.
        xx[0],xx[1],bScale = sp.symbols("xx0 xx1 bScale", positive=True)

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
    
        # xxhat = sp.Matrix([[sp.sin(xx[1]) * sp.cos(xx[2]) * var2 / var1, sp.sin(xx[1]) * sp.sin(xx[2]) * var2 / var1, AA * sp.cos(xx[1]) / var1],
        #                    [AA * sp.cos(xx[1]) * sp.cos(xx[2]) / var1,   AA * sp.cos(xx[1]) * sp.sin(xx[2]) / var1,  -sp.sin(xx[1]) * var2 / var1],
        #                    [-sp.sin(xx[2]),                           sp.cos(xx[2]),                            0                       ]])
    
        xxCart[0] = AA  *sp.sin(xx[1])*sp.cos(xx[2])
        xxCart[1] = AA  *sp.sin(xx[1])*sp.sin(xx[2])
        xxCart[2] = var2*sp.cos(xx[1])
    
        xxSph[0] = sp.sqrt(RHOSYMTP**2 + ZSYMTP**2)
        xxSph[1] = sp.cos(ZSYMTP / xxSph[0])
        xxSph[2] = PHSYMTP
    
        scalefactor_orthog[0] = sp.diff(AA,xx[0]) * var1 / var2
        scalefactor_orthog[1] = var1
        scalefactor_orthog[2] = AA * sp.sin(xx[1])
    
    elif CoordSystem == "Cartesian":
        xmin, xmax, ymin, ymax, zmin, zmax = par.Cparameters("REAL",thismodule,["xmin","xmax","ymin","ymax","zmin","zmax"])
        xxmin = ["params.xmin", "params.ymin", "params.zmin"]
        xxmax = ["params.xmax", "params.ymax", "params.zmax"]
    
        # xxhat = sp.Matrix([[sp.sympify(1), 0, 0],
        #                    [0, sp.sympify(1), 0],
        #                    [0, 0, sp.sympify(1)]])
    
        xxCart[0] = xx[0]
        xxCart[1] = xx[1]
        xxCart[2] = xx[2]

        xxSph[0] = sp.sqrt(xx[0] ** 2 + xx[1] ** 2 + xx[2] ** 2)
        xxSph[1] = sp.cos(xx[2] / xxSph[0])
        xxSph[2] = sp.atan2(xx[1], xx[0])

        scalefactor_orthog[0] = sp.sympify(1)
        scalefactor_orthog[1] = sp.sympify(1)
        scalefactor_orthog[2] = sp.sympify(1)
    else:
        print("CoordSystem == " + CoordSystem + " is not supported.")
        exit(1)

def UnitVectors3D():
    xxhats3D = ixp.zerorank2(3)
    for i in range(3):
        norm = 0
        for j in range(3):
            numer = sp.diff(xxCart[j],xx[i])
            norm += numer**2
            xxhats3D[i][j] = numer
        norm = sp.simplify(sp.sqrt(sp.simplify(norm)))
        for j in range(3):
            xxhats3D[i][j] /= norm
    return xxhats3D

# for i in range(3):
#     set_parameter("REAL", "x"+str(i+1)+"min", xxmin[i], params_type, params_varname, params_value)
#     set_parameter("REAL", "x"+str(i+1)+"max", xxmax[i], params_type, params_varname, params_value)


#     # -={ \gammahat_{ij}: output to codegen_output/NRPy_gammahatDD_*.h }=-
#
# # Set precision to, e.g., double or long double
# PRECISION = get_parameter_value("PRECISION", params_varname, params_value)
#
# # Jacobian matrix needed for ADM integrands inside BSSN_Diagnostics.py
# if RunMode == "Diagnostics" or RunMode == "Everything":
#     NRPy_file_output("common_functions/reference_metric/NRPy_codegen/gammahatDD.h", [], [], [],
#                      protected_varnames, [],
#                      [ghatDD[0][0], "const " + PRECISION + " gammahatDD00",
#                       ghatDD[0][1], "const " + PRECISION + " gammahatDD01",
#                       ghatDD[0][2], "const " + PRECISION + " gammahatDD02",
#                       ghatDD[1][1], "const " + PRECISION + " gammahatDD11",
#                       ghatDD[1][2], "const " + PRECISION + " gammahatDD12",
#                       ghatDD[2][2], "const " + PRECISION + " gammahatDD22"])
#
#     #ROOTDIR = root_directory(params_str,params_val)+"/"
#     NRPy_file_output("common_functions/reference_metric/NRPy_codegen/dist_from_origin.h", [],[],[], protected_varnames,[],[r,"const REAL dist_from_origin"])
#
#     NRPy_file_output("common_functions/reference_metric/NRPy_codegen/scalefactor_orthog.h", [],[],[],protected_varnames,[],
#                      [scalefactor_orthog[0],"scalefactor_orthog[0]",scalefactor_orthog[1], "scalefactor_orthog[1]",scalefactor_orthog[2], "scalefactor_orthog[2]",])
#
#     NRPy_file_output("common_functions/reference_metric/NRPy_codegen/xyz.h", [],[],[], protected_varnames,
#                      [],[xCart,"Cartxyz[0][idx]",yCart,"Cartxyz[1][idx]",zCart,"Cartxyz[2][idx]"])
#
#     NRPy_file_output("common_functions/reference_metric/NRPy_codegen/xxCart.h", [],[],[], protected_varnames,
#                      [],[xCart,"xCart",yCart,"yCart",zCart,"zCart"])
#
#     drdx = [[sp.diff(xCart,xx[0]), sp.diff(xCart,xx[1]), sp.diff(xCart,xx[2])],
#             [sp.diff(yCart,xx[0]), sp.diff(yCart,xx[1]), sp.diff(yCart,xx[2])],
#             [sp.diff(zCart,xx[0]), sp.diff(zCart,xx[1]), sp.diff(zCart,xx[2])]]
#     #pprint(drdx)
#     dxdr = [[sp.sympify(0) for i in range(3)] for j in range(3)]
#     dummyDET = generic_matrix_inverter3xx[2](drdx, dxdr)
#
#     dydx = [sp.sympify(0) for i in range(18)]
#     dxdy = [sp.sympify(0) for i in range(18)]
#     counter = 0
#     for i in range(3):
#         for j in range(3):
#             dydx[counter] = (drdx[i][j])
#             dxdy[counter] = (dxdr[i][j])
#             counter += 1
#             dydx[counter] = "dydx["+str(i)+"]["+str(j)+"]"
#             dxdy[counter] = "dxdy["+str(i)+"]["+str(j)+"]"
#             counter += 1
#
#     NRPy_file_output("common_functions/reference_metric/NRPy_codegen/Jacobian_sp.Matrix.h", [],[],[], protected_varnames,[],dydx)
#     NRPy_file_output("common_functions/reference_metric/NRPy_codegen/Inverse_Jacobian_sp.Matrix.h", [],[],[], protected_varnames,[],dxdy)
#
#     unitVec = [sp.sympify(0) for i in range(18)]
#     counter = 0
#     for i in range(3):
#         for j in range(3):
#             unitVec[counter] = (xxhat[(i,j)])
#             counter += 1
#             unitVec[counter] = "x"+str(i+1)+"hat["+str(j)+"]"
#             counter += 1
#
#     NRPy_file_output("common_functions/reference_metric/NRPy_codegen/xxhat.h", [],[],[], protected_varnames,[],unitVec)
#
# stop1 = time.time()
# print("Defined reference metric quantities in \t\t\t" + str(round(stop1-start1,2)) + " seconds")
# # *****************************************************
