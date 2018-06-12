# reference_metric.py: Define reference metric ghatDD, rescaling "matrix" ReDD, and rescaling "vector" ReU

import time
import sympy as sp

import NRPy_param_funcs as par
from outputC import *



# Here we set protected variable names. By being listed here,
#  they are "protected" against being interpreted as gridfunctions.
protected_varnames = ["Nx1","Nx2","Nx3", "x1","x2","x3", "M_PI"]

t, x1, x2, x3 = symbols('t x1 x2 x3',real=True)

xx = [x1,x2,x3]
xCart, yCart, zCart = symbols('xCart yCart zCart',real=True)
r, th, ph = symbols('r th ph',real=True)
y1, y2, y3 = symbols('y1 y2 y3',real=True)

ReU    =   [ sympify(0) for i in range(3) ]
ReDD   =   [[ sympify(0) for i in range(3)] for j in range(3)]
ghatDD =   [[ sympify(0) for i in range(3)] for j in range(3)]

scalefactor_orthog = [sympify(0) for i in range(3)]

# Read needed parameters from params_str,params_val
CoordSystem = get_parameter_value("CoordSystem", params_varname, params_value)

# Nx123 = total number of points in each direction, EXCLUDING ghostzones. Input at runtime
runtime_parameter("int", "Nx1", runtime_params_type, runtime_params_varname)
runtime_parameter("int", "Nx2", runtime_params_type, runtime_params_varname)
runtime_parameter("int", "Nx3", runtime_params_type, runtime_params_varname)

x123min = []
x123max = []

# Set up hatted metric tensor, rescaling matrix, and rescaling vector
if CoordSystem == "Spherical" or CoordSystem == "SinhSpherical" or CoordSystem == "SinhSphericalv2":
    if CoordSystem == "SinhSpherical" or CoordSystem == "SinhSphericalv2":
        AMPL, SINHW = symbols('AMPL SINHW', positive=True)
        protected_varnames += ["AMPL", "SINHW"]

        runtime_parameter("REAL", "AMPL",  runtime_params_type, runtime_params_varname)
        runtime_parameter("REAL", "SINHW", runtime_params_type, runtime_params_varname)

        x123min = ["0.0", "0.0", "0.0"]
        x123max = ["1.0", "M_PI", "2.0*M_PI"]

        # Set SinhSpherical radial coordinate by default; overwrite later if CoordSystem == "SinhSphericalv2".
        r = AMPL * (exp(x1 / SINHW) - exp(-x1 / SINHW)) / (exp(1 / SINHW) - exp(-1 / SINHW))
        # SinhSphericalv2 adds the parameter "const_dr", which allows for a region near x1=0 to have
        # constant radial resolution of const_dr, provided the sinh() term does not dominate near x1=0.
        if CoordSystem == "SinhSphericalv2":
            const_dr = symbols('const_dr', positive=True)
            protected_varnames += ["const_dr"]
            runtime_parameter("REAL", "const_dr", runtime_params_type, runtime_params_varname)
            r = AMPL*( const_dr*x1 + (exp(x1 / SINHW) - exp(-x1 / SINHW)) / (exp(1 / SINHW) - exp(-1 / SINHW)) )
        th = x2
        ph = x3

    if CoordSystem == "Spherical":
        runtime_parameter("REAL", "RMAX", runtime_params_type, runtime_params_varname)
        r = x1
        th = x2
        ph = x3

        x123min = ["0.0", "0.0", "0.0"]
        x123max = ["params.RMAX", "M_PI", "2.0*M_PI"]

    xxhat = Matrix([[sin(th)*cos(ph), sin(th)*sin(ph), cos(th)],
                    [cos(th)*cos(ph), cos(th)*sin(ph), -sin(th)],
                    [-sin(ph), cos(ph), 0]])

    # The following is used by the Spherical-to-Cartesian Einstein Toolkit (ETK) layer.
    # It solves for x1,x2, and x3 in terms of the Cartesian coordinates x,y,z. The inversion of
    # SinhSphericalv2 has no analytic form, so we cannot use this layer, in its current form, with SinhSphericalv2.
    if (CoordSystem == "Spherical" or CoordSystem == "SinhSpherical") and \
            (RunMode == "Diagnostics" or RunMode == "Everything"):
        # Next define r,th,phi in terms of xCart, yCart, and zCart.
        #  This functionality is only used for interpolating to a Cartesian grid, to use EinsteinToolkit diagnostic thorns
        rCart = sqrt((xCart ** 2 + yCart ** 2 + zCart ** 2))

        if CoordSystem == "Spherical":
            x1Cart = rCart
        elif CoordSystem == "SinhSpherical":
            # r = AMPL/(2*sinh(1/SINHW)) * 2*sinh(x1/SINHW)
            # -> sinh(x1/SINHW) = r*sinh(1/SINHW)/AMPL
            # sinh^{-1}(x) = ln(x + sqrt(1+x^2))
            # -> x1 = SINHW * sinh^{-1}(r*sinh(1/SINHW)/AMPL)
            #       = SINHW * ln( (r*sinh(1/SINHW)/AMPL) + sqrt(1+(r*sinh(1/SINHW)/AMPL)^2)),
            # where r = sqrt(xCart**2 + yCart**2 + zCart**2), and sinh(1/SINHW) = (exp(1/SINHW) - exp(-1/SINHW))/2
            r_temp = rCart * ( (exp(1 / SINHW) - exp(-1 / SINHW))/2 )/AMPL
            x1Cart = SINHW * log(r_temp + sqrt(1 + r_temp**2))
        x2Cart = acos(zCart/rCart)
        x3Cart = atan2(yCart,xCart)

        NRPy_file_output("common_functions/reference_metric/NRPy_codegen/xx_in_terms_of_Cartxyz.h", [], [], [], protected_varnames+["xCart","yCart","zCart"],
                         [], [x1Cart, "xxCart[0]", x2Cart, "xxCart[1]", x3Cart, "xxCart[2]"])
        # atan2() outputs a value between -pi and +pi. However, theta is defined between 0 and pi, and phi between 0 & 2pi
        #  Thus we have to add a couple extra lines of code to ensure we are in bounds.
        with open("common_functions/reference_metric/NRPy_codegen/xx_in_terms_of_Cartxyz.h", "a") as output:
            output.write("if(xxCart[2]<0) xxCart[2]+=2.0*M_PI;\n")
            output.write("if(xxCart[0]==0 || xxCart[1]==0 || xxCart[1]==M_PI) {printf(\"ERROR: tried to interpolate to r*sin(theta)=0, which is a coordinate singularity\\n\");exit(1);}\n")

    # Now define xCart, yCart, and zCart in terms of x1,x2,x3.
    #   Note that the relation between r and x1 is not necessarily trivial in SinhSpherical coordinates. See above.
    xCart = r*sin(th)*cos(ph)
    yCart = r*sin(th)*sin(ph)
    zCart = r*cos(th)

    y1 = r
    y2 = th
    y3 = ph

    scalefactor_orthog[0] = diff(y1,x1)
    scalefactor_orthog[1] = y1
    scalefactor_orthog[2] = y1*sin(y2)
elif CoordSystem == "Cylindrical" or CoordSystem == "SinhCylindrical" or CoordSystem == "SinhCylindricalv2":
    if CoordSystem == "Cylindrical":
        runtime_parameter("REAL", "RHOMAX", runtime_params_type, runtime_params_varname)
        runtime_parameter("REAL", "ZMIN", runtime_params_type, runtime_params_varname)
        runtime_parameter("REAL", "ZMAX", runtime_params_type, runtime_params_varname)

        RHOCYL = x1
        # phi coordinate remains unchanged.
        PHICYL = x2
        ZCYL = x3

        x123min = ["0.0", "0.0", "params.ZMIN"]
        x123max = ["params.RHOMAX", "2.0*M_PI", "params.ZMAX"]
    elif CoordSystem == "SinhCylindrical" or CoordSystem == "SinhCylindricalv2":
        AMPLRHO, SINHWRHO, AMPLZ, SINHWZ = symbols('AMPLRHO SINHWRHO AMPLZ SINHWZ', positive=True)
        protected_varnames += ["AMPLRHO", "SINHWRHO", "AMPLZ", "SINHWZ"]

        runtime_parameter("REAL", "AMPLRHO", runtime_params_type, runtime_params_varname)
        runtime_parameter("REAL", "SINHWRHO", runtime_params_type, runtime_params_varname)
        runtime_parameter("REAL", "AMPLZ", runtime_params_type, runtime_params_varname)
        runtime_parameter("REAL", "SINHWZ", runtime_params_type, runtime_params_varname)

        # Set SinhCylindrical radial & z coordinates by default; overwrite later if CoordSystem == "SinhCylindricalv2".
        RHOCYL = AMPLRHO * (exp(x1 / SINHWRHO) - exp(-x1 / SINHWRHO)) / (exp(1 / SINHWRHO) - exp(-1 / SINHWRHO))
        # phi coordinate remains unchanged.
        PHICYL = x2
        ZCYL   = AMPLZ   * (exp(x3 / SINHWZ)   - exp(-x3 / SINHWZ))   / (exp(1 / SINHWZ)   - exp(-1 / SINHWZ))

        # SinhCylindricalv2 adds the parameters "const_drho", "const_dz", which allows for regions near x1=0
        # and x3=0 to have constant rho and z resolution of const_drho and const_dz, provided the sinh() terms
        # do not dominate near x1=0 and x3=0.
        if CoordSystem == "SinhCylindricalv2":
            const_drho, const_dz = symbols('const_drho const_dz', positive=True)
            protected_varnames += ["const_drho", "const_dz"]
            runtime_parameter("REAL", "const_drho", runtime_params_type, runtime_params_varname)
            runtime_parameter("REAL", "const_dz", runtime_params_type, runtime_params_varname)
            RHOCYL = AMPLRHO * ( const_drho*x1 + (exp(x1 / SINHWRHO) - exp(-x1 / SINHWRHO)) / (exp(1 / SINHWRHO) - exp(-1 / SINHWRHO)) )
            ZCYL   = AMPLZ   * ( const_dz  *x3 + (exp(x3 / SINHWZ  ) - exp(-x3 / SINHWZ  )) / (exp(1 / SINHWZ  ) - exp(-1 / SINHWZ  )) )

        x123min = ["0.0","0.0","-1.0"]
        x123max = ["1.0","2.0*M_PI","1.0"]

    xxhat = Matrix([[cos(PHICYL), sin(PHICYL), 0],
                    [-sin(PHICYL), cos(PHICYL), 0],
                    [0, 0, 1]])

    xCart = RHOCYL*cos(PHICYL)
    yCart = RHOCYL*sin(PHICYL)
    zCart = ZCYL

    r = sqrt(RHOCYL**2 + ZCYL**2)
    th = acos(ZCYL / r)
    ph = PHICYL

    y1 = RHOCYL
    y2 = PHICYL
    y3 = ZCYL

    scalefactor_orthog[0] = diff(y1,x1)
    scalefactor_orthog[1] = y1
    scalefactor_orthog[2] = diff(y3,x3)

elif CoordSystem == "SymTP" or CoordSystem == "SinhSymTP":
    bScale, AW, AA, var1, var2, AMAX = symbols('bScale AW AA var1 var2 AMAX',positive=True)
    protected_varnames += ["bScale", "AW"]
    RHOMAX, ZMIN, ZMAX = symbols('RHOMAX ZMIN ZMAX', positive=True)

    runtime_parameter("REAL", "AMAX",   runtime_params_type, runtime_params_varname)
    runtime_parameter("REAL", "bScale", runtime_params_type, runtime_params_varname)

    x123min = ["0.0","0.0","0.0"]
    x123max = ["params.AMAX","M_PI","2.0*M_PI"]

    AA = x1

    if CoordSystem == "SinhSymTP":
        AA = (exp(x1/AW)-exp(-x1/AW))/2

    var1 = sqrt(AA**2 + (bScale * sin(x2))**2)
    var2 = sqrt(AA**2 + bScale**2)

    RHOSYMTP = AA*sin(x2)
    PHSYMTP = x3
    ZSYMTP = var2*cos(x2)

    xxhat = Matrix([[sin(x2) * cos(x3) * var2 / var1, sin(x2) * sin(x3) * var2 / var1, AA * cos(x2) / var1],
                    [AA * cos(x2) * cos(x3) / var1, AA * cos(x2) * sin(x3) / var1, -sin(x2) * var2 / var1],
                    [-sin(x3), cos(x3), 0]])

    xCart = AA  *sin(x2)*cos(x3)
    yCart = AA  *sin(x2)*sin(x3)
    zCart = var2*cos(x2)

    r = sqrt(RHOSYMTP**2 + ZSYMTP**2)
    th = acos(ZSYMTP / r)
    ph = PHSYMTP

    y1 = RHOSYMTP
    y2 = PHSYMTP
    y3 = ZSYMTP

    scalefactor_orthog[0] = diff(AA,x1) * var1 / var2
    scalefactor_orthog[1] = var1
    scalefactor_orthog[2] = AA * sin(x2)

elif CoordSystem == "Cartesian":
    xmin,xmax,ymin,ymax,zmin,zmax = symbols('xmin xmax ymin ymax zmin zmax',real=True)
    protected_varnames += ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]
    x123min = ["params.xmin", "params.ymin", "params.zmin"]
    x123max = ["params.xmax", "params.ymax", "params.zmax"]

    r = sqrt(x1**2 + x2**2 + x3**2)
    th = acos(x3/r)
    ph = atan2(x2,x1)

    runtime_parameter("REAL", "xmin", runtime_params_type, runtime_params_varname)
    runtime_parameter("REAL", "xmax", runtime_params_type, runtime_params_varname)
    runtime_parameter("REAL", "ymin", runtime_params_type, runtime_params_varname)
    runtime_parameter("REAL", "ymax", runtime_params_type, runtime_params_varname)
    runtime_parameter("REAL", "zmin", runtime_params_type, runtime_params_varname)
    runtime_parameter("REAL", "zmax", runtime_params_type, runtime_params_varname)

    xxhat = Matrix([[sympify(1), 0, 0],
                    [0, sympify(1), 0],
                    [0, 0, sympify(1)]])

    x1Cart = x1
    x2Cart = x2
    x3Cart = x3
    NRPy_file_output("common_functions/reference_metric/NRPy_codegen/xx_in_terms_of_Cartxyz.h", [], [], [], protected_varnames+["xCart","yCart","zCart"],
                     [], [x1Cart, "xxCart[0]", x2Cart, "xxCart[1]", x3Cart, "xxCart[2]"])

    # Now define xCart, yCart, and zCart in terms of x1,x2,x3.
    #   Note that the relation between r and x1 is not necessarily trivial in SinhSpherical coordinates. See above.
    xCart = x1
    yCart = x2
    zCart = x3

    y1 = x1
    y2 = x2
    y3 = x3

    scalefactor_orthog[0] = sympify(1)
    scalefactor_orthog[1] = sympify(1)
    scalefactor_orthog[2] = sympify(1)

else:
    print("CoordSystem == " + CoordSystem + " is not supported.")
    exit(1)

for i in range(3):
    set_parameter("REAL", "x"+str(i+1)+"min", x123min[i], params_type, params_varname, params_value)
    set_parameter("REAL", "x"+str(i+1)+"max", x123max[i], params_type, params_varname, params_value)

for i in range(3):
    ghatDD[i][i] = scalefactor_orthog[i]**2
    ReU[i] = 1/scalefactor_orthog[i]
    for j in range(3):
        ReDD[i][j] = scalefactor_orthog[i]*scalefactor_orthog[j]

    # -={ \gammahat_{ij}: output to codegen_output/NRPy_gammahatDD_*.h }=-

# Set precision to, e.g., double or long double
PRECISION = get_parameter_value("PRECISION", params_varname, params_value)

# Jacobian matrix needed for ADM integrands inside BSSN_Diagnostics.py
if RunMode == "Diagnostics" or RunMode == "Everything":
    NRPy_file_output("common_functions/reference_metric/NRPy_codegen/gammahatDD.h", [], [], [],
                     protected_varnames, [],
                     [ghatDD[0][0], "const " + PRECISION + " gammahatDD00",
                      ghatDD[0][1], "const " + PRECISION + " gammahatDD01",
                      ghatDD[0][2], "const " + PRECISION + " gammahatDD02",
                      ghatDD[1][1], "const " + PRECISION + " gammahatDD11",
                      ghatDD[1][2], "const " + PRECISION + " gammahatDD12",
                      ghatDD[2][2], "const " + PRECISION + " gammahatDD22"])

    #ROOTDIR = root_directory(params_str,params_val)+"/"
    NRPy_file_output("common_functions/reference_metric/NRPy_codegen/dist_from_origin.h", [],[],[], protected_varnames,[],[r,"const REAL dist_from_origin"])

    NRPy_file_output("common_functions/reference_metric/NRPy_codegen/scalefactor_orthog.h", [],[],[],protected_varnames,[],
                     [scalefactor_orthog[0],"scalefactor_orthog[0]",scalefactor_orthog[1], "scalefactor_orthog[1]",scalefactor_orthog[2], "scalefactor_orthog[2]",])

    NRPy_file_output("common_functions/reference_metric/NRPy_codegen/xyz.h", [],[],[], protected_varnames,
                     [],[xCart,"Cartxyz[0][idx]",yCart,"Cartxyz[1][idx]",zCart,"Cartxyz[2][idx]"])

    NRPy_file_output("common_functions/reference_metric/NRPy_codegen/xxCart.h", [],[],[], protected_varnames,
                     [],[xCart,"xCart",yCart,"yCart",zCart,"zCart"])

    drdx = [[diff(xCart,xx[0]), diff(xCart,xx[1]), diff(xCart,xx[2])],
            [diff(yCart,xx[0]), diff(yCart,xx[1]), diff(yCart,xx[2])],
            [diff(zCart,xx[0]), diff(zCart,xx[1]), diff(zCart,xx[2])]]
    #pprint(drdx)
    dxdr = [[sympify(0) for i in range(3)] for j in range(3)]
    dummyDET = generic_matrix_inverter3x3(drdx, dxdr)

    dydx = [sympify(0) for i in range(18)]
    dxdy = [sympify(0) for i in range(18)]
    counter = 0
    for i in range(3):
        for j in range(3):
            dydx[counter] = (drdx[i][j])
            dxdy[counter] = (dxdr[i][j])
            counter += 1
            dydx[counter] = "dydx["+str(i)+"]["+str(j)+"]"
            dxdy[counter] = "dxdy["+str(i)+"]["+str(j)+"]"
            counter += 1

    NRPy_file_output("common_functions/reference_metric/NRPy_codegen/Jacobian_Matrix.h", [],[],[], protected_varnames,[],dydx)
    NRPy_file_output("common_functions/reference_metric/NRPy_codegen/Inverse_Jacobian_Matrix.h", [],[],[], protected_varnames,[],dxdy)

    unitVec = [sympify(0) for i in range(18)]
    counter = 0
    for i in range(3):
        for j in range(3):
            unitVec[counter] = (xxhat[(i,j)])
            counter += 1
            unitVec[counter] = "x"+str(i+1)+"hat["+str(j)+"]"
            counter += 1

    NRPy_file_output("common_functions/reference_metric/NRPy_codegen/xxhat.h", [],[],[], protected_varnames,[],unitVec)

stop1 = time.time()
print("Defined reference metric quantities in \t\t\t" + str(round(stop1-start1,2)) + " seconds")
# *****************************************************
