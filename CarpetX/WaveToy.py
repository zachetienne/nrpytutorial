import os
import sys
import re
import grid
from cactusthorn import CactusThorn, loop
from sympy import sympify, sin, cos, pi
import NRPy_param_funcs as par
from subprocess import call
import numpy as np
from fstr import f

pre_kernel = """
#include <iostream>
#include <cmath>

const int nx = 10, ny = 10;

using std::cos;
using std::sin;

struct PointDesc {
    static constexpr int DI[3]={0,0,0};
    double x,y;
    int I=0;
};

struct GF {
    double data;
    double& operator()(int n) { return data; }
    double& operator()(int n,int p) { return data; }
};

double cctk_delta_space[3] = {.5, .4};
inline double CCTK_DELTA_SPACE(int n) { return cctk_delta_space[n]; }

int main() {
    int VVC_index = 0;
    int CCC_index = 0;
    int VVC_tmp_index = 0;
    int VVC_layout = 0;
    GF uuGF, vvGF, anaGF, rhs_uuGF, rhs_vvGF, tmp1v, tmp0v, regrid_error;
    PointDesc p;
    double cctk_time = 0;
    double wave_speed = .4;

    for(int i=0;i<nx;i++) {
      for(int j=0;j<ny;j++) {
        const double invdx0 = 1 / CCTK_DELTA_SPACE(0);
        const double invdx1 = 1 / CCTK_DELTA_SPACE(1);
        uuGF(0) = .25 + .01*i - .015*j;
        vvGF(0) = .33 - .015*i + .015*j;
    // Begin test.cc
"""
post_kernel="""
    // End test.cc
        std::cout << anaGF(0) << std::endl;
      }
    }
    return 0;
}"""

def before_main():
    # Init environment
    os.environ["CACTUS_HOME"]="Cactus"
    os.environ["CACTUS_DRIVER"]="CarpetX"

    # Trick CactusThorn into believing this is a real install
    try:
        os.makedirs("Cactus/arrangements/CarpetX/CarpetX") #,exist_ok=True)
    except OSError as oe:
        pass

def after_main_fn(fn):
    save = False
    kernel = ""
    with open(fn, "r") as fd:
        for line in fd.readlines():
            if "Begin NRPy+ Kernel" in line:
                save = True
            elif "End NRPy+ Kernel" in line:
                save = False
            elif save:
                kernel += line
    if kernel == "":
        return None

    fc = os.path.join("Cactus","test.cc")
    with open(fc, "w") as fd:
        fd.write(f("""
   {pre_kernel}
   {kernel}
   {post_kernel}
"""))
    test = os.path.join("Cactus","test")
    testc = test + ".cc"
    testo = test + "-out.txt"
    r = call(["g++","-std=c++17","-D_USE_MATH_DEFINES","-o",test,testc])
    if r != 0:
        raise Exception(f("Compile failure of {testc} while processing "+fn))
    with open(testo,"w") as fd:
        call([test],stdout=fd)
    data = np.genfromtxt(testo,encoding="ascii")
    g = re.match(r'.*\bwave_(.*)\.cc', fn)
    if g is not None:
        print("Setting:",g.group(1))
        globals()[g.group(1)] = [float(ff) for ff in data]
        return data

def after_main():
    dn = os.path.join("Cactus","arrangements","TestOne","WaveToyNRPy","src")
    da = []
    for fn in os.listdir(dn):
        data = after_main_fn(os.path.join(dn, fn))
        if data is not None:
            print("File:",fn)
            da += [data]
    return da

def main():
    global par

    # Current options are Carpet and CarpetX
    grid.ET_driver = os.environ.get("CACTUS_DRIVER","Carpet")

    thorn = CactusThorn("TestOne","WaveToyNRPy")

    #FD_order = thorn.declare_param('FD_order',default=2,vmin=2,vmax=8,doc="The finite difference order")

    wave_speed = thorn.declare_param('wave_speed',default=1,vmin=.1,vmax=100,doc="The speed of the wave")
    x0 = thorn.declare_param('x0',default=0,vmin=-100,vmax=100,doc="The x pos of the wave")
    y0 = thorn.declare_param('y0',default=0,vmin=-100,vmax=100,doc="The y pos of the wave")
    z0 = thorn.declare_param('z0',default=0,vmin=-100,vmax=100,doc="The z pos of the wave")
    zero = thorn.declare_param('zero',default=0,vmin=0,vmax=0,doc="zero")

    centering='CVV'

    # AUXEVOL needed for the evo, can be freed after evaluating rhs (1 time level)
    # AUX uu_rhs (1 time level)
    # EVOL evolved gfs (3 time levels)
    uu_rhs, vv_rhs = thorn.register_gridfunctions("AUX", ["rhs_uu", "rhs_vv"], centering=centering)
    uu, vv = thorn.register_gridfunctions("EVOL", ["uu", "vv"], centering=centering)
    ana = thorn.register_gridfunctions("AUX", ["ana"], centering=centering)
    cctk_time = par.Cparameters("CCTK_REAL","Cactus",["cctk_time"],0)
    if grid.ET_driver == "CarpetX":
        tmp0, tmp1 = thorn.register_gridfunctions("TILE_TMP",["tmp0v","tmp1v"],centering=centering)
        regrid_error = thorn.get_regrid_error()
    x,y,z = thorn.get_xyz()

    from outputC import lhrh
    import indexedexp as ixp
    import NRPy_param_funcs as par

    FD_order = 2
    par.set_parval_from_str("finite_difference::FD_CENTDERIVS_ORDER",FD_order)

    uu_dDD = ixp.declarerank2("uu_dDD","sym01")

    if grid.ET_driver == "CarpetX":
        evol_eqns = [
            lhrh(lhs=tmp0, rhs=uu_dDD[0][0]),
            lhrh(lhs=tmp1, rhs=uu_dDD[1][1]),
            loop,
            lhrh(lhs=uu_rhs, rhs=vv),
            lhrh(lhs=vv_rhs, rhs=wave_speed**2*(tmp0 + tmp1))
        ]
    else:
        # Version of evolution equations without temporaries
        evol_eqns = [
            lhrh(lhs=uu_rhs, rhs=vv),
            lhrh(lhs=vv_rhs, rhs=wave_speed**2*(uu_dDD[0][0] + uu_dDD[1][1]))
        ]

    k = sympify(pi/20)
    toff = sympify(pi/2)
    sq2 = sympify(2**.5)

    init_eqns = [
        lhrh(lhs=vv, rhs=sympify(0)),
        lhrh(lhs=uu, rhs=sin(k*(x))*sin(k*(y)))
    ]

    # What is omega?
    # omega / k = wave_speed

    anal_eqns = [
        lhrh(lhs=ana, rhs=sin(k*x)*sin(k*y)*sin(cctk_time*sq2*wave_speed*k+toff)-uu)
    ]

    if grid.ET_driver == "CarpetX":
        thorn.add_func("refine",
            body=[lhrh(lhs=regrid_error, rhs=10/((x-20)**2 + (y-20)**2))],
            doc="do the regrid",
            schedule_bin="ODESolvers_EstimateError")

    # access a variable with a different centering using interpolation
    # looping cell-centered, access vertex-centered, not vice-versa
    # all rhs variables should have the same centering
    # wave toy with fluxes, fluxes are faces
    # schedule something in post-regrid, apply bc's
    thorn.add_func("wave_init",
        body=init_eqns,
        where='everywhere',
        schedule_bin='initial',
        doc='Do the wave init',
        centering=centering)

    #thorn.add_func("wave_bound",
    #    body=bound_eqns,
    #    where='boundary',
    #    schedule_bin='RHS after wave_evol',
    #    doc='Do the b/c',
    #    centering=centering)

    thorn.add_func("wave_evol",
        body=evol_eqns,
        where='interior',
        schedule_bin='RHS',
        doc='Do the wave evol',
        centering=centering)

    thorn.add_func("wave_anal",
        body=anal_eqns,
        where='everywhere',
        schedule_bin='Analysis',
        doc='Check the result',
        centering=centering)

    assert "CACTUS_HOME" in os.environ, "Please set the CACTUS_HOME variable to point to your Cactus installation"
    cactus_home = os.environ["CACTUS_HOME"]
    cactus_sim = os.environ.get("CACTUS_SIM","sim")
    cactus_thornlist = os.environ.get("CACTUS_THORNLIST", None)

    #cactus_home = "/project/sbrandt/release/Cactus"
    thorn.generate(cactus_home,cactus_config=cactus_sim,cactus_thornlist=cactus_thornlist)

def run_all():
    before_main()
    main()
    after_main()

if __name__ == "__main__":
    main()
