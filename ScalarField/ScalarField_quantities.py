# Gridfunction registration for a massless scalar field

# Author: Leonardo R. Werneck
#         wernecklr **at** gmail **dot* com

# This NRPy+ module is used internally by the other ScalarField NRPy+ modules

from sympy import symbols
from grid import glb_gridfcs_list
from grid import register_gridfunctions
from indexedexp import declarerank1
from indexedexp import declarerank2

def declare_scalar_field_gridfunctions_if_not_declared_already():

    # Step 1: Register all needed ScalarField gridfunctions
    # Step 1.a: First check to see if this function has already been called.
    #   If so, do not register the gridfunctions again!
    for i in range(len(glb_gridfcs_list)):
        if "sf" in glb_gridfcs_list[i].name:
            sf, sfM = symbols("sf sfM", real=True)
            return sf, sfM

    # Step 1.b: Register indexed quantities, using ixp.register_... functions
    sf, sfM = register_gridfunctions("EVOL", ["sf", "sfM"])
    return sf, sfM

def ScalarField_quantities():

    # Gridfunctions
    global sf, sfM
    sf, sfM = declare_scalar_field_gridfunctions_if_not_declared_already()

    # Regular derivatives
    global sf_dD, sf_dDD
    sf_dD    = declarerank1("sf_dD")
    sf_dDD   = declarerank2("sf_dDD", "sym01")

    # Upwinded derivatives
    global sf_dupD, sfM_dupD
    sf_dupD  = declarerank1("sf_dupD")
    sfM_dupD = declarerank1("sfM_dupD")

    # Kreiss-Oliger dissipation derivatives
    global sf_dKOD, sfM_dKOD
    sf_dKOD  = declarerank1("sf_dKOD")
    sfM_dKOD = declarerank1("sfM_dKOD")
