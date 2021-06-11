# Gridfunction registration for a massless scalar field

# Author: Leonardo R. Werneck
#         wernecklr **at** gmail **dot* com

# This NRPy+ module is used internally by the other ScalarField NRPy+ modules

import sympy as sp
import grid as gri

def declare_scalar_field_gridfunctions_if_not_declared_already():
    # Step 2: Register all needed BSSN gridfunctions.
    
    global sf, sfM

    # Step 2.a: First check to see if this function has already been called.
    #   If so, do not register the gridfunctions again!
    for i in range(len(gri.glb_gridfcs_list)):
        if "sf" in gri.glb_gridfcs_list[i].name:
            sf, sfM = sp.symbols('sf sfM', real=True)
            return sf, sfM

    # Step 2.b: Register indexed quantities, using ixp.register_... functions
    sf, sfM = gri.register_gridfunctions("EVOL", ["sf", "sfM"])
    return sf, sfM
