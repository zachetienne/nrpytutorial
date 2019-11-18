# This module registers rescaled BSSN T^{mu nu} source term variables
#    as AUX (i.e., not evolved) gridfunctions

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step P1: import all needed modules from NRPy+:
import indexedexp as ixp          # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import grid as gri                # NRPy+: Functions having to do with numerical grids

def define_BSSN_T4UUmunu_rescaled_source_terms():
    # Step 0: First check to see if this function has already been called.
    #   If so, do not register the gridfunctions again!
    for i in range(len(gri.glb_gridfcs_list)):
        if "sDD00" in gri.glb_gridfcs_list[i].name:
            return

    # Step 1: Declare as globals all quantities declared in this function.
    global rho,S,sD,sDD

    # Step 2: Register all needed *evolved* gridfunctions.
    # Step 2a: Register indexed quantities, using ixp.register_... functions
    sDD = ixp.register_gridfunctions_for_single_rank2("AUX", "sDD", "sym01")
    sD  = ixp.register_gridfunctions_for_single_rank1("AUX", "sD")
    # Step 2b: Register scalar quantities, using gri.register_gridfunctions()
    rho, S = gri.register_gridfunctions("AUX",["rho","S"])
