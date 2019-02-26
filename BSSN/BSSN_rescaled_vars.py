# This module registers rescaled BSSN variables as gridfunctions

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step P1: import all needed modules from NRPy+:
import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import reference_metric as rfm

def declare_BSSN_rescaled_gridfunctions_if_not_declared_already():
    # Step 0: First check to see if this function has already been called.
    #   If so, do not register the gridfunctions again!
    for i in range(len(gri.glb_gridfcs_list)):
        if "hDD00" in gri.glb_gridfcs_list[i].name:
            return

    # Step 1: Declare as globals all quantities declared in this function.
    global hDD,cf,aDD,trK,lambdaU,alpha,vetU,betU

    # Step 2: Register all needed *evolved* gridfunctions.
    # Step 2a: Register indexed quantities, using ixp.register_... functions
    hDD = ixp.register_gridfunctions_for_single_rank2("EVOL", "hDD", "sym01")
    aDD = ixp.register_gridfunctions_for_single_rank2("EVOL", "aDD", "sym01")
    lambdaU = ixp.register_gridfunctions_for_single_rank1("EVOL", "lambdaU")
    vetU = ixp.register_gridfunctions_for_single_rank1("EVOL", "vetU")
    betU = ixp.register_gridfunctions_for_single_rank1("EVOL", "betU")
    # Step 2b: Register scalar quantities, using gri.register_gridfunctions()
    trK, cf, alpha = gri.register_gridfunctions("EVOL",["trK", "cf", "alpha"])