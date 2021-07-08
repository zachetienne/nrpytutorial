# As documented in the NRPy+ tutorial module
#   Tutorial-RK_Butcher_Table_Generating_C_Code.ipynb,
#   this module will produce the required C codes for
#   allocating required memory Method of Lines (MoL) timestepping,
#   implementing MoL timestepping, and deallocating memory

# Authors: Brandon Clark
#          Zachariah B. Etienne
#          zachetie **at** gmail **dot* com

# Step 1: Initialize needed Python/NRPy+ modules

import sympy as sp  # Import SymPy, a computer algebra system written entirely in Python
import os           # Standard Python module for multiplatform OS-level functions
from MoLtimestepping.RK_Butcher_Table_Dictionary import Butcher_dict

# Step 2: Checking if Butcher Table is Diagonal
def diagonal(key):
    diagonal = True #  Start with the Butcher table is diagonal
    Butcher = Butcher_dict[key][0]
    L = len(Butcher)-1 # Establish the number of rows to check for diagonal trait, all bust last row
    row_idx = 0 # Initialize the Butcher table row index
    for i in range(L): # Check all the desired rows
        for j in range(1,row_idx): # Check each element before the diagonal element in a row
            if Butcher[i][j] != sp.sympify(0): # If any non-diagonal coeffcient is non-zero,
                                               # then the table is not diagonal
                diagonal = False
                return diagonal
        row_idx += 1  # Update to check the next row
    return diagonal

# Step 3.a: When allocating memory, we populate a list malloced_gridfunctions,
#         which is used here to determine which gridfunctions need memory freed,
#         via the free() command. Free the mallocs!
def free_allocated_memory(outdir, RK_method, malloced_gridfunctions):
    # This step is made extremely easy, as we had to
    with open(os.path.join(outdir, "RK_Free_Memory.h"), "w") as file:
        file.write("// Code snippet freeing gridfunction memory for \"" + RK_method + "\" method:\n")

        for gridfunction in malloced_gridfunctions:
            file.write("free(" + gridfunction + ");\n")


# # State whether each Butcher table is diagonal or not
# for key, value in Butcher_dict.items():
#     if diagonal(key) == True:
#         print("The RK method "+str(key)+" is diagonal! \n")
#     else:
#         print("The RK method "+str(key)+" is NOT diagonal! \n")
# #################################################################

# Step 3.b: Main driver function for outputting all the MoL C Code
def MoL_C_Code_Generation(RK_method = "RK4", RHS_string = "", post_RHS_string = "",outdir="MoLtimestepping/",
                          MemAllocOnly=False):

    # MoL gridfunctions fall into 3 overlapping categories:
    #           1) y_n=y_i(t_n) gridfunctions y_n_gfs, which stores data for the vector of gridfunctions y_i at t_n,
    #              the start of each MoL timestep.
    #           2) non-y_n gridfunctions, needed to compute the data at t_{n+1}. Often labeled with k_i in the name,
    #              these gridfunctions are *not* needed at the start of each timestep, so are available for temporary
    #              storage when gridfunctions needed for diagnostics are computed at the start of each timestep.
    #              These gridfunctions can also be freed during a regrid, to enable storage for the post-regrid
    #              destination y_n_gfs.
    #           3) Diagnostic output gridfunctions diagnostic_output_gfs, which simply uses the memory from auxiliary
    #              gridfunctions at one auxiliary time to compute diagnostics at t_n.

####### Step 3.b.i: Allocating Memory
    malloc_str = "// Code snippet allocating gridfunction memory for \"" + RK_method + "\" method:\n"

    # Loop over grids
    malloced_gridfunctions = []
    # Set gridfunction type
    type_str = "REAL *restrict "
    # Define a couple useful functions for outputting the needed C code for allocating memory
    def malloc_gfs_str(varname):
        malloced_gridfunctions.append(varname)
        memory_alloc_str = " = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot"+")"
        return type_str + varname + memory_alloc_str + ";\n"
    def diagnostic_output_gfs_equal_to(gfs):
        return type_str + "diagnostic_output_gfs"+" = "+gfs + ";\n"

    # No matter the method we define gridfunctions "y_n_gfs" to store the initial data
    malloc_str += malloc_gfs_str("y_n_gfs")
    if diagonal(RK_method) == True and "RK3" in RK_method:
        malloc_str += malloc_gfs_str("k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs")
        malloc_str += malloc_gfs_str("k2_or_y_nplus_a32_k2_gfs")
        malloc_str += diagnostic_output_gfs_equal_to("k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs")
    else:
        if not diagonal(RK_method):  # Allocate memory for non-diagonal Butcher tables
            # Determine the number of k_i steps based on length of Butcher Table
            num_k = len(Butcher_dict[RK_method][0])-1
            # For non-diagonal tables an intermediate gridfunction "next_y_input" is used for rhs evaluations
            malloc_str += malloc_gfs_str("next_y_input_gfs")
            for i in range(num_k): # Need to allocate all k_i steps for a given method
                malloc_str += malloc_gfs_str("k"+str(i+1)+"_gfs")
            malloc_str += diagnostic_output_gfs_equal_to("k1_gfs")
        else: # Allocate memory for diagonal Butcher tables, which use a "y_nplus1_running_total gridfunction"
            malloc_str += malloc_gfs_str("y_nplus1_running_total_gfs")
            if RK_method != 'Euler': # Allocate memory for diagonal Butcher tables that aren't Euler
                # Need k_odd for k_1,3,5... and k_even for k_2,4,6...
                malloc_str += malloc_gfs_str("k_odd_gfs")
                malloc_str += malloc_gfs_str("k_even_gfs")
            malloc_str += diagnostic_output_gfs_equal_to("y_nplus1_running_total_gfs")
    with open(os.path.join(outdir,"RK_Allocate_Memory.h"), "w") as file:
        file.write(malloc_str)

    if MemAllocOnly:
        free_allocated_memory(outdir,RK_method,malloced_gridfunctions)
        return
########################################################################################################################
# EXAMPLE
# ODE: y' = f(t,y), y(t_0) = y_0
# Starting at time t_n with solution having value y_n and trying to update to y_nplus1 with timestep dt

# Example of scheme for RK4 with k_1, k_2, k_3, k_4 (Using non-diagonal algorithm) Notice this requires storage of
# y_n, y_nplus1, k_1 through k_4

# k_1      = dt*f(t_n, y_n)
# k_2      = dt*f(t_n + 1/2*dt, y_n + 1/2*k_1)
# k_3      = dt*f(t_n + 1/2*dt, y_n + 1/2*k_2)
# k_4      = dt*f(t_n + dt, y_n + k_3)
# y_nplus1 = y_n + 1/3k_1 + 1/6k_2 + 1/6k_3 + 1/3k_4

# Example of scheme RK4 using only k_odd and k_even (Diagonal algroithm) Notice that this only requires storage

# k_odd     = dt*f(t_n, y_n)
# y_nplus1  = 1/3*k_odd
# k_even    = dt*f(t_n + 1/2*dt, y_n + 1/2*k_odd)
# y_nplus1 += 1/6*k_even
# k_odd     = dt*f(t_n + 1/2*dt, y_n + 1/2*k_even)
# y_nplus1 += 1/6*k_odd
# k_even    = dt*f(t_n + dt, y_n + k_odd)
# y_nplus1 += 1/3*k_even
########################################################################################################################

####### Step 3.b.ii: Implementing the Runge Kutta Scheme for Method of Lines Timestepping
    Butcher = Butcher_dict[RK_method][0] # Get the desired Butcher table from the dictionary
    num_steps = len(Butcher)-1 # Specify the number of required steps to update solution
    # Diagonal RK3 only!!!

    def single_RK_substep(commentblock, RHS_str, RHS_input_str, RHS_output_str, RK_lhss_list, RK_rhss_list,
                          post_RHS_list, post_RHS_output_list, indent = "  "):
        return_str  = commentblock + "\n"
        if not isinstance(RK_lhss_list,list):
            RK_lhss_list = [RK_lhss_list]
        if not isinstance(RK_rhss_list,list):
            RK_rhss_list = [RK_rhss_list]

        if not isinstance(post_RHS_list,list):
            post_RHS_list = [post_RHS_list]
        if not isinstance(post_RHS_output_list,list):
            post_RHS_output_list = [post_RHS_output_list]

        # Part 1: RHS evaluation:
        return_str += RHS_str.replace("RK_INPUT_GFS", RHS_input_str).\
                              replace("RK_OUTPUT_GFS",RHS_output_str)+"\n"

        # Part 2: RK update
        return_str += "LOOP_ALL_GFS_GPS"+"(i) {\n"
        for lhs,rhs in zip(RK_lhss_list,RK_rhss_list):
            return_str += indent + lhs + "[i] = " + rhs.replace("_gfs","_gfs") + ";\n"
        return_str += "}\n"

        # Part 3: Call post-RHS functions
        for post_RHS,post_RHS_output in zip(post_RHS_list,post_RHS_output_list):
            return_str += post_RHS.replace("RK_OUTPUT_GFS",post_RHS_output)+"\n"

        return return_str+"\n"

    RK_str = "// C code implementation of " + RK_method + " Method of Lines timestepping.\n"

    # Diagonal RK3 only!!!
    if diagonal(RK_method) and "RK3" in RK_method:
        #  In a diagonal RK3 method, only 3 gridfunctions need be defined. Below implements this approach.

        # k_1
        RK_str += """
// In a diagonal RK3 method like this one, only 3 gridfunctions need be defined. Below implements this approach.
// Using y_n_gfs as input, k1 and apply boundary conditions\n"""

        RK_str += single_RK_substep(
            commentblock = """
// ***k1 substep:***
//  1. We will store k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs now as
//     ...  the update for the next rhs evaluation y_n + a21*k1*dt
// Post-RHS evaluation:
//  1. Apply post-RHS to y_n + a21*k1*dt""",
            RHS_str = RHS_string,
            RHS_input_str = "y_n_gfs", RHS_output_str = "k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs",
            RK_lhss_list = ["k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs"],
            RK_rhss_list = ["("+sp.ccode(Butcher[1][1]).replace("L","")+")*k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i]*dt + y_n_gfs[i]"],
            post_RHS_list = [post_RHS_string], post_RHS_output_list = ["k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs"])

        # k_2
        RK_str += single_RK_substep(
            commentblock="""
// ***k2 substep:***
//    1. Reassign k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs to be the running total y_{n+1}; a32*k2*dt to the running total
//    2. Store k2_or_y_nplus_a32_k2_gfs now as y_n + a32*k2*dt

// Post-RHS evaluation:
//    1. Apply post-RHS to both y_n + a32*k2 (stored in k2_or_y_nplus_a32_k2_gfs)
//       ... and the y_{n+1} running total, as they have not been applied yet to k2-related gridfunctions""",
            RHS_str=RHS_string,
            RHS_input_str="k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs", RHS_output_str="k2_or_y_nplus_a32_k2_gfs",
            RK_lhss_list=["k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs","k2_or_y_nplus_a32_k2_gfs"],
            RK_rhss_list=["("+sp.ccode(Butcher[3][1]).replace("L","")+")*(k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i] - y_n_gfs[i])/("+sp.ccode(Butcher[1][1]).replace("L","")+") + y_n_gfs[i] + ("+sp.ccode(Butcher[3][2]).replace("L","")+")*k2_or_y_nplus_a32_k2_gfs[i]*dt",
                          "("+sp.ccode(Butcher[2][2]).replace("L","")+")*k2_or_y_nplus_a32_k2_gfs[i]*dt + y_n_gfs[i]"],
            post_RHS_list=[post_RHS_string,post_RHS_string],
            post_RHS_output_list=["k2_or_y_nplus_a32_k2_gfs","k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs"])

        # k_3
        RK_str += single_RK_substep(
            commentblock="""
// ***k3 substep:***
//    1. Add k3 to the running total and save to y_n

// Post-RHS evaluation:
//    1. Apply post-RHS to y_n""",
            RHS_str=RHS_string,
            RHS_input_str="k2_or_y_nplus_a32_k2_gfs", RHS_output_str="y_n_gfs",
            RK_lhss_list=["y_n_gfs","k2_or_y_nplus_a32_k2_gfs"],
            RK_rhss_list=["k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i] + ("+sp.ccode(Butcher[3][3]).replace("L","")+")*y_n_gfs[i]*dt"],
            post_RHS_list=[post_RHS_string],
            post_RHS_output_list=["y_n_gfs"])
    else:
        y_n = "y_n_gfs"
        if not diagonal(RK_method):
            for s in range(num_steps):
                next_y_input  = "next_y_input_gfs"

                # If we're on the first step (s=0), we use y_n gridfunction as input.
                #      Otherwise next_y_input is input. Output is just the reverse.
                if s == 0:  # If on first step:
                    RHS_input = y_n
                else:    # If on second step or later:
                    RHS_input = next_y_input
                RHS_output = "k" + str(s + 1) + "_gfs"
                if s == num_steps-1: # If on final step:
                    RK_lhs = y_n
                    RK_rhs = y_n + "[i] + dt*("
                else:                # If on anything but the final step:
                    RK_lhs = next_y_input
                    RK_rhs = y_n + "[i] + dt*("
                for m in range(s+1):
                    if Butcher[s+1][m+1] != 0:
                        if Butcher[s+1][m+1] != 1:
                            RK_rhs += " + k"+str(m+1)+"_gfs[i]*("+sp.ccode(Butcher[s+1][m+1]).replace("L","")+")"
                        else:
                            RK_rhs += " + k"+str(m+1)+"_gfs[i]"
                RK_rhs += " )"

                post_RHS = post_RHS_string
                if s == num_steps-1: # If on final step:
                    post_RHS_output = y_n
                else:                # If on anything but the final step:
                    post_RHS_output = next_y_input

                RK_str += single_RK_substep(
                    commentblock="// ***k" + str(s + 1) + " substep:***",
                    RHS_str=RHS_string,
                    RHS_input_str=RHS_input, RHS_output_str=RHS_output,
                    RK_lhss_list=[RK_lhs],   RK_rhss_list=[RK_rhs],
                    post_RHS_list=[post_RHS],
                    post_RHS_output_list=[post_RHS_output])
        else:  # diagonal case:
            y_nplus1_running_total = "y_nplus1_running_total_gfs"
            if RK_method == 'Euler': # Euler's method doesn't require any k_i, and gets its own unique algorithm
                RK_str += single_RK_substep(
                    commentblock="// ***Euler timestepping only requires one RHS evaluation***",
                    RHS_str=RHS_string,
                    RHS_input_str=y_n, RHS_output_str=y_nplus1_running_total,
                    RK_lhss_list=[y_n],   RK_rhss_list=[y_n+"[i] + "+y_nplus1_running_total+"[i]*dt"],
                    post_RHS_list=[post_RHS_string],
                    post_RHS_output_list=[y_n])
            else:
                for s in range(num_steps):
                    # If we're on the first step (s=0), we use y_n gridfunction as input.
                    # and k_odd as output.
                    if s == 0:
                        RHS_input  = "y_n_gfs"
                        RHS_output = "k_odd_gfs"
                    # For the remaining steps the inputs and ouputs alternate between k_odd and k_even
                    elif s % 2 == 0:
                        RHS_input = "k_even_gfs"
                        RHS_output = "k_odd_gfs"
                    else:
                        RHS_input = "k_odd_gfs"
                        RHS_output = "k_even_gfs"

                    RK_lhs_list = []
                    RK_rhs_list = []
                    if s != num_steps-1:  # For anything besides the final step
                        if s == 0:  # The first RK step
                            RK_lhs_list.append(y_nplus1_running_total)
                            RK_rhs_list.append(RHS_output+"[i]*dt*("+sp.ccode(Butcher[num_steps][s+1]).replace("L","")+")")

                            RK_lhs_list.append(RHS_output)
                            RK_rhs_list.append(y_n+"[i] + "+RHS_output+"[i]*dt*("+sp.ccode(Butcher[s+1][s+1]).replace("L","")+")")
                        else:
                            if Butcher[num_steps][s+1] != 0:
                                RK_lhs_list.append(y_nplus1_running_total)
                                if Butcher[num_steps][s+1] != 1:
                                    RK_rhs_list.append(y_nplus1_running_total+"[i] + "+RHS_output+"[i]*dt*("+sp.ccode(Butcher[num_steps][s+1]).replace("L","")+")")
                                else:
                                    RK_rhs_list.append(y_nplus1_running_total+"[i] + "+RHS_output+"[i]*dt")
                            if Butcher[s+1][s+1] != 0:
                                RK_lhs_list.append(RHS_output)
                                if Butcher[s+1][s+1] != 1:
                                    RK_rhs_list.append(y_n+"[i] + "+RHS_output+"[i]*dt*("+sp.ccode(Butcher[s+1][s+1]).replace("L","")+")")
                                else:
                                    RK_rhs_list.append(y_n+"[i] + "+RHS_output+"[i]*dt")
                        post_RHS_output = RHS_output
                    if s == num_steps-1:  # If on the final step
                        if Butcher[num_steps][s+1] != 0:
                            RK_lhs_list.append(y_n)
                            if Butcher[num_steps][s+1] != 1:
                                RK_rhs_list.append(y_n+"[i] + "+y_nplus1_running_total+"[i] + "+RHS_output+"[i]*dt*("+sp.ccode(Butcher[num_steps][s+1]).replace("L","")+")")
                            else:
                                RK_rhs_list.append(y_n+"[i] + "+y_nplus1_running_total+"[i] + "+RHS_output+"[i]*dt)")
                        post_RHS_output = y_n

                    RK_str += single_RK_substep(
                        commentblock="// ***k" + str(s + 1) + " substep:***",
                        RHS_str=RHS_string,
                        RHS_input_str=RHS_input, RHS_output_str=RHS_output,
                        RK_lhss_list=RK_lhs_list, RK_rhss_list=RK_rhs_list,
                        post_RHS_list=[post_RHS_string],
                        post_RHS_output_list=[post_RHS_output])

    with open(os.path.join(outdir,"RK_MoL.h"), "w") as file:
        file.write(RK_str)

####### Step 3.b.iii: Freeing Allocated Memory
    free_allocated_memory(outdir,RK_method,malloced_gridfunctions)
