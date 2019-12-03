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
                break
        row_idx += 1 # Update to check the next row
    return diagonal

# # State whether each Butcher table is diagonal or not
# for key, value in Butcher_dict.items():
#     if diagonal(key) == True:
#         print("The RK method "+str(key)+" is diagonal! \n")
#     else:
#         print("The RK method "+str(key)+" is NOT diagonal! \n")
# #################################################################

# Step 3: Generating the C Code
def MoL_C_Code_Generation(RK_method = "RK4", RHS_string = "", post_RHS_string = "",outdir="MoLtimestepping/"):

####### Step 3a:Allocating Memory
    with open(os.path.join(outdir,"RK_Allocate_Memory.h"), "w") as file:
        file.write("// Code snippet allocating gridfunction memory for \""+str(RK_method)+"\" method:\n")
        # No matter the method we define gridfunctions "y_n_gfs" to store the initial data    
        file.write("REAL *restrict y_n_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);\n")
        if diagonal(RK_method) == True and "RK3" in RK_method:
            file.write("""REAL *restrict k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
REAL *restrict k2_or_y_nplus_a32_k2_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
REAL *restrict diagnostic_output_gfs = k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs;""")
        else:    
            if diagonal(RK_method) == False: #  Allocate memory for non-diagonal Butcher tables 
                # Determine the number of k_i steps based on length of Butcher Table
                num_k = len(Butcher_dict[RK_method][0])-1
                # For non-diagonal tables an intermediate gridfunction "next_y_input" is needed for rhs evaluations
                file.write("REAL *restrict next_y_input_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);\n")
                for i in range(num_k): # Need to allocate all k_i steps for a given method 
                    file.write("REAL *restrict k"+str(i+1)+"_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);\n")
                file.write("REAL *restrict diagnostic_output_gfs = k1_gfs;\n")
            else: # Allocate memory for diagonal Butcher tables, which use a "y_nplus1_running_total gridfunction"
                file.write("REAL *restrict y_nplus1_running_total_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);\n")               
                if RK_method != 'Euler': # Allocate memory for diagonal Butcher tables that aren't Euler
                    # Need k_odd for k_1,3,5... and k_even for k_2,4,6...
                    file.write("REAL *restrict k_odd_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);\n")
                    file.write("REAL *restrict k_even_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);\n")
                file.write("REAL *restrict diagnostic_output_gfs = y_nplus1_running_total_gfs;\n")

######################################################################################################################## 
# EXAMPLE
# ODE: y' = f(t,y), y(t_0) = y_0
# Starting at time t_n with solution having value y_n and trying to update to y_nplus1 with timestep dt

# Example of scheme for RK4 with k_1, k_2, k_3, k_4 (Using non-diagonal algortihm) Notice this requires storage of
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

####### Step 3b: Implementing the Runge Kutta Scheme for Method of Lines Timestepping
    Butcher = Butcher_dict[RK_method][0] # Get the desired Butcher table from the dictionary
    num_steps = len(Butcher)-1 # Specify the number of required steps to update solution
    indent = "  "
    with open(os.path.join(outdir,"RK_MoL.h"), "w") as file:
        file.write("// Code snippet implementing "+RK_method+" algorithm for Method of Lines timestepping\n")
        # Diagonal RK3 only!!!
        if diagonal(RK_method) == True and "RK3" in RK_method:
            #  In a diagonal RK3 method, only 3 gridfunctions need be defined. Below implements this approach.
            file.write("""
// In a diagonal RK3 method like this one, only 3 gridfunctions need be defined. Below implements this approach.
// Using y_n_gfs as input, compute k1 and apply boundary conditions
"""+RHS_string.replace("RK_INPUT_GFS" ,"y_n_gfs").
               replace("RK_OUTPUT_GFS","k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs")+"""
LOOP_ALL_GFS_GPS(i) {
    // Store k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs now as
    //  the update for the next rhs evaluation y_n + a21*k1*dt:
    k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i] = ("""+sp.ccode(Butcher[1][1]).replace("L","")+""")*k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i]*dt + y_n_gfs[i];
}
// Apply boundary conditions to y_n + a21*k1*dt:
"""+post_RHS_string.replace("RK_OUTPUT_GFS","k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs")+"""

// Compute k2 using yn + a21*k1*dt
"""+RHS_string.replace("RK_INPUT_GFS" ,"k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs").
               replace("RK_OUTPUT_GFS","k2_or_y_nplus_a32_k2_gfs")+"""
LOOP_ALL_GFS_GPS(i) {
    // Reassign k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs to be
    //    the running total y_{n+1}
    k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i] = ("""+sp.ccode(Butcher[3][1]).replace("L","")+""")*(k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i] - y_n_gfs[i])/("""+sp.ccode(Butcher[1][1]).replace("L","")+""") + y_n_gfs[i];

    // Add a32*k2*dt to the running total
    k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i]+= ("""+sp.ccode(Butcher[3][2]).replace("L","")+""")*k2_or_y_nplus_a32_k2_gfs[i]*dt;

    // Store k2_or_y_nplus_a32_k2_gfs now as y_n + a32*k2*dt
    k2_or_y_nplus_a32_k2_gfs[i] = ("""+sp.ccode(Butcher[2][2]).replace("L","")+""")*k2_or_y_nplus_a32_k2_gfs[i]*dt + y_n_gfs[i];
}
// Apply boundary conditions to both y_n + a32*k2 (stored in k2_or_y_nplus_a32_k2_gfs)
//    ... and the y_{n+1} running total, as they have not been applied yet to k2-related gridfunctions:
"""+post_RHS_string.replace("RK_OUTPUT_GFS","k2_or_y_nplus_a32_k2_gfs")+"""
"""+post_RHS_string.replace("RK_OUTPUT_GFS","k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs")+"""

// Compute k3
"""+RHS_string.replace("RK_INPUT_GFS" ,"k2_or_y_nplus_a32_k2_gfs").
               replace("RK_OUTPUT_GFS","y_n_gfs")+"""
LOOP_ALL_GFS_GPS(i) {
    // Add k3 to the running total and save to y_n
    y_n_gfs[i] = k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs[i] + ("""+sp.ccode(Butcher[3][3]).replace("L","")+""")*y_n_gfs[i]*dt;
}
// Apply boundary conditions to the running total
"""+post_RHS_string.replace("RK_OUTPUT_GFS","y_n_gfs")+"\n")
        else:    
            y_n           = "y_n_gfs"
            if diagonal(RK_method) == False:
                for s in range(num_steps):
                    next_y_input  = "next_y_input_gfs"

                    # If we're on the first step (s=0), we use y_n gridfunction as input. 
                    #      Otherwise next_y_input is input. Output is just the reverse.
                    if s==0: # If on first step:
                        file.write(RHS_string.replace("RK_INPUT_GFS",y_n).replace("RK_OUTPUT_GFS","k"+str(s+1)+"_gfs")+"\n")
                    else:    # If on second step or later:
                        file.write(RHS_string.replace("RK_INPUT_GFS",next_y_input).replace("RK_OUTPUT_GFS","k"+str(s+1)+"_gfs")+"\n")
                    file.write("LOOP_ALL_GFS_GPS(i) {\n")
                    RK_update_string = ""
                    if s == num_steps-1: # If on final step:
                        RK_update_string += indent + y_n+"[i] += dt*("
                    else:                # If on anything but the final step:
                        RK_update_string += indent + next_y_input+"[i] = "+y_n+"[i] + dt*("
                    for m in range(s+1):
                        if Butcher[s+1][m+1] != 0:
                            if Butcher[s+1][m+1] != 1:
                                RK_update_string += " + k"+str(m+1)+"_gfs[i]*("+sp.ccode(Butcher[s+1][m+1]).replace("L","")+")"
                            else:
                                RK_update_string += " + k"+str(m+1)+"_gfs[i]"
                    RK_update_string += " );\n}\n"
                    file.write(RK_update_string)
                    if s == num_steps-1: # If on final step:
                        file.write(post_RHS_string.replace("RK_OUTPUT_GFS",y_n)+"\n")
                    else:                # If on anything but the final step:
                        file.write(post_RHS_string.replace("RK_OUTPUT_GFS",next_y_input)+"\n")
            else:
                y_nplus1_running_total = "y_nplus1_running_total_gfs"
                if RK_method == 'Euler': # Euler's method doesn't require any k_i, and gets its own unique algorithm
                    file.write(RHS_string.replace("RK_INPUT_GFS",y_n).replace("RK_OUTPUT_GFS",y_nplus1_running_total)+"\n")
                    file.write("LOOP_ALL_GFS_GPS(i) {\n")
                    file.write(indent + y_n+"[i] +=  "+y_nplus1_running_total+"[i]*dt;\n")
                    file.write("}\n")
                    file.write(post_RHS_string.replace("RK_OUTPUT_GFS",y_n)+"\n")
                else:
                    for s in range(num_steps):
                        # If we're on the first step (s=0), we use y_n gridfunction as input. 
                        # and k_odd as output.
                        if s == 0:
                            rhs_input  = "y_n_gfs"
                            rhs_output = "k_odd_gfs"
                        # For the remaining steps the inputs and ouputs alternate between k_odd and k_even
                        elif s%2 == 0:
                            rhs_input = "k_even_gfs"
                            rhs_output = "k_odd_gfs"
                        else:
                            rhs_input = "k_odd_gfs"
                            rhs_output = "k_even_gfs"
                        file.write(RHS_string.replace("RK_INPUT_GFS",rhs_input).replace("RK_OUTPUT_GFS",rhs_output)+"\n")
                        file.write("LOOP_ALL_GFS_GPS(i) {\n")
                        if s == num_steps-1: # If on the final step
                            if Butcher[num_steps][s+1] !=0:
                                if Butcher[num_steps][s+1] !=1:  
                                    file.write(indent+y_n+"[i] += "+y_nplus1_running_total+"[i] + "+rhs_output+"[i]*dt*("+sp.ccode(Butcher[num_steps][s+1]).replace("L","")+");\n")
                                else: 
                                    file.write(indent+y_n+"[i] += "+y_nplus1_running_total+"[i] + "+rhs_output+"[i]*dt;\n")     
                            file.write("}\n")
                            file.write(post_RHS_string.replace("RK_OUTPUT_GFS",y_n)+"\n")
                        else: # For anything besides the final step
                            if s == 0:
                                file.write(indent+y_nplus1_running_total+"[i] = "+rhs_output+"[i]*dt*("+sp.ccode(Butcher[num_steps][s+1]).replace("L","")+");\n")
                                file.write(indent+rhs_output+"[i] = "+y_n+"[i] + "+rhs_output+"[i]*dt*("+sp.ccode(Butcher[s+1][s+1]).replace("L","")+");\n")
                            else:
                                if Butcher[num_steps][s+1] !=0:
                                    if Butcher[num_steps][s+1] !=1:
                                        file.write(indent+y_nplus1_running_total+"[i] += "+rhs_output+"[i]*dt*("+sp.ccode(Butcher[num_steps][s+1]).replace("L","")+");\n")
                                    else: 
                                        file.write(indent+y_nplus1_running_total+"[i] += "+rhs_output+"[i]*dt;\n")
                                if Butcher[s+1][s+1] !=0:
                                    if Butcher[s+1][s+1] !=1:
                                        file.write(indent+rhs_output+"[i] = "+y_n+"[i] + "+rhs_output+"[i]*dt*("+sp.ccode(Butcher[s+1][s+1]).replace("L","")+");\n")
                                    else:
                                        file.write(indent+rhs_output+"[i] = "+y_n+"[i] + "+rhs_output+"[i]*dt;\n")
                            file.write("}\n")
                            file.write(post_RHS_string.replace("RK_OUTPUT_GFS",rhs_output)+"\n")

####### Step 3c: Freeing Allocated Memory
    L = len(Butcher_dict[RK_method][0])-1 # Useful when freeing k_i gridfunctions

    with open(os.path.join(outdir,"RK_Free_Memory.h"), "w") as file:

        file.write("// CODE SNIPPET FOR FREEING ALL ALLOCATED MEMORY FOR "+str(RK_method)+" METHOD:\n")
        if diagonal(RK_method) == True and "RK3" in RK_method:
            file.write("""
free(k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs);
free(k2_or_y_nplus_a32_k2_gfs);
free(y_n_gfs);""")
        else:
            file.write("free(y_n_gfs);\n") 
            if diagonal(RK_method) == False:
                file.write("free(next_y_input_gfs);\n")
                for i in range(L):
                    file.write("free(k"+str(i+1)+"_gfs);\n")
            else:
                file.write("free(y_nplus1_running_total_gfs);\n")
                if RK_method != 'Euler':       
                    file.write("free(k_odd_gfs);\n")
                    file.write("free(k_even_gfs);\n")