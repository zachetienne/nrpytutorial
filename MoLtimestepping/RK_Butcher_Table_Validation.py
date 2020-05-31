# As documented in the NRPy+ tutorial module
#   Tutorial-RK_Butcher_Table_Validation.ipynb ,
#   this module will validate the Python dictionary
#   of Butcher tables.

# Authors: Brandon Clark
#          Zachariah B. Etienne
#          zachetie **at** gmail **dot* com


# Step 1: Initialize needed Python/NRPy+ modules
import sympy as sp              # SymPy: The Python computer algebra package upon which NRPy+ depends
import numpy as np              # NumPy: A numerical methods module for Python

# Step 2a:  Defining the right-hand side of the ODE
rhs_dict = {}
def fypt(y,t):
    return y+t
rhs_dict['ypt'] = fypt

def fy(y,t):
    return y
rhs_dict['y'] = fy

def feypt(y,t):
    return sp.exp(1.0*(y+t))
rhs_dict['eypt'] = feypt

def ftpoly6(y,t):
    return 2*t**6-389*t**5+15*t**4-22*t**3+81*t**2-t+42
rhs_dict['tpoly6'] = ftpoly6

def ftpoly5(y,t):
    return t**5 + t**4 + t**3 + t**2 + t + 1
rhs_dict["tpoly5"] = ftpoly5

def fty2pt(y,t):
    return t*y**2+t
rhs_dict["ty2pt"] = fty2pt

def fymcost(y,t):
    return y-sp.cos(t)
rhs_dict["ymcost"] = fymcost

def fypsint(y,t):
    return y+sp.sin(t)
rhs_dict['ypsint'] = fypsint




# Step 2b: Defining a Validation Function

def Validate(Butcher_dict, Butcher_key, yn, tn, rhs_key):
    # Set needed symbolic expressions
    t, dt = sp.symbols('t dt')

    # 1. First we solve the ODE exactly
    y = sp.Function('y')
    sol = sp.dsolve(sp.Eq(y(t).diff(t), rhs_dict[rhs_key](y(t), t)), y(t)).rhs
    constants = sp.solve([sol.subs(t,tn)-yn])
    exact = sol.subs(constants)

    # 2. Now we solve the ODE numerically using specified Butcher table

    # Access the requested Butcher table
    Butcher = Butcher_dict[Butcher_key][0]
    # Determine number of predictor-corrector steps
    L = len(Butcher)-1
    # Set a temporary array for update values
    k = np.zeros(L, dtype=object)
    # Initialize intermediate variable
    yhat = 0
    # Initialize the updated solution
    ynp1 = 0
    for i in range(L):
        #Initialize and approximate update for solution
        yhat = yn
        for j in range(i):
            # Update yhat for solution using a_ij Butcher table coefficients
            yhat += Butcher[i][j+1]*k[j]
            if Butcher_key == "DP8" or Butcher_key == "L6":
                yhat = 1.0*sp.N(yhat,20) # Otherwise the adding of fractions kills performance.
        # Determine the next corrector variable k_i using c_i Butcher table coefficients
        k[i] = dt*rhs_dict[rhs_key](yhat, tn + Butcher[i][0]*dt)
        # Update the solution at the next iteration ynp1 using Butcher table coefficients
        ynp1 += Butcher[L][i+1]*k[i]
    # Finish determining the solution for the next iteration
    ynp1 += yn

    # Determine the order of the RK method
    order = Butcher_dict[Butcher_key][1]+2
    # Produces Taylor series of exact solution at t=tn about t = 0 with the specified order
    exact_series = sp.series(exact.subs(t, dt),dt, 0, order)
    num_series = sp.series(ynp1, dt, 0, order)
    diff = exact_series-num_series
    return diff

# # Step 2d: Validating Convergence with ScalarWave PDE in Cartesian Coordinates
# def fd_order(RK_order):
#     if (RK_order+1)%2==0: # If RK_order is odd, then set FD_order (must be even) to RK_order+1
#         return RK_order+1
#     else:
#         return RK_order+2 # If RK_order is even, then set FD_order (must be even) to RK_order+2

# # CFL_FACTOR = 0.5
# for key,value in Butcher_dict.items():
#     for CFL_FACTOR in [0.5]:
#         with open("ScalarWave/params.txt", "w") as file:
#             # DON'T ADD SPACES TO THE FOLLOWING!
#             file.write(
# """RK_method="""+str(key)+"""
# FD_order="""+str(fd_order(value[1]))+"""
# REAL=double
# CFL_FACTOR="""+str(CFL_FACTOR))
#         outfilename = "Out"+str(key)+"_CFL_"+str(CFL_FACTOR)+".png"
#         os.system("rm -f "+outfilename)
#         os.system("jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=-1 Tutorial-Start_to_Finish-ScalarWave.ipynb")
#         if os.path.isfile(outfilename):
#             print("Successfully generated figure "+outfilename)
#             #display(Image(filename=outfilename,width=600))
#         else:
#             print("ERROR: Failed to generate figure "+outfilename)
#             sys.exit(1)
