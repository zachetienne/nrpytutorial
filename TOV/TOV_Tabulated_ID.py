from sys import argv, exit
from numpy import log10, fabs
from TOV.TOV_Solver_tabulated import TOV_Solver
from TOV.tabulated_eos_fcts import read_and_slice_EOS_table_at_given_temperature
# from TOV.Polytropic_EOSs import set_up_EOS_parameters__Read_et_al_input_variables
from time import time

def solve_tov_equations(tov_mass, tov_temperature, eos_table_path=None, tol=1e-7, outfile="./"):
  eos_type="tabulated"

#   outfile = "outputTOV_@name@.txt"

#   if "ls220" in eos_table_path.lower():
#     outfile = outfile.replace("@name@", "ls220")
#     print(f"Auto-detected LS220 EOS from table name. Output file: {outfile}")
#   elif "sfho" in eos_table_path.lower():
#     outfile = outfile.replace("@name@", "sfho")
#     print(f"Auto-detected SFHo EOS from table name. Output file: {outfile}")
#   elif "sly4" in eos_table_path.lower():
#     outfile = outfile.replace("@name@", "sly4")
#     print(f"Auto-detected SLy4 EOS from table name. Output file: {outfile}")
#   else:
#     outfile = outfile.replace("@name@", "tabulated")
#     print(f"Could not auto-detect EOS name. Output file: {outfile}")

  eos = read_and_slice_EOS_table_at_given_temperature(eos_table_path, tov_temperature)
  # eos = set_up_EOS_parameters__Read_et_al_input_variables(eos_table_path, units="G=c=Msun=1")

  # Compute M0
  rho0 = 8e-4
  while True:
    M0, _, _ = TOV_Solver(eos,
                          eos_type=eos_type,
                          outfile=outfile,
                          rho_baryon_central=rho0,
                          verbose = False,
                          return_M_RSchw_and_Riso=True,
                          accuracy="medium",
                          integrator_type="default")

    if M0 < tov_mass:
      break
    else:
      rho0 /= 1.1

  # Compute M1
  rho1 = 1.2e-3
  while True:
    M1, _, _ = TOV_Solver(eos,
                          eos_type=eos_type,
                          outfile=outfile,
                          rho_baryon_central=rho1,
                          verbose = False,
                          return_M_RSchw_and_Riso=True,
                          accuracy="medium",
                          integrator_type="default")
    if M1 > tov_mass:
      break
    else:
      rho1 *= 1.1

  if fabs(M1 - tov_mass) > fabs(M0 - tov_mass):
    rho0, rho1 = rho1, rho0
    M0  , M1   = M1  , M0

  # Now start bisection:
  print(".----.----------.----------.")
  print("| It |    dM    |   drho   |")
  print(".----.----------.----------.")
  counter = 0
  while True:
    counter += 1
    rhoc    = 10**( 0.5*( log10(rho0) + log10(rho1) ) )
    M, R_Schw_TOV, R_iso_TOV = TOV_Solver(eos,
                        eos_type=eos_type,
                        outfile=outfile,
                        rho_baryon_central=rhoc,
                        verbose = False,
                        return_M_RSchw_and_Riso=True,
                        accuracy="medium",
                        integrator_type="default")

    dM = fabs(M-tov_mass)
    if dM < fabs(M1-tov_mass):
      rho0, rho1 = rho1, rhoc
      M0  , M1   = M1  , M
    else:
      rho0, M0 = rhoc, M

    drho = fabs(rho1 - rho0)
    print(f"| {counter:02d} | {dM:8.2e} | {drho:8.2e} |")
    if dM < tol or drho < 1e-15:
      break

  print(".----.----------.----------.")
  print(f"Converged! central density = {rho1:.6f}, TOV mass: {M1:.6f}, TOV radius: {R_Schw_TOV:.6f}")

  return rho1, M1, R_Schw_TOV

# if __name__ == '__main__':

#   if len(argv) != 4:
#     print(f"Usage: python {__file__.split('/')[-1]} <target mass> <temperature> <eos table>")
#     exit(1)

#   tov_mass        = float(argv[1])
#   tov_temperature = float(argv[2])
#   eos_table_path  = argv[3]

#   start = time()
#   solve_tov_equations(tov_mass, tov_temperature, eos_table_path=eos_table_path)
#   print(f"Found TOV solution in {time()-start:.2f} seconds")