#
# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com
import outputC as outC
import os

# Output data on the plane closest to xy or yz
def output_plane_yz_or_xy_body(plane="yz", include_ghosts=True):
    out_str = r"""  // First unpack:
  const paramstruct params = griddata->params;
  REAL *restrict xx[3] = { griddata->xx[0],griddata->xx[1],griddata->xx[2] };
  const int Nx0 = params.Nxx0, Nx1 = params.Nxx1, Nx2 = params.Nxx2;
  const int Nxx_plus_2NGHOSTS0 = params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params.Nxx_plus_2NGHOSTS2;

  // Next set loop bounds, based on coordinate system (params.CoordSystemName)
"""
    if include_ghosts:
        out_str += r"""
  int numpts_i0=Nx0 + 2*NGHOSTS, numpts_i1=Nx1 + 2*NGHOSTS, numpts_i2=Nx2 + 2*NGHOSTS;

  int i0_pts[Nx0 + 2*NGHOSTS], i1_pts[Nx1 + 2*NGHOSTS], i2_pts[Nx2 + 2*NGHOSTS];
#pragma omp parallel for
  for(int i0=0; i0<Nx0 + 2*NGHOSTS; i0++) i0_pts[i0] = i0;
#pragma omp parallel for
  for(int i1=0; i1<Nx1 + 2*NGHOSTS; i1++) i1_pts[i1] = i1;
#pragma omp parallel for
  for(int i2=0; i2<Nx2 + 2*NGHOSTS; i2++) i2_pts[i2] = i2;
"""
    else:
        out_str += r"""
  int numpts_i0=Nx0, numpts_i1=Nx1, numpts_i2=Nx2;

  int i0_pts[Nx0], i1_pts[Nx1], i2_pts[Nx2];
#pragma omp parallel for
  for(int i0=NGHOSTS; i0<Nx0 + NGHOSTS; i0++) i0_pts[i0-NGHOSTS] = i0;
#pragma omp parallel for
  for(int i1=NGHOSTS; i1<Nx1 + NGHOSTS; i1++) i1_pts[i1-NGHOSTS] = i1;
#pragma omp parallel for
  for(int i2=NGHOSTS; i2<Nx2 + NGHOSTS; i2++) i2_pts[i2-NGHOSTS] = i2;
"""
    if plane == "yz":
        out_str += r"""
  if(strstr(params.CoordSystemName, "Cartesian") != NULL) {
    // yz-plane == { x_mid }, where x index is i0
    numpts_i0 = 1;
    i0_pts[0] = (Nx0 + 2*NGHOSTS) / 2;
  } else if(strstr(params.CoordSystemName, "Cylindrical") != NULL) {
    // yz-plane == { phi_min or phi_mid }, where phi index is i1; note that phi_min=-PI and phi_max=+PI, modulo ghostzones
    numpts_i1 = 2;
    // See documentation below.
    i1_pts[0] = (int)(NGHOSTS + 0.25*params.Nxx1 - 0.5);
    i2_pts[1] = (int)(NGHOSTS + 0.75*params.Nxx1 - 0.5);
  } else if(strstr(params.CoordSystemName, "Spherical") != NULL ||
            strstr(params.CoordSystemName, "SymTP")     != NULL) {
    // In Spherical/SymTP coordinates, phi_min=-PI and phi_max=+PI
    //   The yz plane is at -PI/2 and +PI/2.
    //   When Nphi=2, NGHOSTS and NGHOSTS+1 correspond to -PI/2 and +PI/2 *exactly*
    //   WARNING:
    //   When Nphi=4, we don't sample -PI/2 and +PI/2 exactly. Instead we sample:
    //              {-3/4, -1/4, +1/4, +3/4}*PI
    //                      ^           ^  <--- closest planes; at 1 and 3
    //                ^           ^        <--- closest planes; at 0 and 2
    //    The same applies for Nphi=4, 8, 12, etc. Best to choose Nphi a multiple of 2 but not 4
    // General pattern:
    // phi_i = -PI + [(i-NGHOSTS) + (0.5)]*(2PI)/Nphi; // Cell-centered grid.
    /* In isympy:
     i_PIo2, PI, Nphi, NGHOSTS = symbols('i_PIo2 PI Nphi NGHOSTS', real=True)
     expr = -PI + ( (i_PIo2-NGHOSTS) + Rational(1,2) )*(2*PI)/Nphi
     print(solve(expr - PI/2, i_PIo2))
     # Output:
     # [NGHOSTS + 3*Nphi/4 - 1/2]
     print(solve(expr - (-PI/2), i_PIo2))
     # Output:
     # [NGHOSTS + Nphi/4 - 1/2]
     */
    numpts_i2 = 2;
    i2_pts[0] = (int)(NGHOSTS + 0.25*params.Nxx2 - 0.5);
    i2_pts[1] = (int)(NGHOSTS + 0.75*params.Nxx2 - 0.5);
  } else {
    fprintf(stderr, "output_"""+plane+r"""...() ERROR: params.CoordSystemName == %s unsupported.\n", params.CoordSystemName);
    exit(1);
  }
"""
    elif plane == "xy":
        out_str += r"""
  if(strstr(params.CoordSystemName, "Cartesian") != NULL) {
    // xy-plane == { z_mid }, where z index is i2
    numpts_i2 = 1;
    i2_pts[0] = (Nx2 + 2*NGHOSTS) / 2;
  } else if(strstr(params.CoordSystemName, "Cylindrical") != NULL) {
    // xy-plane == { z_mid }, where z index is i2
    numpts_i2 = 1;
    i2_pts[0] = (Nx2 + 2*NGHOSTS) / 2;
  } else if(strstr(params.CoordSystemName, "Spherical") != NULL ||
            strstr(params.CoordSystemName, "SymTP")     != NULL) {
    // xy-plane == { theta_mid }, where theta index is i1
    numpts_i1 = 1;
    i1_pts[0] = (Nx1 + 2*NGHOSTS) / 2;;
  } else {
    fprintf(stderr, "output_"""+plane+r"""...() ERROR: params.CoordSystemName == %s unsupported.\n", params.CoordSystemName);
    exit(1);
  }
"""
    return out_str


def add_to_Cfunction_dict__plane_diagnostics(plane="xy", include_ghosts=False,
                                               list_of_outputs=None, num_sig_figs=8):
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = """Output various quantities at points closest to
the given coordinate system's """+plane+""" plane.
"""
    c_type = "void"
    name = plane + "_plane_diagnostics"
    params = """FILE *fp, const griddata_struct *restrict griddata,
    const REAL *restrict y_n_gfs, const REAL *restrict diagnostic_output_gfs"""
    body  = output_plane_yz_or_xy_body(plane=plane, include_ghosts=include_ghosts)
    printf_char_string = "\"%e %e %e"
    output_quantities = ",xCart[0],xCart[1],xCart[2]"
    for output in list_of_outputs:
        printf_char_string += " %e"
        output_quantities += ",\n           " + output
    printf_char_string += "\\n\""
    printf_char_string = printf_char_string.replace("%e", "%."+str(num_sig_figs)+"e")
    body += r"""
  LOOP_NOOMP(i0_pt,0,numpts_i0, i1_pt,0,numpts_i1, i2_pt,0,numpts_i2) {
    const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
    REAL xCart[3];
    xx_to_Cart(&params, xx, i0,i1,i2, xCart);
    int idx = IDX3S(i0,i1,i2);
    fprintf(fp, """+printf_char_string+output_quantities+r""");
  }
"""
    outC.add_to_Cfunction_dict(
        includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        rel_path_to_Cparams=os.path.join("."), enableCparameters=False)
