//          Boundary condtion driver routine: Apply BCs to all
//          boundary faces of the 3D numerical domain, filling in the
//          outer boundary ghost zone layers, starting with the innermost
//          layer and working outward.

#include "sommerfeld_params.h"
#include <string.h>

void apply_bcs_sommerfeld(const paramstruct *restrict params, REAL *restrict xx[3],
                          const bc_struct *restrict bcstruct, const int NUM_GFS,
                          const int8_t *restrict gfs_parity, REAL *restrict gfs,
                          REAL *restrict rhs_gfs)
{

#include "RELATIVE_PATH__set_Cparameters.h" /* Header file containing correct #include for set_Cparameters.h;
                                          * accounting for the relative path */
  if (strcmp(coord, "Cartesian") == 0)
  {
#pragma omp parallel for
    for (int which_gf = 0; which_gf < NUM_GFS; which_gf++)
    {
      REAL var_at_infinity = evolgf_at_inf[which_gf];
      REAL radpower = evolgf_radpower[which_gf];
      REAL char_speed = evolgf_speed[which_gf];

      for (int which_gz = 0; which_gz < NGHOSTS; which_gz++)
      {
        for (int pt = 0; pt < bcstruct->num_ob_gz_pts[which_gz]; pt++)
        {
          int i0 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i0;
          int i1 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i1;
          int i2 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i2;
          int8_t FACEX0 = bcstruct->outer[which_gz][pt].FACEi0;
          int8_t FACEX1 = bcstruct->outer[which_gz][pt].FACEi1;
          int8_t FACEX2 = bcstruct->outer[which_gz][pt].FACEi2;

          REAL xx0 = xx[0][i0];
          REAL xx1 = xx[1][i1];
          REAL xx2 = xx[2][i2];
          REAL dfdx = 0.;
          REAL dfdy = 0.;
          REAL dfdz = 0.;

          // On a +x or -x face, do up/down winding as appropriate:
          if (abs(FACEX0) == 1 || i0 + NGHOSTS >= Nxx_plus_2NGHOSTS0 || i0 - NGHOSTS <= 0)
          {
            int8_t FACE0PARITY = FACEX0;
            if (i0 + NGHOSTS >= Nxx_plus_2NGHOSTS0)
              FACE0PARITY = -1;
            if (i0 - NGHOSTS <= 0)
              FACE0PARITY = +1;
            dfdx = FACE0PARITY * (-3 * gfs[IDX4S(which_gf, i0, i1, i2)] + 4 * gfs[IDX4S(which_gf, i0 + 1 * FACE0PARITY, i1, i2)] - 1 * gfs[IDX4S(which_gf, i0 + 2 * FACE0PARITY, i1, i2)]) * invdx0 * 0.5;

            // Not on a +x or -x face, using centered difference:
          }
          else
          {
            dfdx = (gfs[IDX4S(which_gf, i0 + 1, i1, i2)] - gfs[IDX4S(which_gf, i0 - 1, i1, i2)]) * invdx0 * 0.5;
          }

          // On a +y or -y face, do up/down winding as appropriate:
          if (abs(FACEX1) == 1 || i1 + NGHOSTS >= Nxx_plus_2NGHOSTS1 || i1 - NGHOSTS <= 0)
          {
            int8_t FACE1PARITY = FACEX1;
            if (i1 + NGHOSTS >= Nxx_plus_2NGHOSTS1)
              FACE1PARITY = -1;
            if (i1 - NGHOSTS <= 0)
              FACE1PARITY = +1;
            dfdy = FACE1PARITY * (-3 * gfs[IDX4S(which_gf, i0, i1, i2)] + 4 * gfs[IDX4S(which_gf, i0, i1 + 1 * FACE1PARITY, i2)] - 1 * gfs[IDX4S(which_gf, i0, i1 + 2 * FACE1PARITY, i2)]) * invdx1 * 0.5;

            // Not on a +y or -y face, using centered difference:
          }
          else
          {
            dfdy = (gfs[IDX4S(which_gf, i0, i1 + 1, i2)] - gfs[IDX4S(which_gf, i0, i1 - 1, i2)]) * invdx1 * 0.5;
          }

          // On a +z or -z face, do up/down winding as appropriate:
          if (abs(FACEX2) == 1 || i2 + NGHOSTS >= Nxx_plus_2NGHOSTS2 || i2 - NGHOSTS <= 0)
          {
            int8_t FACE2PARITY = FACEX2;
            if (i2 + NGHOSTS >= Nxx_plus_2NGHOSTS2)
              FACE2PARITY = -1;
            if (i2 - NGHOSTS <= 0)
              FACE2PARITY = +1;
            dfdz = FACE2PARITY * (-3 * gfs[IDX4S(which_gf, i0, i1, i2)] + 4 * gfs[IDX4S(which_gf, i0, i1, i2 + 1 * FACE2PARITY)] - 1 * gfs[IDX4S(which_gf, i0, i1, i2 + 2 * FACE2PARITY)]) * invdx2 * 0.5;

            // Not on a +z or -z face, using centered difference:
          }
          else
          {
            dfdz = (gfs[IDX4S(which_gf, i0, i1, i2 + 1)] - gfs[IDX4S(which_gf, i0, i1, i2 - 1)]) * invdx2 * 0.5;
          }

          REAL invr = 1. / sqrt(xx0 * xx0 + xx1 * xx1 + xx2 * xx2);
          REAL source_rhs = -invr * char_speed * (xx0 * dfdx + xx1 * dfdy + xx2 * dfdz + gfs[IDX4S(which_gf, i0, i1, i2)] - var_at_infinity);
          rhs_gfs[IDX4S(which_gf, i0, i1, i2)] = source_rhs;

          /************* For radial falloff and the extrapolated h'(t) term *************/
          if (radpower > 0)
          {

            // Move one point away from gz point to compare pure advection to df/dt|interior
            int ip0 = i0 + FACEX0;
            int ip1 = i1 + FACEX1;
            int ip2 = i2 + FACEX2;

            REAL xx0 = xx[0][ip0];
            REAL xx1 = xx[1][ip1];
            REAL xx2 = xx[2][ip2];
            REAL dfdx = 0.;
            REAL dfdy = 0.;
            REAL dfdz = 0.;

            // On a +x or -x face, do up/down winding as appropriate:
            if (abs(FACEX0) == 1 || ip0 + NGHOSTS >= Nxx_plus_2NGHOSTS0 || ip0 - NGHOSTS <= 0)
            {
              int8_t FACE0PARITY = FACEX0;
              if (ip0 + NGHOSTS >= Nxx_plus_2NGHOSTS0)
                FACE0PARITY = -1;
              if (ip0 - NGHOSTS <= 0)
                FACE0PARITY = +1;
              dfdx = FACE0PARITY * (-3 * gfs[IDX4S(which_gf, ip0, ip1, ip2)] + 4 * gfs[IDX4S(which_gf, ip0 + 1 * FACE0PARITY, ip1, ip2)] - 1 * gfs[IDX4S(which_gf, ip0 + 2 * FACE0PARITY, ip1, ip2)]) * invdx0 * 0.5;

              // Not on a +x or -x face, using centered difference:
            }
            else
            {
              dfdx = (gfs[IDX4S(which_gf, ip0 + 1, ip1, ip2)] - gfs[IDX4S(which_gf, ip0 - 1, ip1, ip2)]) * invdx0 * 0.5;
            }

            // On a +y or -y face, do up/down winding as appropriate:
            if (abs(FACEX1) == 1 || ip1 + NGHOSTS >= Nxx_plus_2NGHOSTS1 || ip1 - NGHOSTS <= 0)
            {
              int8_t FACE1PARITY = FACEX1;
              if (ip1 + NGHOSTS >= Nxx_plus_2NGHOSTS1)
                FACE1PARITY = -1;
              if (ip1 - NGHOSTS <= 0)
                FACE1PARITY = +1;
              dfdy = FACE1PARITY * (-3 * gfs[IDX4S(which_gf, ip0, ip1, ip2)] + 4 * gfs[IDX4S(which_gf, ip0, ip1 + 1 * FACE1PARITY, ip2)] - 1 * gfs[IDX4S(which_gf, ip0, ip1 + 2 * FACE1PARITY, ip2)]) * invdx1 * 0.5;

              // Not on a +y or -y face, using centered difference:
            }
            else
            {
              dfdy = (gfs[IDX4S(which_gf, ip0, ip1 + 1, ip2)] - gfs[IDX4S(which_gf, ip0, ip1 - 1, ip2)]) * invdx1 * 0.5;
            }

            // On a +z or -z face, do up/down winding as appropriate:
            if (abs(FACEX2) == 1 || ip2 + NGHOSTS >= Nxx_plus_2NGHOSTS2 || ip2 - NGHOSTS <= 0)
            {
              int8_t FACE2PARITY = FACEX2;
              if (ip2 + NGHOSTS >= Nxx_plus_2NGHOSTS2)
                FACE2PARITY = -1;
              if (ip2 - NGHOSTS <= 0)
                FACE2PARITY = +1;
              dfdz = FACE2PARITY * (-3 * gfs[IDX4S(which_gf, ip0, ip1, ip2)] + 4 * gfs[IDX4S(which_gf, ip0, ip1, ip2 + 1 * FACE2PARITY)] - 1 * gfs[IDX4S(which_gf, ip0, ip1, ip2 + 2 * FACE2PARITY)]) * invdx2 * 0.5;

              // Not on a +z or -z face, using centered difference:
            }
            else
            {
              dfdz = (gfs[IDX4S(which_gf, ip0, ip1, ip2 + 1)] - gfs[IDX4S(which_gf, ip0, ip1, ip2 - 1)]) * invdx2 * 0.5;
            }

            REAL rp = sqrt(xx0 * xx0 + xx1 * xx1 + xx2 * xx2);
            REAL invrp = 1. / rp;

            // Pure advection
            REAL extrap_rhs = invrp * char_speed * (xx0 * dfdx + xx1 * dfdy + xx2 * dfdz + gfs[IDX4S(which_gf, ip0, ip1, ip2)] - var_at_infinity);

            // Take difference between pure advection and df/dt|interior
            REAL aux = rhs_gfs[IDX4S(which_gf, ip0, ip1, ip2)] + extrap_rhs;

            // Solve for h'(t)/(r_gz)^n term
            rhs_gfs[IDX4S(which_gf, i0, i1, i2)] += aux * pow(rp * invr, radpower);
          }
        } // END for(int pt=0;pt<num_ob_gz_pts[which_gz];pt++)

        // Then apply INNER (parity) boundary conditions:
        for (int pt = 0; pt < bcstruct->num_ib_gz_pts[which_gz]; pt++)
        {
          const int i0dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i0;
          const int i1dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i1;
          const int i2dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i2;
          const int i0src = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i0;
          const int i1src = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i1;
          const int i2src = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i2;
          const int8_t *prty = bcstruct->inner[which_gz][pt].parity;
          //                printf("%d\n",bcstruct->inner_bc_parity[which_gz][pt].parity[gfs_parity[which_gf]]);
          gfs[IDX4S(which_gf, i0dest, i1dest, i2dest)] =
              bcstruct->inner[which_gz][pt].parity[gfs_parity[which_gf]] * gfs[IDX4S(which_gf, i0src, i1src, i2src)];
        } // END for(int pt=0;pt<num_ib_gz_pts[which_gz];pt++)
      }   // END for(int which_gz = 0; which_gz < NGHOSTS; which_gz++)
    }     // END for(int which_gf=0;which_gf<NUM_GFS;which_gf++)
  }       // END if coord = Cartesian
  else if (strcmp(coord, "Spherical") == 0)
  {
#pragma omp parallel for
    for (int which_gf = 0; which_gf < NUM_GFS; which_gf++)
    {
      REAL var_at_infinity = evolgf_at_inf[which_gf];
      REAL radpower = evolgf_radpower[which_gf];
      REAL char_speed = evolgf_speed[which_gf];

      for (int which_gz = 0; which_gz < NGHOSTS; which_gz++)
      {
        for (int pt = 0; pt < bcstruct->num_ob_gz_pts[which_gz]; pt++)
        {
          int i0 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i0;
          int i1 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i1;
          int i2 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i2;
          // int8_t FACEX0 = bcstruct->outer[which_gz][pt].FACEi0;
          // int8_t FACEX1 = bcstruct->outer[which_gz][pt].FACEi1;
          // int8_t FACEX2 = bcstruct->outer[which_gz][pt].FACEi2;

          REAL invr = 1. / (xx[0][i0]);
          REAL dfdr = 0.;

          // Backwards finite difference stencil
          dfdr = (3 * gfs[IDX4S(which_gf, i0, i1, i2)] - 4 * gfs[IDX4S(which_gf, i0 - 1, i1, i2)] + 1 * gfs[IDX4S(which_gf, i0 - 2, i1, i2)]) * invdx0 * 0.5;

          REAL source_rhs = -char_speed * (dfdr + invr * (gfs[IDX4S(which_gf, i0, i1, i2)] - var_at_infinity));
          rhs_gfs[IDX4S(which_gf, i0, i1, i2)] = source_rhs;

          /////////For radial falloff and the extrapolated h'(t) term////////
          if (radpower > 0)
          {

            int ip0 = i0 - 1; // Move towards the interior of the grid

            REAL invrp = 1. / (xx[0][ip0]);
            REAL dfdr = 0.;

            // Use centered finite difference stencils if possible, otherwise use backwards derivatives
            if (ip0 + NGHOSTS >= Nxx_plus_2NGHOSTS2)
            {
              // Backwards derivative - notice the sign change
              dfdr = (3 * gfs[IDX4S(which_gf, ip0, i1, i2)] - 4 * gfs[IDX4S(which_gf, ip0 - 1, i1, i2)] + 1 * gfs[IDX4S(which_gf, ip0 - 2, i1, i2)]) * invdx0 * 0.5;
            }
            else
            {
              // Standard centered derivative
              dfdr = (gfs[IDX4S(which_gf, ip0 + 1, i1, i2)] - gfs[IDX4S(which_gf, ip0 - 1, i1, i2)]) * invrp * 0.5;
            }

            // Same as before
            REAL extrap_rhs = char_speed * (dfdr + invrp * (gfs[IDX4S(which_gf, ip0, i1, i2)] - var_at_infinity));
            REAL aux = rhs_gfs[IDX4S(which_gf, ip0, i1, i2)] + extrap_rhs;
            rhs_gfs[IDX4S(which_gf, i0, i1, i2)] += aux * pow(xx[0][ip0] * invr, radpower);
          }
        } // END for(int pt=0;pt<num_ob_gz_pts[which_gz];pt++)
          // Then apply INNER (parity) boundary conditions:
        for (int pt = 0; pt < bcstruct->num_ib_gz_pts[which_gz]; pt++)
        {
          const int i0dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i0;
          const int i1dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i1;
          const int i2dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i2;
          const int i0src = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i0;
          const int i1src = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i1;
          const int i2src = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i2;
          const int8_t *prty = bcstruct->inner[which_gz][pt].parity;
          //                printf("%d\n",bcstruct->inner_bc_parity[which_gz][pt].parity[gfs_parity[which_gf]]);
          gfs[IDX4S(which_gf, i0dest, i1dest, i2dest)] =
              bcstruct->inner[which_gz][pt].parity[gfs_parity[which_gf]] * gfs[IDX4S(which_gf, i0src, i1src, i2src)];
        } // END for(int pt=0;pt<num_ib_gz_pts[which_gz];pt++)
      }   // END for(int which_gz = 0; which_gz < NGHOSTS; which_gz++)
    }     // END for(int which_gf=0;which_gf<NUM_GFS;which_gf++)
  }
  else
  {
    printf("ERROR: Sommerfeld boundary conditions are currently only enabled for Cartesian coordinates.\n");
    exit(1);
  } // END coord != Cartesian
} // END function
