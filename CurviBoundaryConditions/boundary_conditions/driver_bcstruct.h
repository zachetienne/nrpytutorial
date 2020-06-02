
// Step 1: Allocate memory storage for bc_gz_map, which
//         in the case a boundary point is a *parity*
//         boundary, is set to the interior, non-
//         boundary point corresponding to the same
//         Cartesian gridpoint. Otherwise bc_gz_map
//         is set to (i0,i1,i2) = (-1,-1,-1).
gz_map *bc_gz_map = (gz_map *)malloc(sizeof(gz_map)*Nxx_plus_2NGHOSTS_tot);

// Step 2: Allocate memory storage for bc_parity_conditions,
//         which store parity conditions for all 10
//         gridfunction types at all grid points.
parity_condition *bc_parity_conditions = (parity_condition *)malloc(sizeof(parity_condition)*Nxx_plus_2NGHOSTS_tot);

// Step 3: Set bc_gz_map and bc_parity_conditions at *all*
//         points; on the boundary and otherwise.
set_up__bc_gz_map_and_parity_condns(&params, xx, bc_gz_map,
                                    bc_parity_conditions
                                   );

// Step 4: Declare and allocate memory for bcstruct,
//         which will store all information needed for
//         applying the boundary conditions.
bcstruct.outer = (outer_bc **)malloc(sizeof(outer_bc *)*NGHOSTS);
bcstruct.inner = (inner_bc **)malloc(sizeof(inner_bc *)*NGHOSTS);
bcstruct.num_ob_gz_pts = (    int *)malloc(sizeof(int)*NGHOSTS);
bcstruct.num_ib_gz_pts = (    int *)malloc(sizeof(int)*NGHOSTS);

// Step 4: Store all information needed to quickly and
//         efficiently apply boundary conditions. This
//         function transfers all information from
//         bc_gz_map (defined at *all gridpoints*) into
//         bcstruct (defined only at boundary points).
//         Thus when this function has finished,
//         bc_gz_map is no longer needed.
set_bcstruct(&params,bc_gz_map,
             bc_parity_conditions,
             &bcstruct);

// Step 5: As described in Step 4, bc_gz_map is no
//         longer needed at this point, so we free its
//         memory. Farewell, friend!
free(bc_gz_map);
free(bc_parity_conditions);
