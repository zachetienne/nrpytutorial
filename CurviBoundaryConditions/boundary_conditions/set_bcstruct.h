
// set_bcstruct() loops from the innermost boundary
//      ghostzones on the cube ("which_gz==0",
//      corresponding to the single layer of ghostzones
//      closest to the interior data), and at each
//      ghostzone layer, we apply the following 5-step
//      algorithm:
// Step 1: Count the number of outer and inner
//         boundary points, store to
//         num_ob_pts and num_ib_pts, respectively.
// Step 2: Now that we know the number of outer
//         boundary points on this ghostzone layer,
//         allocate memory needed for storing the
//         outer and inner boundary condition data.
// Step 2.a: At all outer boundary ghost zones, allocate
//           memory for a single member of the outer_bc
//           data type.
// Step 2.b: At all inner boundary ghost zones, allocate
//           memory for a single member of the inner_bc
//           data type.
// Step 3: Store the number of outer and inner boundary
//         points on each ghostzone layer, where e.g.,
//         which_gz==0 corresponds to the innermost
//         ghostzones on the numerical domain.
// Step 4: Store information needed for outer boundary
//         conditions, to outer_bc_dest_pt and
//         outer_bc_face arrays.
// Step 5: Store information needed for inner boundary
//         conditions, including interior point to which
//         inner ghost zone maps, and parity conditions
//         for all 10 gridfunction types.
void set_bcstruct(const paramstruct *restrict params,
                  gz_map *restrict bc_gz_map,
                  parity_condition *bc_parity_conditions,
                  bc_struct *restrict bcstruct) {
#include "RELATIVE_PATH__set_Cparameters.h" /* Header file containing correct #include for set_Cparameters.h;
                                             * accounting for the relative path */
    int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
    int imax[3] = { Nxx_plus_2NGHOSTS0-NGHOSTS, Nxx_plus_2NGHOSTS1-NGHOSTS, Nxx_plus_2NGHOSTS2-NGHOSTS };

    // Loop from the innermost ghostzone on the cube (which_gz==0) and work outward.
    //      This ordering is necessary, as ghostzones at which_gz==1 will generally
    //      depend on ghostzones at which_gz==0 being already set.
    for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {

    // Step 1: Count the number of outer and inner
    //         boundary points, store to
    //         num_ob_pts and num_ib_pts, respectively.
#define COUNT_INNER_OR_OUTER if(bc_gz_map[IDX3S(i0,i1,i2)].i0==-1) { num_ob_pts++;} else { num_ib_pts++; }
        int num_ob_pts = 0;
        int num_ib_pts = 0;
        LOOP_REGION(imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2]) { COUNT_INNER_OR_OUTER } imin[0]--;
        LOOP_REGION(imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2]) { COUNT_INNER_OR_OUTER } imax[0]++;
        LOOP_REGION(imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2]) { COUNT_INNER_OR_OUTER } imin[1]--;
        LOOP_REGION(imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2]) { COUNT_INNER_OR_OUTER } imax[1]++;
        LOOP_REGION(imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2]) { COUNT_INNER_OR_OUTER } imin[2]--;
        LOOP_REGION(imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1) { COUNT_INNER_OR_OUTER } imax[2]++;

        // Step 2: Now that we know the number of outer boundary points on this ghostzone
        //    layer, we allocate memory needed for storing the outer and inner boundary
        //     condition data.

        // Step 2.a: At all outer boundary ghost zones, allocate memory for a single member of the outer_bc
        //           data type.
        bcstruct->outer[which_gz] = (outer_bc *)malloc(sizeof(outer_bc)*num_ob_pts);
        // Step 2.b: At all inner boundary ghost zones, allocate memory for a single member of the inner_bc
        //           data type.
        bcstruct->inner[which_gz] = (inner_bc *)malloc(sizeof(inner_bc)*num_ib_pts);

        // Step 3: Store the number of outer and inner boundary points on each ghostzone layer, where e.g.,
        //         which_gz==0 corresponds to the innermost ghostzones on the numerical domain.
        bcstruct->num_ob_gz_pts[which_gz] = num_ob_pts;
        bcstruct->num_ib_gz_pts[which_gz] = num_ib_pts;

        // Reset imin[] and imax[], to prepare for the next step.
        for(int ii=0;ii<3;ii++) {imin[ii]++; imax[ii]--;}

        // Step 4: Store information needed for outer boundary conditions, to outer_bc_dest_pt[which_gz][]
        //         and outer_bc_face[which_gz][] arrays:
#define OB_SET(facei0,facei1,facei2) if(bc_gz_map[IDX3S(i0,i1,i2)].i0==-1) {  \
    bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i0 = i0;     \
    bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i1 = i1;     \
    bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i2 = i2;     \
    bcstruct->outer[which_gz][pt].FACEi0= facei0; \
    bcstruct->outer[which_gz][pt].FACEi1= facei1; \
    bcstruct->outer[which_gz][pt].FACEi2= facei2; \
    pt++; }

        int pt = 0;
        LOOP_REGION(imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2]) {OB_SET(MINFACE,NUL,NUL)} imin[0]--;
        LOOP_REGION(imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2]) {OB_SET(MAXFACE,NUL,NUL)} imax[0]++;
        LOOP_REGION(imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2]) {OB_SET(NUL,MINFACE,NUL)} imin[1]--;
        LOOP_REGION(imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2]) {OB_SET(NUL,MAXFACE,NUL)} imax[1]++;
        LOOP_REGION(imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2]) {OB_SET(NUL,NUL,MINFACE)} imin[2]--;
        LOOP_REGION(imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1) {OB_SET(NUL,NUL,MAXFACE)} imax[2]++;
        // fprintf(stderr,"num OB points with which_gz = %d: %d | should be: %d\n",which_gz,pt,num_ob_gz_pts[which_gz]);

        // Reset imin[] and imax[], to prepare for the next step.
        for(int ii=0;ii<3;ii++) {imin[ii]++; imax[ii]--;}

        // Step 5: Store information needed for inner boundary conditions, including interior point to which
        //         inner ghost zone maps, and parity conditions for all 10 gridfunction types.
#define IB_SET if(bc_gz_map[IDX3S(i0,i1,i2)].i0!=-1) { \
    bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i0=i0; \
    bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i1=i1; \
    bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i2=i2; \
    bcstruct->inner[which_gz][pt].inner_bc_src_pt.i0 =bc_gz_map[IDX3S(i0,i1,i2)].i0; \
    bcstruct->inner[which_gz][pt].inner_bc_src_pt.i1 =bc_gz_map[IDX3S(i0,i1,i2)].i1; \
    bcstruct->inner[which_gz][pt].inner_bc_src_pt.i2 =bc_gz_map[IDX3S(i0,i1,i2)].i2; \
    for(int ii=0;ii<10;ii++) { \
      bcstruct->inner[which_gz][pt].parity[ii] = \
                              (int8_t)bc_parity_conditions[IDX3S(i0,i1,i2)].parity[ii]; } \
    pt++; }

        pt = 0;
        LOOP_REGION(imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2]) {IB_SET} imin[0]--;
        LOOP_REGION(imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2]) {IB_SET} imax[0]++;
        LOOP_REGION(imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2]) {IB_SET} imin[1]--;
        LOOP_REGION(imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2]) {IB_SET} imax[1]++;
        LOOP_REGION(imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2]) {IB_SET} imin[2]--;
        LOOP_REGION(imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1) {IB_SET} imax[2]++;

    } // END for(int which_gz = 0; which_gz < NGHOSTS; which_gz++)
} // END function
