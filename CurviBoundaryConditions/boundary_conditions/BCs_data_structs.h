
typedef struct __ghostzone_map__ {
  short i0,i1,i2; // i0,i1,i2 stores values from -1 (used to indicate outer boundary)
                  // to Nxx_plus_2NGHOSTS*. We assume that grid extents beyond the
                  // limits of short (i.e., beyond about 32,000) are unlikely. This
                  // can be easily extended if needed, though.
} gz_map;

const int8_t MAXFACE = -1;
const int8_t NUL     = +0;
const int8_t MINFACE = +1;

typedef struct __parity__ {
  int8_t parity[10]; // We store the 10 parity conditions in 10 int8_t integers,
                     // one for each condition. Note that these conditions can
                     // only take one of two values: +1 or -1, hence the use of
                     // int8_t, the smallest C data type.
} parity_condition;

typedef struct __inner_bc__ {
    gz_map inner_bc_dest_pt;
    gz_map inner_bc_src_pt;
    int8_t parity[10]; // We store the 10 parity conditions in 10 int8_t integers,
                       // one for each condition. Note that these conditions can
                       // only take one of two values: +1 or -1, hence the use of
                       // int8_t, the smallest C data type.
} inner_bc;

typedef struct __outer_bc__ {
    gz_map outer_bc_dest_pt;
    int8_t FACEi0,FACEi1,FACEi2; // FACEi* takes values of -1, 0, and +1 only,
                                 // corresponding to MAXFACE, NUL, and MINFACE
                                 // respectively.
                                 // Thus int8_t (one byte each, the smallest C
                                 // type) is sufficient.
} outer_bc;

typedef struct __bcstruct__ {
    outer_bc **outer; // Array of 1D arrays, of length
                      //   [NGHOSTS][num_ob_gz_pts[which_outer_ghostzone_point]]

    inner_bc **inner; // Array of 1D arrays, of length
                      //   [NGHOSTS][num_ib_gz_pts[which_inner_ghostzone_point]]

    // Arrays storing number of outer/inner boundary ghostzone points at each ghostzone,
    //   of length NGHOSTS:
    int     *num_ob_gz_pts;
    int     *num_ib_gz_pts;
} bc_struct;
