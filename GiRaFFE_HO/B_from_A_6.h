{
   /* 
    * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
    */
   const double gammaDD00 = aux_gfs[IDX4(GAMMADD00GF, i0,i1,i2)];
   const double gammaDD01 = aux_gfs[IDX4(GAMMADD01GF, i0,i1,i2)];
   const double gammaDD02 = aux_gfs[IDX4(GAMMADD02GF, i0,i1,i2)];
   const double gammaDD11 = aux_gfs[IDX4(GAMMADD11GF, i0,i1,i2)];
   const double gammaDD12 = aux_gfs[IDX4(GAMMADD12GF, i0,i1,i2)];
   const double gammaDD22 = aux_gfs[IDX4(GAMMADD22GF, i0,i1,i2)];
   const double AD0_i0_i1_i2m3 = in_gfs[IDX4(AD0GF, i0,i1,i2-3)];
   const double AD0_i0_i1_i2m2 = in_gfs[IDX4(AD0GF, i0,i1,i2-2)];
   const double AD0_i0_i1_i2m1 = in_gfs[IDX4(AD0GF, i0,i1,i2-1)];
   const double AD0_i0_i1m3_i2 = in_gfs[IDX4(AD0GF, i0,i1-3,i2)];
   const double AD0_i0_i1m2_i2 = in_gfs[IDX4(AD0GF, i0,i1-2,i2)];
   const double AD0_i0_i1m1_i2 = in_gfs[IDX4(AD0GF, i0,i1-1,i2)];
   const double AD0_i0_i1p1_i2 = in_gfs[IDX4(AD0GF, i0,i1+1,i2)];
   const double AD0_i0_i1p2_i2 = in_gfs[IDX4(AD0GF, i0,i1+2,i2)];
   const double AD0_i0_i1p3_i2 = in_gfs[IDX4(AD0GF, i0,i1+3,i2)];
   const double AD0_i0_i1_i2p1 = in_gfs[IDX4(AD0GF, i0,i1,i2+1)];
   const double AD0_i0_i1_i2p2 = in_gfs[IDX4(AD0GF, i0,i1,i2+2)];
   const double AD0_i0_i1_i2p3 = in_gfs[IDX4(AD0GF, i0,i1,i2+3)];
   const double AD1_i0_i1_i2m3 = in_gfs[IDX4(AD1GF, i0,i1,i2-3)];
   const double AD1_i0_i1_i2m2 = in_gfs[IDX4(AD1GF, i0,i1,i2-2)];
   const double AD1_i0_i1_i2m1 = in_gfs[IDX4(AD1GF, i0,i1,i2-1)];
   const double AD1_i0m3_i1_i2 = in_gfs[IDX4(AD1GF, i0-3,i1,i2)];
   const double AD1_i0m2_i1_i2 = in_gfs[IDX4(AD1GF, i0-2,i1,i2)];
   const double AD1_i0m1_i1_i2 = in_gfs[IDX4(AD1GF, i0-1,i1,i2)];
   const double AD1_i0p1_i1_i2 = in_gfs[IDX4(AD1GF, i0+1,i1,i2)];
   const double AD1_i0p2_i1_i2 = in_gfs[IDX4(AD1GF, i0+2,i1,i2)];
   const double AD1_i0p3_i1_i2 = in_gfs[IDX4(AD1GF, i0+3,i1,i2)];
   const double AD1_i0_i1_i2p1 = in_gfs[IDX4(AD1GF, i0,i1,i2+1)];
   const double AD1_i0_i1_i2p2 = in_gfs[IDX4(AD1GF, i0,i1,i2+2)];
   const double AD1_i0_i1_i2p3 = in_gfs[IDX4(AD1GF, i0,i1,i2+3)];
   const double AD2_i0_i1m3_i2 = in_gfs[IDX4(AD2GF, i0,i1-3,i2)];
   const double AD2_i0_i1m2_i2 = in_gfs[IDX4(AD2GF, i0,i1-2,i2)];
   const double AD2_i0_i1m1_i2 = in_gfs[IDX4(AD2GF, i0,i1-1,i2)];
   const double AD2_i0m3_i1_i2 = in_gfs[IDX4(AD2GF, i0-3,i1,i2)];
   const double AD2_i0m2_i1_i2 = in_gfs[IDX4(AD2GF, i0-2,i1,i2)];
   const double AD2_i0m1_i1_i2 = in_gfs[IDX4(AD2GF, i0-1,i1,i2)];
   const double AD2_i0p1_i1_i2 = in_gfs[IDX4(AD2GF, i0+1,i1,i2)];
   const double AD2_i0p2_i1_i2 = in_gfs[IDX4(AD2GF, i0+2,i1,i2)];
   const double AD2_i0p3_i1_i2 = in_gfs[IDX4(AD2GF, i0+3,i1,i2)];
   const double AD2_i0_i1p1_i2 = in_gfs[IDX4(AD2GF, i0,i1+1,i2)];
   const double AD2_i0_i1p2_i2 = in_gfs[IDX4(AD2GF, i0,i1+2,i2)];
   const double AD2_i0_i1p3_i2 = in_gfs[IDX4(AD2GF, i0,i1+3,i2)];
   const double AD_dD01 = invdx1*(-3.0/4.0*AD0_i0_i1m1_i2 + (3.0/20.0)*AD0_i0_i1m2_i2 - 1.0/60.0*AD0_i0_i1m3_i2 + (3.0/4.0)*AD0_i0_i1p1_i2 - 3.0/20.0*AD0_i0_i1p2_i2 + (1.0/60.0)*AD0_i0_i1p3_i2);
   const double AD_dD02 = invdx2*(-3.0/4.0*AD0_i0_i1_i2m1 + (3.0/20.0)*AD0_i0_i1_i2m2 - 1.0/60.0*AD0_i0_i1_i2m3 + (3.0/4.0)*AD0_i0_i1_i2p1 - 3.0/20.0*AD0_i0_i1_i2p2 + (1.0/60.0)*AD0_i0_i1_i2p3);
   const double AD_dD10 = invdx0*(-3.0/4.0*AD1_i0m1_i1_i2 + (3.0/20.0)*AD1_i0m2_i1_i2 - 1.0/60.0*AD1_i0m3_i1_i2 + (3.0/4.0)*AD1_i0p1_i1_i2 - 3.0/20.0*AD1_i0p2_i1_i2 + (1.0/60.0)*AD1_i0p3_i1_i2);
   const double AD_dD12 = invdx2*(-3.0/4.0*AD1_i0_i1_i2m1 + (3.0/20.0)*AD1_i0_i1_i2m2 - 1.0/60.0*AD1_i0_i1_i2m3 + (3.0/4.0)*AD1_i0_i1_i2p1 - 3.0/20.0*AD1_i0_i1_i2p2 + (1.0/60.0)*AD1_i0_i1_i2p3);
   const double AD_dD20 = invdx0*(-3.0/4.0*AD2_i0m1_i1_i2 + (3.0/20.0)*AD2_i0m2_i1_i2 - 1.0/60.0*AD2_i0m3_i1_i2 + (3.0/4.0)*AD2_i0p1_i1_i2 - 3.0/20.0*AD2_i0p2_i1_i2 + (1.0/60.0)*AD2_i0p3_i1_i2);
   const double AD_dD21 = invdx1*(-3.0/4.0*AD2_i0_i1m1_i2 + (3.0/20.0)*AD2_i0_i1m2_i2 - 1.0/60.0*AD2_i0_i1m3_i2 + (3.0/4.0)*AD2_i0_i1p1_i2 - 3.0/20.0*AD2_i0_i1p2_i2 + (1.0/60.0)*AD2_i0_i1p3_i2);
   /* 
    * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
    */
   const double tmp0 = pow(gammaDD00*gammaDD11*gammaDD22 - gammaDD00*pow(gammaDD12, 2) - pow(gammaDD01, 2)*gammaDD22 + 2*gammaDD01*gammaDD02*gammaDD12 - pow(gammaDD02, 2)*gammaDD11, -1.0/2.0);
   aux_gfs[IDX4(BU0GF, i0, i1, i2)] = -AD_dD12*tmp0 + AD_dD21*tmp0;
   aux_gfs[IDX4(BU1GF, i0, i1, i2)] = AD_dD02*tmp0 - AD_dD20*tmp0;
   aux_gfs[IDX4(BU2GF, i0, i1, i2)] = -AD_dD01*tmp0 + AD_dD10*tmp0;
}
