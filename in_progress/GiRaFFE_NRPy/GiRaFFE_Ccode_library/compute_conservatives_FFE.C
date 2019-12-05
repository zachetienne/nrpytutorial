void GiRaFFE_NRPy_compute_conservatives(const REAL gxxL,const REAL gxyL,const REAL gxzL,const REAL gyyL,const REAL gyzL,const REAL gzzL,
                                      const REAL BxL, const REAL ByL, const REAL BzL, const REAL vxL, const REAL vyL, const REAL vzL,
                                      //const REAL betaxL, const REAL betayL, const REAL betazL, const REAL alpL,
                                      const REAL sqrtg,REAL *StildeD0L, REAL *StildeD1L, REAL *StildeD2L) {

  //const REAL fourpialpha_inv = 1.0/( 4.0*M_PI*(METRIC[LAPM1] + 1.0) );
  const REAL fourpi_inv = 1.0/( 4.0*M_PI );

  const REAL B2 = gxxL*BxL*BxL + gyyL*ByL*ByL + gzzL*BzL*BzL
    + 2.0*(gxyL*BxL*ByL + gxzL*BxL*BzL + gyzL*ByL*BzL);


  // NOTE: SIGNIFICANTLY MODIFIED FROM ILLINOISGRMHD VERSION:
  //       velocities in GiRaFFE are defined to be "drift" velocity.
  //       cf. Eqs 47 and 85 in http://arxiv.org/pdf/1310.3274.pdf 
  // Modified again from the original GiRaFFE to use Valencia velocity

  const REAL v_xL = gxxL*vxL + gxyL*vyL + gxzL*vzL;
  const REAL v_yL = gxyL*vxL + gyyL*vyL + gyzL*vzL;
  const REAL v_zL = gxzL*vxL + gyzL*vyL + gzzL*vzL;
  
  /*
   * Comments:
   * Eq. 85 in https://arxiv.org/pdf/1310.3274.pdf:
   * v^i = 4 pi alpha * (gamma^{ij} tilde{S}_j) / (sqrtgamma * B^2) - beta^i
   * which implies that
   * (v^i + beta^i)*(sqrtgamma * B^2)/(4 pi alpha) = gamma^{ij} tilde{S}_j
   * Multiply both sides by gamma_{ik}:
   * gamma_{ik} (v^i + beta^i)*(sqrtgamma * B^2)/(4 pi alpha) = gamma_{ik} gamma^{ij} tilde{S}_j
   * 
   * -> tilde{S}_k = gamma_{ik} (v^i + beta^i)*(sqrtgamma * B^2)/(4 pi alpha)
   */

  *StildeD0L = v_xL * sqrtg * B2 * fourpi_inv;
  *StildeD1L = v_yL * sqrtg * B2 * fourpi_inv;
  *StildeD2L = v_zL * sqrtg * B2 * fourpi_inv;
}
