// Code snippet implementing RK4 algorithm for Method of Lines timestepping
rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, aux_gfs, y_n_gfs, k_odd_gfs);
LOOP_ALL_GFS_GPS(i) {
  y_nplus1_running_total_gfs[i] = k_odd_gfs[i]*dt*(1.0/6.0);
  k_odd_gfs[i] = y_n_gfs[i] + k_odd_gfs[i]*dt*(1.0/2.0);
}

apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, k_odd_gfs);
enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx,                                         k_odd_gfs);

rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, aux_gfs, k_odd_gfs, k_even_gfs);
LOOP_ALL_GFS_GPS(i) {
  y_nplus1_running_total_gfs[i] += k_even_gfs[i]*dt*(1.0/3.0);
  k_even_gfs[i] = y_n_gfs[i] + k_even_gfs[i]*dt*(1.0/2.0);
}

apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, k_even_gfs);
enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx,                                         k_even_gfs);

rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, aux_gfs, k_even_gfs, k_odd_gfs);
LOOP_ALL_GFS_GPS(i) {
  y_nplus1_running_total_gfs[i] += k_odd_gfs[i]*dt*(1.0/3.0);
  k_odd_gfs[i] = y_n_gfs[i] + k_odd_gfs[i]*dt;
}

apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, k_odd_gfs);
enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx,                                         k_odd_gfs);

rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, aux_gfs, k_odd_gfs, k_even_gfs);
LOOP_ALL_GFS_GPS(i) {
  y_n_gfs[i] += y_nplus1_running_total_gfs[i] + k_even_gfs[i]*dt*(1.0/6.0);
}

apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, y_n_gfs);
enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx,                                         y_n_gfs);

