from mpmath import mpf, mp, mpc
from UnitTesting.standard_constants import precision

mp.dps = precision
trusted_values_dict = {}

# Generated on: 2019-08-09
trusted_values_dict['ScalarWave_RHSs__ScalarWave_RHSs__globals'] = {'wavespeed': mpf('0.668616460564278702882745619717753'), 'uu_rhs': mpf('0.65261108714780824424650518267299'), 'vv_rhs': mpf('0.368531803639439351490746662269146')}

# Generated on: 2019-08-09
trusted_values_dict['ScalarWaveCurvilinear_RHSs__ScalarWaveCurvilinear_RHSs__globals'] = {'uu_rhs': mpf('0.65261108714780824424650518267299'), 'vv_rhs': mpf('7.91102076566947763349564936041703')}

# Generated on: 2020-05-23
trusted_values_dict['InitialData__InitialData__Type__PlaneWave___globals'] = {'uu_ID': mpf('2.52358067004083477110217386122'), 'vv_ID': mpf('-0.569645247209539485189673979502')}

# Generated on: 2020-07-18
# Notes: Added 1 to uu, so that relative error is well-defined everywhere.
trusted_values_dict['InitialData__InitialData__Type__SphericalGaussian___globals'] = {'uu_ID': mpf('1.365134057130369776890221067476'), 'vv_ID': mpf('0.0563512459589975646763098472872')}
