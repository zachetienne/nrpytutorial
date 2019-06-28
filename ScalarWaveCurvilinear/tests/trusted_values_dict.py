from mpmath import mpf,mp
from UnitTesting.standard_constants import precision

# Dictionary of trusted values to be used throughout files.
# Standard precision and seed values are precision: 30, seed: 1234.
# Note that changing these may drastically change the calculated values.

mp.dps = precision
trusted_values_dict = dict()

# Generated on: 2019-06-19 15:58:07.085418
trusted_values_dict['ScalarWaveCurvilinear_RHSsGlobals'] = {'uu_rhs': mpf('0.763038382428399942157724162730687'), 'vv_rhs': mpf('29.8343379048485150514807183285427')}