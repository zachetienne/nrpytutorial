from mpmath import mpf,mp
from UnitTesting.standard_constants import precision

# Dictionary of trusted values to be used throughout files.
# Standard precision and seed values are precision: 30, seed: 1234.
# Note that changing these may drastically change the calculated values.

mp.dps = precision
trusted_values_dict = dict()

# Generated on: 2019-06-19 15:55:32.335570
trusted_values_dict['InitialData_PlaneWaveGlobals'] = {'uu_ID': mpf('2.18311491393033683821392419836456'), 'vv_ID': mpf('-0.952771904111410843220904372551045')}
# Generated on: 2019-06-19 15:55:32.358916
trusted_values_dict['ScalarWave_RHSsGlobals'] = {'wavespeed': mpf('0.969158912337787880328707563437431'), 'uu_rhs': mpf('0.954450607652970251254358778361204'), 'vv_rhs': mpf('1.62823587362693526874974669214808')}