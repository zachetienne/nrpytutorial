from mpmath import mpf,mp
from UnitTesting.standard_constants import precision

# Dictionary of trusted values to be used throughout files.
# Standard precision and seed values are precision: 30, seed: 1234.
# Note that changing these may drastically change the calculated values.

mp.dps = precision

trusted_values_dict = dict()

# Generated on: 2019-06-20 13:12:44.704524
trusted_values_dict['GiRaFFEfood_HOGlobals'] = {'AD[0]': mpf('-0.94115744530567775184433060969805'), 'AD[1]': mpf('1.69321905451215397667294884427835'), 'AD[2]': mpf('0.0'), 'ValenciavU[0]': mpf('0.294521454719030203077048644263695'), 'ValenciavU[1]': mpf('0.163706555966521359089651866391384'), 'ValenciavU[2]': mpf('-0.975070156720603120338029391237849')}
