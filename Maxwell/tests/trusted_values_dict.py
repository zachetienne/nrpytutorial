from mpmath import mpf,mp
from UnitTesting.standard_constants import precision

# Dictionary of trusted values to be used throughout files.
# Standard precision and seed values are precision: 30, seed: 1234.
# Note that changing these may drastically change the calculated values.

mp.dps = precision
trusted_values_dict = dict()

# Generated on: 2019-06-24 12:08:40.504117
trusted_values_dict['MaxwellCartesian_EvolGlobals'] = {'ArhsD[0]': mpf('-0.950915101018727879507020176139165'), 'ArhsD[1]': mpf('-1.10339392620588513150233683326111'), 'ArhsD[2]': mpf('-0.769335868654741065500887471807761'), 'ErhsD[0]': mpf('1.71087347555010495881703991789342'), 'ErhsD[1]': mpf('-2.75117124738722802644005371722078'), 'ErhsD[2]': mpf('1.1111903263907933510550915059574'), 'psi_rhs': mpf('-0.232755662973305208943682954171212'), 'Gamma_rhs': mpf('-0.0942848938593562703491925868017232'), 'Cviola': mpf('2.8019895707527610013891389176519')}

# Generated on: 2019-06-24 12:08:40.589134
trusted_values_dict['MaxwellCartesian_IDGlobals'] = {'AidD[0]': mpf('0.0'), 'AidD[1]': mpf('0.0'), 'AidD[2]': mpf('0.0'), 'EidD[0]': mpf('-1.20505683837101607802241844855587'), 'EidD[1]': mpf('0.444218089464664718238492841056562'), 'EidD[2]': mpf('-1.33160376093988613357396497777671'), 'psi_ID': mpf('0.0'), 'Gamma_ID': mpf('0.0')}
