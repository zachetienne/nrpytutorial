import unittest
# First we import needed core NRPy+ modules
import NRPy_param_funcs as par
import reference_metric as rfm
import grid as gri

import BSSN.BrillLindquist as bl
import BSSN.ADM_Exact_Spherical_or_Cartesian_to_BSSNCurvilinear as AtoB

import hashlib
import sys

class TestStringMethods(unittest.TestCase):
    # Set spatial dimension (must be 3 for BSSN)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)

    # Then we set the coordinate system for the numerical grid
    par.set_parval_from_str("reference_metric::CoordSystem","Spherical")
    rfm.reference_metric() # Create ReU, ReDD needed for rescaling B-L initial data, generating BSSN RHSs, etc.
    
    bl.BrillLindquist(ComputeADMGlobalsOnly = True)

    def test_BL_ID_ADM(self):
        everything = bl.alphaCart
        for i in range(3):
            everything += bl.betaCartU[i] + bl.BCartU[i]
            for j in range(3):
                everything += bl.gammaCartDD[i][j] + bl.KCartDD[i][j]
        md5sum = "empty"
        if sys.version_info[0]==2:
            md5sum = hashlib.md5(str(everything)).hexdigest()
        elif sys.version_info[0]==3:
            md5sum = hashlib.md5(str(everything).encode('utf-8')).hexdigest()

        #print(md5sum)
        self.assertEqual(md5sum, 'd1d9fb9bd0c0ce61c06e6e2f6e386040')

    def test_BL_ID(self):
        cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU = \
                                                 AtoB.Convert_Spherical_or_Cartesian_ADM_to_BSSN_curvilinear("Cartesian", bl.Cartxyz, bl.gammaCartDD, 
                                                                                                             bl.KCartDD, bl.alphaCart, bl.betaCartU, bl.BCartU)
        everything = cf+alpha+trK
        for i in range(3):
            everything += lambdaU[i]+vetU[i]+betU[i]
            for j in range(i,3):
                everything += hDD[i][j] + aDD[i][j]
        md5sum = "empty"
        if sys.version_info[0]==2:
            md5sum = hashlib.md5(str(everything)).hexdigest()
        elif sys.version_info[0]==3:
            md5sum = hashlib.md5(str(everything).encode('utf-8')).hexdigest()

        #print(md5sum)
        self.assertEqual(md5sum, 'eb6f859a0016fc1679c45af662ea32ef')
        # OLD ONE: self.assertEqual(md5sum, 'cb8bc25380b8fd5bd0e2e10579259ac9')

if __name__ == '__main__':
    unittest.main()
    
