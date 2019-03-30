import unittest
import sympy as sp
import NRPy_param_funcs as par
import BSSN.BSSN_RHSs_new as rhs
import BSSN.BSSN_gauge_RHSs as gaugerhs
import hashlib
import sys

def get_md5sum(sympy_expr):
    if sys.version_info[0]==2:
        return hashlib.md5(str(sympy_expr)).hexdigest()
    elif sys.version_info[0]==3:
        return hashlib.md5(str(sympy_expr).encode('utf-8')).hexdigest()
    sys.exit(1)

class TestStringMethods(unittest.TestCase):
    par.set_parval_from_str("BSSN.BSSN_gauge_RHSs::ShiftEvolutionOption", "GammaDriving2ndOrder_Covariant")
    rhs.BSSN_RHSs()
    gaugerhs.BSSN_gauge_RHSs()

    def test_BSSN_RHSs_scalars(self):
        everyscalar = gaugerhs.alpha_rhs + rhs.cf_rhs + rhs.trK_rhs
        md5sum = get_md5sum(everyscalar)
#        print(md5sum)
        self.assertEqual(md5sum, 'b7b76505a53a38a507f45417857ebb7c')

    def test_BSSN_RHSs_vectors(self):
        everyvector = sp.sympify(0)
        for i in range(3):
            everyvector += gaugerhs.bet_rhsU[i] + gaugerhs.vet_rhsU[i] + rhs.lambda_rhsU[i]
        md5sum = get_md5sum(everyvector)
#        print("vector: ",md5sum)
        self.assertEqual(md5sum, '24b47058967ebe5501790f3d813d6c84')

    def test_BSSN_RHSs_tensors(self):
        everytensor = sp.sympify(0)
        for i in range(3):
            for j in range(i,3):
                everytensor += rhs.a_rhsDD[i][j] + rhs.h_rhsDD[i][j]
        md5sum = get_md5sum(everytensor)
#        print("tensor: ",md5sum)
        self.assertEqual(md5sum, 'd84fc94358305b7135dc18680089dff9')

if __name__ == '__main__':
    unittest.main()
    
