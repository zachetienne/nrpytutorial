""" Unit Testing for Parsing BSSN """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

# pylint: disable = import-error, protected-access, exec-used
from UnitTesting.assert_equal import assert_equal
from nrpylatex import parse_latex
import unittest, sys

import NRPy_param_funcs as par, reference_metric as rfm
import BSSN.BSSN_RHSs as Brhs, BSSN.BSSN_quantities as Bq
import BSSN.BSSN_gauge_RHSs as gaugerhs
import BSSN.BSSN_constraints as bssncon

class TestParser(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None

    @staticmethod
    def test_example_BSSN():
        parse_latex(r"""
            \begin{align}
                % keydef basis [x, y, z]
                % ignore "\\%", "\qquad"

                % vardef -kron 'deltaDD'
                % parse \hat{\gamma}_{ij} = \delta_{ij}
                % assign -diff_type=symbolic -metric 'gammahatDD'
                % vardef -diff_type=dD -symmetry=sym01 'hDD'
                % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
                % assign -diff_type=dD -metric 'gammabarDD'

                % srepl "\beta" -> "\text{vet}"
                % vardef -diff_type=dD 'vetU'
                %% upwind pattern inside Lie derivative expansion
                % srepl -persist "\text{vet}^{<1>} \partial_{<1>}" -> "\text{vet}^{<1>} \vphantom{dupD} \partial_{<1>}"
                %% substitute tensor identity (see appropriate BSSN notebook)
                % srepl "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

                % srepl "\bar{A}" -> "\text{a}"
                % vardef -diff_type=dD -symmetry=sym01 'aDD'
                % assign -metric='gammabar' 'aDD'
                % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
                \partial_t \bar{\gamma}_{ij} &= \mathcal{L}_\beta \bar{\gamma}_{ij}
                    + \frac{2}{3} \bar{\gamma}_{ij} \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{ij} \\

                % srepl "K" -> "\text{trK}"
                % vardef -diff_type=dD 'cf', 'trK'
                %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
                % srepl "e^{-4\phi}" -> "\text{cf}^2"
                % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
                % srepl -persist "\partial_{<1>} \phi" -> "\partial_{<1>} \text{cf} \frac{-1}{2 \text{cf}}"
                % srepl "\partial_<1> \phi" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"
                \partial_t \phi &= \mathcal{L}_\beta \phi + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

                % vardef -diff_type=dD 'alpha'
                % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
                \partial_t K &= \mathcal{L}_\beta K + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{ij} \bar{A}^{ij}
                    - e^{-4\phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi\right) \\

                % srepl "\bar{\Lambda}" -> "\text{lambda}"
                % vardef -diff_type=dD 'lambdaU'
                % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
                % assign -metric='gammabar' 'DeltaUDD'
                % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
                % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
                \partial_t \bar{\Lambda}^i &= \mathcal{L}_\beta \bar{\Lambda}^i + \bar{\gamma}^{jk} \hat{D}_j \hat{D}_k \beta^i
                    + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
                    &\qquad- 2 \bar{A}^{ij} \left(\partial_j \alpha - 6 \alpha \partial_j \phi\right)
                    + 2 \alpha \bar{A}^{jk} \Delta^i_{jk} - \frac{4}{3} \alpha \bar{\gamma}^{ij} \partial_j K \\

                % vardef -diff_type=dD -symmetry=sym01 'RbarDD'
                X_{ij} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi
                    + 2 \bar{D}_i \alpha \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi
                    - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{ij} \\
                \hat{X}_{ij} &= X_{ij} - \frac{1}{3} \bar{\gamma}_{ij} \bar{\gamma}^{kl} X_{kl} \\
                % srepl "\partial_t \text{a}" -> "\text{a_rhs}"
                \partial_t \bar{A}_{ij} &= \mathcal{L}_\beta \bar{A}_{ij} - \frac{2}{3} \bar{A}_{ij} \bar{D}_k \beta^k
                    - 2 \alpha \bar{A}_{ik} \bar{A}^k_j + \alpha \bar{A}_{ij} K + e^{-4\phi} \hat{X}_{ij} \\

                % srepl "\partial_t \alpha" -> "\text{alpha_rhs}"
                \partial_t \alpha &= \mathcal{L}_\beta \alpha - 2 \alpha K \\

                % srepl "B" -> "\text{bet}"
                % vardef -diff_type=dD 'betU'
                % srepl "\partial_t \text{vet}" -> "\text{vet_rhs}"
                \partial_t \beta^i &= \left[\beta^j \vphantom{dupD} \bar{D}_j \beta^i\right] + B^i \\

                % vardef -const 'eta'
                % srepl "\partial_t \text{bet}" -> "\text{bet_rhs}"
                \partial_t B^i &= \left[\beta^j \vphantom{dupD} \bar{D}_j B^i\right]
                    + \frac{3}{4} \left(\partial_t \bar{\Lambda}^i - \left[\beta^j \vphantom{dupD} \bar{D}_j \bar{\Lambda}^i\right]\right) - \eta B^i \\

                % parse \bar{R} = \bar{\gamma}^{ij} \bar{R}_{ij}
                % srepl "\bar{D}^2" -> "\bar{D}^i \bar{D}_i", "\mathcal{<1>}" -> "<1>"
                \mathcal{H} &= \frac{2}{3} K^2 - \bar{A}_{ij} \bar{A}^{ij}
                    + e^{-4\phi} \left(\bar{R} - 8 \bar{D}^i \phi \bar{D}_i \phi - 8 \bar{D}^2 \phi\right) \\

                \mathcal{M}^i &= e^{-4\phi} \left(\bar{D}_j \bar{A}^{ij} + 6 \bar{A}^{ij} \partial_j \phi
                    - \frac{2}{3} \bar{\gamma}^{ij} \partial_j K\right) \\

                \bar{R}_{ij} &= -\frac{1}{2} \bar{\gamma}^{kl} \hat{D}_k \hat{D}_l \bar{\gamma}_{ij}
                    + \frac{1}{2} \left(\bar{\gamma}_{ki} \hat{D}_j \bar{\Lambda}^k + \bar{\gamma}_{kj} \hat{D}_i \bar{\Lambda}^k\right)
                    + \frac{1}{2} \Delta^k \left(\Delta_{ijk} + \Delta_{jik}\right) \\%
                    &\qquad+ \bar{\gamma}^{kl} \left(\Delta^m_{ki} \Delta_{jml} + \Delta^m_{kj} \Delta_{iml} + \Delta^m_{ik} \Delta_{mjl}\right)
            \end{align}
        """, ignore_warning=True)
        par.set_parval_from_str('reference_metric::CoordSystem', 'Cartesian')
        par.set_parval_from_str('BSSN.BSSN_quantities::LeaveRicciSymbolic', 'True')
        rfm.reference_metric()
        Brhs.BSSN_RHSs()
        gaugerhs.BSSN_gauge_RHSs()
        bssncon.BSSN_constraints()
        par.set_parval_from_str('BSSN.BSSN_quantities::LeaveRicciSymbolic', 'False')
        Bq.RicciBar__gammabarDD_dHatD__DGammaUDD__DGammaU()
        assert_equal({'h_rhsDD': h_rhsDD,
                      'cf_rhs': cf_rhs,
                      'trK_rhs': trK_rhs,
                      'Lambdabar_rhsU': Lambdabar_rhsU,
                      'a_rhsDD': a_rhsDD,
                      'alpha_rhs': alpha_rhs,
                      'vet_rhsU': vet_rhsU,
                      'bet_rhsU': bet_rhsU,
                      'H': H,
                      'MU': MU,
                      'RbarDD': RbarDD},
                     {'h_rhsDD': Brhs.h_rhsDD,
                      'cf_rhs': Brhs.cf_rhs,
                      'trK_rhs': Brhs.trK_rhs,
                      'Lambdabar_rhsU': Brhs.Lambdabar_rhsU,
                      'a_rhsDD': Brhs.a_rhsDD,
                      'alpha_rhs': gaugerhs.alpha_rhs,
                      'vet_rhsU': gaugerhs.vet_rhsU,
                      'bet_rhsU': gaugerhs.bet_rhsU,
                      'H': bssncon.H,
                      'MU': bssncon.MU,
                      'RbarDD': Bq.RbarDD},
                    suppress_message=True)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestParser('test_example_BSSN'))
    result = unittest.TextTestRunner().run(suite)
    sys.exit(not result.wasSuccessful())
