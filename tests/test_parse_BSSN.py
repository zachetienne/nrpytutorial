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
            % coord [x, y, z]
            % ignore "\begin{align}" "\end{align}" "\\%" "\qquad"

            \begin{align}
                % define gammahatDD --dim 3 --zeros
                % \hat{\gamma}_{ii} = 1 % noimpsum
                % assign gammahatDD --metric
                % define hDD --dim 3 --suffix dD --sym sym01
                % \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
                % assign gammabarDD --suffix dD --metric

                % srepl "\beta" -> "\mathrm{vet}"
                % define vetU --dim 3 --suffix dD
                %% upwind pattern inside Lie derivative expansion
                % srepl "\mathrm{vet}^{<1..>} \partial_{<1..>}" -> "\mathrm{vet}^{<1..>} % suffix dupD
                \partial_{<1..>}" --persist
                %% substitute tensor identity (see appropriate BSSN notebook)
                % srepl "\bar{D}_k \mathrm{vet}^k" -> "(\partial_k \mathrm{vet}^k + \frac{\partial_k \mathrm{gammahatdet} \mathrm{vet}^k}{2 \mathrm{gammahatdet}})"

                % srepl "\bar{A}" -> "\mathrm{a}"
                % define aDD --dim 3 --suffix dD --sym sym01
                % assign aDD --metric gammabar
                % srepl "\partial_t \bar{\gamma}" -> "\mathrm{h_rhs}"
                \partial_t \bar{\gamma}_{ij} &= \mathcal{L}_\beta \bar{\gamma}_{ij} + \frac{2}{3} \bar{\gamma}_{ij} \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{ij} \\

                % define cf trK --dim 3 --suffix dD
                % srepl "K" -> "\mathrm{trK}"
                %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
                % srepl "e^{-4\phi}" -> "\mathrm{cf}^2"
                % srepl "\partial_t \phi = <1..> \\" -> "\mathrm{cf_rhs} = -2 \mathrm{cf} (<1..>) \\"
                % srepl "\partial_{<1..>} \phi" -> "\partial_{<1..>} \mathrm{cf} \frac{-1}{2 \mathrm{cf}}" --persist
                % srepl "\partial_<1> \phi" -> "\partial_<1> \mathrm{cf} \frac{-1}{2 \mathrm{cf}}"
                \partial_t \phi &= \mathcal{L}_\beta \phi + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

                % define alpha --dim 3 --suffix dD
                % srepl "\partial_t \mathrm{trK}" -> "\mathrm{trK_rhs}"
                \partial_t K &= \mathcal{L}_\beta K + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{ij} \bar{A}^{ij}
                    - e^{-4\phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi\right) \\

                % srepl "\bar{\Lambda}" -> "\mathrm{lambda}"
                % define lambdaU --dim 3 --suffix dD
                % \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
                % assign DeltaUDD --metric gammabar
                % \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
                % srepl "\partial_t \mathrm{lambda}" -> "\mathrm{Lambdabar_rhs}"
                \partial_t \bar{\Lambda}^i &= \mathcal{L}_\beta \bar{\Lambda}^i + \bar{\gamma}^{jk} \hat{D}_j \hat{D}_k \beta^i
                    + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
                    &\qquad- 2 \bar{A}^{ij} \left(\partial_j \alpha - 6 \alpha \partial_j \phi\right)
                    + 2 \alpha \bar{A}^{jk} \Delta^i_{jk} - \frac{4}{3} \alpha \bar{\gamma}^{ij} \partial_j K \\

                % define RbarDD --dim 3 --suffix dD --sym sym01
                X_{ij} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi
                    + 2 \bar{D}_i \alpha \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi
                    - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{ij} \\
                \hat{X}_{ij} &= X_{ij} - \frac{1}{3} \bar{\gamma}_{ij} \bar{\gamma}^{kl} X_{kl} \\
                % srepl "\partial_t \mathrm{a}" -> "\mathrm{a_rhs}"
                \partial_t \bar{A}_{ij} &= \mathcal{L}_\beta \bar{A}_{ij} - \frac{2}{3} \bar{A}_{ij} \bar{D}_k \beta^k
                    - 2 \alpha \bar{A}_{ik} \bar{A}^k_j + \alpha \bar{A}_{ij} K + e^{-4\phi} \hat{X}_{ij} \\

                % srepl "\partial_t \alpha" -> "\mathrm{alpha_rhs}"
                \partial_t \alpha &= \mathcal{L}_\beta \alpha - 2 \alpha K \\

                % srepl "B" -> "\mathrm{bet}"
                % define betU --dim 3 --suffix dD
                % srepl "\partial_t \mathrm{vet}" -> "\mathrm{vet_rhs}"
                \partial_t \beta^i &= \left[\beta^j % suffix dupD
                \bar{D}_j \beta^i\right] + B^i \\

                % define eta --const
                % srepl "\partial_t \mathrm{bet}" -> "\mathrm{bet_rhs}"
                \partial_t B^i &= \left[\beta^j % suffix dupD
                \bar{D}_j B^i\right]
                    + \frac{3}{4} \left(\partial_t \bar{\Lambda}^i - \left[\beta^j % suffix dupD
                    \bar{D}_j \bar{\Lambda}^i\right]\right) - \eta B^i \\

                % \bar{R} = \bar{\gamma}^{ij} \bar{R}_{ij}
                % srepl "\bar{D}^2" -> "\bar{D}^i \bar{D}_i"
                % srepl "\mathcal{<1>}" -> "<1>"
                \mathcal{H} &= \frac{2}{3} K^2 - \bar{A}_{ij} \bar{A}^{ij}
                    + e^{-4\phi} \left(\bar{R} - 8 \bar{D}^i \phi \bar{D}_i \phi - 8 \bar{D}^2 \phi\right) \\

                \mathcal{M}^i &= e^{-4\phi} \left(\bar{D}_j \bar{A}^{ij} + 6 \bar{A}^{ij} \partial_j \phi
                    - \frac{2}{3} \bar{\gamma}^{ij} \partial_j K\right) \\

                \bar{R}_{ij} &= -\frac{1}{2} \bar{\gamma}^{kl} \hat{D}_k \hat{D}_l \bar{\gamma}_{ij}
                    + \frac{1}{2} \left(\bar{\gamma}_{ki} \hat{D}_j \bar{\Lambda}^k + \bar{\gamma}_{kj} \hat{D}_i \bar{\Lambda}^k\right)
                    + \frac{1}{2} \Delta^k \left(\Delta_{ijk} + \Delta_{jik}\right) \\%
                    &\qquad+ \bar{\gamma}^{kl} \left(\Delta^m_{ki} \Delta_{jml} + \Delta^m_{kj} \Delta_{iml} + \Delta^m_{ik} \Delta_{mjl}\right)
            \end{align}
        """)
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
                    suppress_message=False)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestParser('test_example_BSSN'))
    result = unittest.TextTestRunner().run(suite)
    sys.exit(not result.wasSuccessful())
