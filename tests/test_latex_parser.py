""" Unit Testing for LaTeX Parser """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

# pylint: disable=import-error
import sys; sys.path.append('..')
from latex_parser import OverrideWarning, parse_expr, parse
from warnings import filterwarnings
import unittest, sys

class TestParser(unittest.TestCase):
    """ Unit Testing for LaTeX Parser """

    def setUp(self):
        self.maxDiff = None
        filterwarnings('ignore', category=OverrideWarning)

    def test_expression_1(self):
        expr = r'-(\frac{2}{3} + 2\sqrt[5]{x + 3})'
        self.assertEqual(
            str(parse_expr(expr)),
            '-2*(x + 3)**(1/5) - 2/3'
        )

    def test_expression_2(self):
        expr = r'e^{\ln x} + \sin(\sin^{-1} y) - \tanh(xy)'
        self.assertEqual(
            str(parse_expr(expr)),
            'x + y - tanh(x*y)'
        )

    def test_example_1(self):
        self.assertEqual(
            parse(r"""
                % hUD [4]: nosym;
                h = h^\mu{}_\mu
            """),
            ['hUD', 'h']
        )
        self.assertEqual(str(h),
            'hUD00 + hUD11 + hUD22 + hUD33'
        )

    def test_example_2(self):
        self.assertEqual(
            parse(r"""
                % gUU [3]: metric, vD [3]: nosym;
                v^\mu = g^{\mu\nu} v_\nu
            """),
            ['gDD', 'gdet', 'gUU', 'vD', 'vU']
        )
        self.assertEqual(str(vU),
            '[gUU00*vD0 + gUU01*vD1 + gUU02*vD2, gUU01*vD0 + gUU11*vD1 + gUU12*vD2, gUU02*vD0 + gUU12*vD1 + gUU22*vD2]'
        )

    def test_example_3(self):
        self.assertEqual(
            parse(r"""
                % gDD [2]: metric, RUU [2]: sym01;
                \begin{align*}
                    R &= g_{ab} R^{ab} \\
                    G^{ab} &= R^{ab} - \frac{1}{2}g^{ab}R
                \end{align*}
            """),
            ['gUU', 'gdet', 'gDD', 'RUU', 'R', 'GUU']
        )
        self.assertEqual(str(R),
            'RUU00*gDD00 + 2*RUU01*gDD01 + RUU11*gDD11'
        )
        self.assertEqual(str(GUU),
            '[[RUU00 + gDD11*(-RUU00*gDD00 - 2*RUU01*gDD01 - RUU11*gDD11)/(2*(gDD00*gDD11 - gDD01**2)), RUU01 - gDD01*(-RUU00*gDD00 - 2*RUU01*gDD01 - RUU11*gDD11)/(2*(gDD00*gDD11 - gDD01**2))], [RUU01 - gDD01*(-RUU00*gDD00 - 2*RUU01*gDD01 - RUU11*gDD11)/(2*(gDD00*gDD11 - gDD01**2)), RUU11 + gDD00*(-RUU00*gDD00 - 2*RUU01*gDD01 - RUU11*gDD11)/(2*(gDD00*gDD11 - gDD01**2))]]'
        )

    def test_example_4(self):
        self.assertEqual(
            parse(r"""
                % epsilonDDD [3]: permutation, vU [3]: nosym, wU [3]: nosym;
                u_i = \epsilon_{ijk} v^j w^k
            """),
            ['epsilonDDD', 'vU', 'wU', 'uD']
        )
        self.assertEqual(str(uD),
            '[vU1*wU2 - vU2*wU1, -vU0*wU2 + vU2*wU0, vU0*wU1 - vU1*wU0]'
        )

    def test_example_5(self):
        self.assertEqual(
            parse(r"""
                % FUU [4]: anti01;
                J^\mu = (4\pi k)^{-1} \nabla_\nu F^{\mu\nu}
            """, debug=True),
            ['FUU', 'FUU_dD', 'gDD', 'gdet', 'gUU', 'gDD_dD', 'GammaUDD', 'FUU_cdD', 'JU']
        )
        self.assertEqual(GammaUDD,
            '[[[sum([gUU[mu][c]*gDD_dD[a][c][b]/2 for c in range(4)]) - sum([gUU[mu][c]*gDD_dD[b][a][c]/2 for c in range(4)]) + sum([gUU[mu][c]*gDD_dD[c][b][a]/2 for c in range(4)]) for a in range(4)] for b in range(4)] for mu in range(4)]'
        )
        self.assertEqual(FUU_cdD,
            '[[[sum([FUU[b][nu]*GammaUDD[mu][b][a] for b in range(4)]) + sum([FUU[mu][b]*GammaUDD[nu][b][a] for b in range(4)]) + FUU_dD[mu][nu][a] for a in range(4)] for mu in range(4)] for nu in range(4)]'
        )
        self.assertEqual(JU,
            '[sum([FUU_cdD[mu][nu][nu]/(4*pi*k) for nu in range(4)]) for mu in range(4)]'
        )

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestParser('test_expression_1'))
    suite.addTest(TestParser('test_expression_2'))
    suite.addTest(TestParser('test_example_1'))
    suite.addTest(TestParser('test_example_2'))
    suite.addTest(TestParser('test_example_3'))
    suite.addTest(TestParser('test_example_4'))
    suite.addTest(TestParser('test_example_5'))
    result = unittest.TextTestRunner().run(suite)
    sys.exit(not result.wasSuccessful())
