""" Unit Testing for LaTeX Parser """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

# pylint: disable=import-error,protected-access
# import sys; sys.path.append('..')
from latex_parser import Tensor, OverrideWarning
from latex_parser import Parser, parse_expr, parse
from sympy import Function, Symbol
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

    def test_expression_3(self):
        tensor = Tensor(Function('Tensor')(Symbol('TUU'), Symbol('mu'), Symbol('nu')), 4)
        self.assertEqual(
            Parser._generate_covdrv(tensor, [('beta', 'D')]),
            r'\nabla_\beta T^{\mu \nu} = \partial_\beta (T^{\mu \nu}) + \Gamma^\mu_{a \beta} (T^{a \nu}) + \Gamma^\nu_{a \beta} (T^{\mu a})'
        )
        tensor = Tensor(Function('Tensor')(Symbol('TUD'), Symbol('mu'), Symbol('nu')), 4)
        self.assertEqual(
            Parser._generate_covdrv(tensor, [('beta', 'D')]),
            r'\nabla_\beta T^\mu_\nu = \partial_\beta (T^\mu_\nu) + \Gamma^\mu_{a \beta} (T^a_\nu) - \Gamma^a_{\nu \beta} (T^\mu_a)'
        )
        tensor = Tensor(Function('Tensor')(Symbol('TDD'), Symbol('mu'), Symbol('nu')), 4)
        self.assertEqual(
            Parser._generate_covdrv(tensor, [('beta', 'D')]),
            r'\nabla_\beta T_{\mu \nu} = \partial_\beta (T_{\mu \nu}) - \Gamma^a_{\mu \beta} (T_{a \nu}) - \Gamma^a_{\nu \beta} (T_{\mu a})'
        )

    def test_expression_4(self):
        tensor = Tensor(Function('Tensor')(Symbol('vU'), Symbol('mu')), 4)
        self.assertEqual(
            Parser._generate_covdrv(tensor, [('a', 'U'), ('b', 'D')]),
            r'\nabla_a \nabla_b v^\mu = \partial_a (\partial_b (v^\mu) + \Gamma^\mu_{c b} (v^c)) + \Gamma^\mu_{c a} (\partial_b (v^c) + \Gamma^c_{d b} (v^d))'
        )

    def test_expression_5(self):
        self.assertEqual(
            parse(r"""
                % define vU (4);
                T^{\mu\nu} = \nabla^\nu v^\mu
            """),
            ['vU', 'gDD', 'gdet', 'gUU', 'vU_dD', 'gDD_dD', 'GammaUDD', 'vU_cdD', 'vU_cdU', 'TUU']
        )

    def test_assignment_1(self):
        self.assertEqual(
            parse(r"""
                % define vU (2), wU (2);
                T^{ij}_k = (v^i w^j)_{,k}
            """),
            ['vU', 'wU', 'vU_dD', 'wU_dD', 'TUUD']
        )
        self.assertEqual(str(TUUD[0][0][0]),
            'vU0*wU_dD00 + vU_dD00*wU0'
        )

    def test_assignment_2(self):
        self.assertEqual(
            parse(r"""
                % define vU (2), const w;
                T^i_k = (v^i w)_{,k}
            """),
            ['vU', 'w', 'vU_dD', 'TUD']
        )
        self.assertEqual(str(TUD),
            '[[vU_dD00*w, vU_dD10*w], [vU_dD01*w, vU_dD11*w]]'
        )

    def test_example_1(self):
        self.assertEqual(
            parse(r"""
                % define nosym hUD (4);
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
                % define metric gUU (3), vD (3);
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
                % define metric gDD (2), sym01 RUU (2);
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
                % define permutation epsilonDDD (3);
                % define vU (3), wU (3);
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
                % define metric gUU (4), kronecker deltaUU (4);
                % parse g^{\mu\nu} = \delta^{\mu\nu};
                \begin{align*}
                    g^{0 0} &= -\left(1 - \frac{\mathop{r_s}}{r}\right) \\
                    g^{1 1} &=  \left(1 - \frac{\mathop{r_s}}{r}\right)^{-1} \\
                    g^{2 2} &= r^{{2}} \\
                    g^{3 3} &= r^{{2}} \sin^2(\theta)
                \end{align*}
                % update metric gUU
            """),
            ['gDD', 'gdet', 'gUU', 'deltaUU', 'r_s', 'r', 'theta']
        )
        self.assertEqual(str(gUU),
            '[[-1 + r_s/r, 0, 0, 0], [0, 1/(1 - r_s/r), 0, 0], [0, 0, r**2, 0], [0, 0, 0, r**2*sin(theta)**2]]'
        )
        self.assertEqual(str(gdet),
            '(1 - r_s/r)/(r**4*(-1 + r_s/r)*sin(theta)**2)'
        )

    def test_example_6_1(self):
        self.assertEqual(
            parse(r"""
                % define anti01 FUU (4), const k;
                J^\mu = (4\pi k)^{-1} \nabla_\nu F^{\mu\nu}
            """),
            ['FUU', 'k', 'FUU_dD', 'gDD', 'gdet', 'gUU', 'gDD_dD', 'GammaUDD', 'FUU_cdD', 'JU']
        )

    def test_example_6_2(self):
        self.assertEqual(
            parse(r"""
                % define anti01 FUU (4), const k;
                J^\mu = (4\pi k)^{-1} \hat{\nabla}_\nu F^{\mu\nu}
            """),
            ['FUU', 'k', 'FUU_dD', 'ghatDD', 'ghatdet', 'ghatUU', 'ghatDD_dD', 'GammahatUDD', 'FUU_cdhatD', 'JU']
        )

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestParser('test_expression_1'))
    suite.addTest(TestParser('test_expression_2'))
    suite.addTest(TestParser('test_expression_3'))
    suite.addTest(TestParser('test_expression_4'))
    # suite.addTest(TestParser('test_expression_5'))
    # suite.addTest(TestParser('test_assignment_1'))
    # suite.addTest(TestParser('test_assignment_2'))
    suite.addTest(TestParser('test_example_1'))
    suite.addTest(TestParser('test_example_2'))
    suite.addTest(TestParser('test_example_3'))
    suite.addTest(TestParser('test_example_4'))
    # suite.addTest(TestParser('test_example_5'))
    # suite.addTest(TestParser('test_example_6_1'))
    # suite.addTest(TestParser('test_example_6_2'))
    result = unittest.TextTestRunner().run(suite)
    sys.exit(not result.wasSuccessful())
