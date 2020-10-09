""" Unit Testing for LaTeX Parser """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

# pylint: disable=import-error,protected-access
# import sys; sys.path.append('..')
from latex_parser import Tensor, OverrideWarning
from latex_parser import Parser, parse_expr, parse
from sympy import Function, Symbol, simplify
from warnings import filterwarnings
import unittest, sys

class TestParser(unittest.TestCase):
    """ Unit Testing for LaTeX Parser """

    def setUp(self):
        self.maxDiff = None
        Parser.ignore_override()

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
        expr = r'\partial_x (x^{{2}} + 2x)'
        self.assertEqual(
            str(parse_expr(expr).doit()),
            '2*x + 2'
        )

    def test_expression_4(self):
        tensor = Tensor(Function('Tensor')(Symbol('TUU'), Symbol('mu'), Symbol('nu')), 4)
        self.assertEqual(
            Parser._generate_covdrv(tensor, [('beta', 'D')], 'symbolic'),
            r'\nabla_\beta T^{\mu \nu} = \partial_\beta (T^{\mu \nu}) + \Gamma^\mu_{a \beta} (T^{a \nu}) + \Gamma^\nu_{a \beta} (T^{\mu a})'
        )
        tensor = Tensor(Function('Tensor')(Symbol('TUD'), Symbol('mu'), Symbol('nu')), 4)
        self.assertEqual(
            Parser._generate_covdrv(tensor, [('beta', 'D')], 'symbolic'),
            r'\nabla_\beta T^\mu_\nu = \partial_\beta (T^\mu_\nu) + \Gamma^\mu_{a \beta} (T^a_\nu) - \Gamma^a_{\nu \beta} (T^\mu_a)'
        )
        tensor = Tensor(Function('Tensor')(Symbol('TDD'), Symbol('mu'), Symbol('nu')), 4)
        self.assertEqual(
            Parser._generate_covdrv(tensor, [('beta', 'D')], 'symbolic'),
            r'\nabla_\beta T_{\mu \nu} = \partial_\beta (T_{\mu \nu}) - \Gamma^a_{\mu \beta} (T_{a \nu}) - \Gamma^a_{\nu \beta} (T_{\mu a})'
        )

    def test_expression_5(self):
        tensor = Tensor(Function('Tensor')(Symbol('vU'), Symbol('mu')), 4)
        self.assertEqual(
            Parser._generate_covdrv(tensor, [('a', 'U'), ('b', 'D')], 'symbolic'),
            r'\nabla_a \nabla_b v^\mu = \partial_a (\partial_b (v^\mu) + \Gamma^\mu_{c b} (v^c)) + \Gamma^\mu_{c a} (\partial_b (v^c) + \Gamma^c_{d b} (v^d))'
        )

    def test_assignment_1(self):
        self.assertEqual(
            parse(r"""
                % define vU (2), wU (2);
                T^{ij}_k = \vphantom{_d} \partial_k (v^i w^j)
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
                T^i_k = \vphantom{_d} \partial_k (v^i w)
            """),
            ['vU', 'w', 'vU_dD', 'TUD']
        )
        self.assertEqual(str(TUD),
            '[[vU_dD00*w, vU_dD01*w], [vU_dD10*w, vU_dD11*w]]'
        )

    def test_assignment_3(self):
        self.assertEqual(
            parse(r"""
                % define vU (4);
                T^{\mu\nu} = \vphantom{_d} \nabla^\nu v^\mu
            """),
            ['vU', 'gDD', 'gdet', 'gUU', 'vU_dD', 'gDD_dD', 'GammaUDD', 'vU_cdD', 'vU_cdU', 'TUU']
        )

    def test_assignment_4(self):
        parse(r"""
            % define basis [x, y];
            % define deriv _d;
            % define vD (2);
            v_0 = x^{{2}} + 2x;
            v_1 = y\sqrt{x};
            T_{\mu\nu} = (\partial_\nu v_\mu)\,\partial^2_x (y\sqrt{x})
        """)
        # T_{\mu\nu} = (\partial_\nu v_\mu)\,\partial^2_x (x^{{2}} + 2x)
        self.assertEqual(str(TDD),
            '[[-vD_dD00*y/(4*x**(3/2)), -vD_dD01*y/(4*x**(3/2))], [-vD_dD10*y/(4*x**(3/2)), -vD_dD11*y/(4*x**(3/2))]]'
        )

    def test_assignment_5(self):
        parse(r"""
            % define basis [\theta, \phi];
            % define metric gDD (2), kronecker deltaDD (2);
            % define const r;
            % parse g_{\mu\nu} = \delta_{\mu\nu};
            \begin{align*}
                g_{0 0} &= r^{{2}} \\
                g_{1 1} &= r^{{2}} \sin^2(\theta)
            \end{align*}
            % update metric gDD;
            \begin{align*}
                R^\alpha_{\beta\mu\nu} &= \partial_\mu \Gamma^\alpha_{\beta\nu} - \partial_\nu \Gamma^\alpha_{\beta\mu} + \Gamma^\alpha_{\mu\gamma}\Gamma^\gamma_{\beta\nu} - \Gamma^\alpha_{\nu\sigma}\Gamma^\sigma_{\beta\mu} \\
                R_{\alpha\beta\mu\nu} &= g_{\alpha a} R^a_{\beta\mu\nu} \\
                R_{\beta\nu} &= R^\alpha_{\beta\alpha\nu} \\
                R &= g^{\beta\nu} R_{\beta\nu}
            \end{align*}
        """)
        self.assertEqual(0,
            simplify(GammaUDD[0][1][1] - parse_expr(r'-\sin(\theta)\cos(\theta)'))
        )
        self.assertEqual(0,
            simplify(GammaUDD[1][0][1] - GammaUDD[1][1][0])
        )
        self.assertEqual(0,
            simplify(GammaUDD[1][0][1] - parse_expr(r'\cos(\theta)/\sin(\theta)'))
        )
        self.assertEqual(0,
            simplify(RDDDD[0][1][0][1] - (-RDDDD[0][1][1][0]) + (-RDDDD[1][0][0][1]) - RDDDD[1][0][1][0])
        )
        self.assertEqual(0,
            simplify(RDDDD[0][1][0][1] - parse_expr(r'r^{{2}} \sin^2(\theta)'))
        )
        self.assertEqual(0,
            simplify(RDD[0][0] - 1)
        )
        self.assertEqual(0,
            simplify(RDD[1][1] - parse_expr(r'\sin^2(\theta)'))
        )
        self.assertEqual(0,
            simplify(RDD[0][1] - RDD[1][0])
        )
        self.assertEqual(0,
            simplify(RDD[0][1] - 0)
        )
        self.assertEqual(0,
            simplify(R - parse_expr(r'2/r^{{2}}'))
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
                % define permutation epsilonDDD (3);
                % define vU (3), wU (3);
                u_i = \epsilon_{ijk} v^j w^k
            """),
            ['epsilonDDD', 'vU', 'wU', 'uD']
        )
        self.assertEqual(str(uD),
            '[vU1*wU2 - vU2*wU1, -vU0*wU2 + vU2*wU0, vU0*wU1 - vU1*wU0]'
        )

    def test_example_4_1(self):
        self.assertEqual(
            parse(r"""
                % define anti01 FUU (4), const k;
                J^\mu = (4\pi k)^{-1} \vphantom{_d} \nabla_\nu F^{\mu\nu}
            """),
            ['FUU', 'k', 'FUU_dD', 'gDD', 'gdet', 'gUU', 'gDD_dD', 'GammaUDD', 'FUU_cdD', 'JU']
        )

    def test_example_4_2(self):
        self.assertEqual(
            parse(r"""
                % define anti01 FUU (4), const k;
                J^\mu = (4\pi k)^{-1} \vphantom{_d} \hat{\nabla}_\nu F^{\mu\nu}
            """),
            ['FUU', 'k', 'FUU_dD', 'ghatDD', 'ghatdet', 'ghatUU', 'ghatDD_dD', 'GammahatUDD', 'FUU_cdhatD', 'JU']
        )

    def test_example_5(self):
        parse(r"""
            % define basis [t, r, \theta, \phi];
            % define metric gDD (4), kronecker deltaDD (4);
            % define const G, const M;
            % parse g_{\mu\nu} = \delta_{\mu\nu};
            \begin{align}
                g_{0 0} &= -\left(1 - \frac{2GM}{r}\right) \\
                g_{1 1} &=  \left(1 - \frac{2GM}{r}\right)^{-1} \\
                g_{2 2} &= r^{{2}} \\
                g_{3 3} &= r^{{2}} \sin^2\theta
            \end{align}
            % update metric gDD;
            \begin{align}
                R^\alpha{}_{\beta\mu\nu} &= \partial_\mu \Gamma^\alpha_{\beta\nu} - \partial_\nu \Gamma^\alpha_{\beta\mu} + \Gamma^\alpha_{\mu\gamma}\Gamma^\gamma_{\beta\nu} - \Gamma^\alpha_{\nu\sigma}\Gamma^\sigma_{\beta\mu} \\
                R^{\alpha\beta\mu\nu} &= g^{\beta a} g^{\mu b} g^{\nu c} R^\alpha_{a b c} \\
                R_{\alpha\beta\mu\nu} &= g_{\alpha a} R^a_{\beta\mu\nu} \\
                K &= R^{\alpha\beta\mu\nu} R_{\alpha\beta\mu\nu}
            \end{align}
        """)
        self.assertEqual(0,
            simplify(K - parse_expr(r'48(GM)^2/r^{{6}}'))
        )
        self.assertEqual(0,
            simplify(GammaUDD[0][0][1] - GammaUDD[0][1][0])
        )
        self.assertEqual(0,
            simplify(GammaUDD[0][0][1] - parse_expr(r'-GM/(r^{{2}}(2GM/r - 1))'))
        )
        self.assertEqual(0,
            simplify(GammaUDD[1][0][0] - parse_expr(r'GM(1 - 2GM/r)/r^{{2}}'))
        )
        self.assertEqual(0,
            simplify(GammaUDD[1][1][1] - parse_expr(r'-GM/(r^{{2}}(1 - 2GM/r))'))
        )
        self.assertEqual(0,
            simplify(GammaUDD[1][3][3] - parse_expr(r'-r(1 - 2GM/r)\sin^2\theta'))
        )
        self.assertEqual(0,
            simplify(GammaUDD[2][1][2] - GammaUDD[2][2][1])
        )
        self.assertEqual(0,
            simplify(GammaUDD[2][1][2] - parse_expr(r'1/r'))
        )
        self.assertEqual(0,
            simplify(GammaUDD[2][3][3] - parse_expr(r'-\sin\theta\cos\theta'))
        )
        self.assertEqual(0,
            simplify(GammaUDD[2][1][3] - GammaUDD[2][3][1])
        )
        self.assertEqual(0,
            simplify(GammaUDD[3][1][3] - parse_expr(r'1/r'))
        )
        self.assertEqual(0,
            simplify(GammaUDD[3][2][3] - GammaUDD[3][3][2])
        )
        self.assertEqual(0,
            simplify(GammaUDD[3][2][3] - parse_expr(r'\cos\theta/\sin\theta'))
        )


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestParser('test_expression_1'))
    suite.addTest(TestParser('test_expression_2'))
    suite.addTest(TestParser('test_expression_3'))
    suite.addTest(TestParser('test_expression_4'))
    suite.addTest(TestParser('test_expression_5'))
    suite.addTest(TestParser('test_assignment_1'))
    suite.addTest(TestParser('test_assignment_2'))
    suite.addTest(TestParser('test_assignment_3'))
    suite.addTest(TestParser('test_assignment_4'))
    suite.addTest(TestParser('test_assignment_5'))
    suite.addTest(TestParser('test_example_1'))
    suite.addTest(TestParser('test_example_2'))
    suite.addTest(TestParser('test_example_3'))
    suite.addTest(TestParser('test_example_4_1'))
    suite.addTest(TestParser('test_example_4_2'))
    suite.addTest(TestParser('test_example_5'))
    result = unittest.TextTestRunner().run(suite)
    sys.exit(not result.wasSuccessful())
