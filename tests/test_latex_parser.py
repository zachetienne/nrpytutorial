""" Unit Testing for LaTeX Parser """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

# pylint: disable=import-error,protected-access
# import sys; sys.path.append('..')
from latex_parser import Tensor, OverrideWarning
from latex_parser import Parser, parse_expr, parse
from sympy import Function, Symbol, Matrix, simplify
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
        Parser.clear_namespace()
        self.assertEqual(
            parse(r"""
                % define vU (2D), wU (2D);
                T^{ab}_c = \vphantom{_d} \partial_c (v^a w^b)
            """),
            ('vU', 'wU', 'vU_dD', 'wU_dD', 'TUUD')
        )
        self.assertEqual(str(TUUD[0][0][0]),
            'vU0*wU_dD00 + vU_dD00*wU0'
        )

    def test_assignment_2(self):
        Parser.clear_namespace()
        self.assertEqual(
            parse(r"""
                % define vU (2D), const w;
                T^a_c = \vphantom{_d} \partial_c (v^a w)
            """),
            ('vU', 'w', 'vU_dD', 'TUD')
        )
        self.assertEqual(str(TUD),
            '[[vU_dD00*w, vU_dD01*w], [vU_dD10*w, vU_dD11*w]]'
        )

    def test_assignment_3(self):
        Parser.clear_namespace()
        self.assertEqual(
            parse(r"""
                % define metric gDD (4D), vU (4D);
                T^{\mu\nu} = \vphantom{_d} \nabla^\nu v^\mu
            """),
            ('gUU', 'gdet', 'gDD', 'vU', 'vU_dD', 'gDD_dD', 'GammaUDD', 'vU_cdD', 'vU_cdU', 'TUU')
        )

    def test_assignment_4(self):
        Parser.clear_namespace()
        self.assertEqual(
            parse(r"""
                % define basis [x, y], deriv _d;
                % define vD (2D);
                v_0 = x^{{2}} + 2x;
                v_1 = y\sqrt{x};
                T_{\mu\nu} = (\partial_\nu v_\mu)\,\partial^2_x v_0
            """),
            ('vD', 'x', 'y', 'vD_dD', 'TDD')
        )
        self.assertEqual(str(TDD),
            '[[2*vD_dD00, 2*vD_dD01], [2*vD_dD10, 2*vD_dD11]]'
        )

    def test_assignment_5(self):
        Parser.clear_namespace()
        parse(r"""
            % define basis [\theta, \phi];
            % define const r, kronecker deltaDD (2D);
            % parse g_{\mu\nu} = \delta_{\mu\nu};
            \begin{align*}
                g_{0 0} &= r^{{2}} \\
                g_{1 1} &= r^{{2}} \sin^2(\theta)
            \end{align*}
            % assign metric gDD;
            \begin{align*}
                R^\alpha_{\beta\mu\nu} &= \partial_\mu \Gamma^\alpha_{\beta\nu} - \partial_\nu \Gamma^\alpha_{\beta\mu} + \Gamma^\alpha_{\mu\gamma}\Gamma^\gamma_{\beta\nu} - \Gamma^\alpha_{\nu\sigma}\Gamma^\sigma_{\beta\mu} \\
                R_{\alpha\beta\mu\nu} &= g_{\alpha a} R^a_{\beta\mu\nu} \\
                R_{\beta\nu} &= R^\alpha_{\beta\alpha\nu} \\
                R &= g^{\beta\nu} R_{\beta\nu}
            \end{align*}
        """)
        self.assertEqual(str(GammaUDD[0][1][1]),
            '-sin(theta)*cos(theta)'
        )
        self.assertEqual(0,
            simplify(GammaUDD[1][0][1] - GammaUDD[1][1][0])
        )
        self.assertEqual(str(GammaUDD[1][0][1]),
            'cos(theta)/sin(theta)'
        )
        self.assertEqual(0,
            simplify(RDDDD[0][1][0][1] - (-RDDDD[0][1][1][0]) + (-RDDDD[1][0][0][1]) - RDDDD[1][0][1][0])
        )
        self.assertEqual(str(RDDDD[0][1][0][1]),
            'r**2*sin(theta)**2'
        )
        self.assertEqual(RDD[0][0], 1)
        self.assertEqual(str(RDD[1][1]),
            'sin(theta)**2'
        )
        self.assertEqual(0,
            simplify(RDD[0][1] - RDD[1][0])
        )
        self.assertEqual(RDD[0][1], 0)
        self.assertEqual(str(R),
            '2/r**2'
        )

    def test_assignment_6(self):
        Parser.clear_namespace()
        self.assertEqual(
            parse(r"""
                % define metric gDD (4D);
                \gamma_{ij} = g_{ij}
            """),
            ('gUU', 'gdet', 'gDD', 'gammaDD')
        )
        self.assertEqual(str(gammaDD),
            '[[gDD11, gDD12, gDD13], [gDD12, gDD22, gDD23], [gDD13, gDD23, gDD33]]'
        )
        Parser.clear_namespace()
        self.assertEqual(
            parse(r"""
                % define index [i-j] = 0:1;
                % define nosym TUU (3D), vD (2D);
                w^\mu = T^{\mu i} v_i
            """),
            ('TUU', 'vD', 'wU')
        )
        self.assertEqual(str(wU),
            '[TUU01*vD0 + TUU02*vD1, TUU11*vD0 + TUU12*vD1, TUU21*vD0 + TUU22*vD1]'
        )

    def test_example_1(self):
        Parser.clear_namespace()
        self.assertEqual(
            parse(r"""
                % define nosym hUD (4D);
                h = h^\mu{}_\mu
            """),
            ('hUD', 'h')
        )
        self.assertEqual(str(h),
            'hUD00 + hUD11 + hUD22 + hUD33'
        )

    def test_example_2(self):
        Parser.clear_namespace()
        self.assertEqual(
            parse(r"""
                % define metric gUU (3D), vD (3D);
                v^\mu = g^{\mu\nu} v_\nu
            """),
            ('gDD', 'gdet', 'gUU', 'vD', 'vU')
        )
        self.assertEqual(str(vU),
            '[gUU00*vD0 + gUU01*vD1 + gUU02*vD2, gUU01*vD0 + gUU11*vD1 + gUU12*vD2, gUU02*vD0 + gUU12*vD1 + gUU22*vD2]'
        )

    def test_example_3(self):
        Parser.clear_namespace()
        self.assertEqual(
            parse(r"""
                % define permutation epsilonDDD (3D);
                % define vU (3D), wU (3D);
                u_i = \epsilon_{ijk} v^j w^k
            """),
            ('epsilonDDD', 'vU', 'wU', 'uD')
        )
        self.assertEqual(str(uD),
            '[vU1*wU2 - vU2*wU1, -vU0*wU2 + vU2*wU0, vU0*wU1 - vU1*wU0]'
        )

    def test_example_4_1(self):
        Parser.clear_namespace()
        self.assertEqual(
            parse(r"""
                % define anti01 FUU (4D), metric gDD (4D), const k;
                J^\mu = (4\pi k)^{-1} \vphantom{_d} \nabla_\nu F^{\mu\nu}
            """),
            ('FUU', 'gUU', 'gdet', 'gDD', 'k', 'FUU_dD', 'gDD_dD', 'GammaUDD', 'FUU_cdD', 'JU')
        )

    def test_example_4_2(self):
        Parser.clear_namespace()
        self.assertEqual(
            parse(r"""
                % define anti01 FUU (4D), metric ghatDD (4D), const k;
                J^\mu = (4\pi k)^{-1} \vphantom{_d} \hat{\nabla}_\nu F^{\mu\nu}
            """),
            ('FUU', 'ghatUU', 'ghatdet', 'ghatDD', 'k', 'FUU_dD', 'ghatDD_dD', 'GammahatUDD', 'FUU_cdhatD', 'JU')
        )

    def test_example_5(self):
        Parser.clear_namespace()
        parse(r"""
            % define kronecker deltaDD (4D);
            % define const G, const M;
            % parse g_{\mu\nu} = \delta_{\mu\nu};
            \begin{align}
                g_{0 0} &= -\left(1 - \frac{2GM}{r}\right) \\
                g_{1 1} &=  \left(1 - \frac{2GM}{r}\right)^{-1} \\
                g_{2 2} &= r^{{2}} \\
                g_{3 3} &= r^{{2}} \sin^2\theta
            \end{align}
            % assign metric gDD
        """)
        self.assertEqual(str(gdet),
            'r**4*(2*G*M/r - 1)*sin(theta)**2/(-2*G*M/r + 1)'
        )

    def test_example_6(self):
        parse(r"""
            % define basis [t, r, \theta, \phi];
            \begin{align}
                R^\alpha{}_{\beta\mu\nu} &= \partial_\mu \Gamma^\alpha_{\beta\nu} - \partial_\nu \Gamma^\alpha_{\beta\mu} + \Gamma^\alpha_{\mu\gamma}\Gamma^\gamma_{\beta\nu} - \Gamma^\alpha_{\nu\sigma}\Gamma^\sigma_{\beta\mu} \\
                R^{\alpha\beta\mu\nu} &= g^{\beta a} g^{\mu b} g^{\nu c} R^\alpha_{a b c} \\
                R_{\alpha\beta\mu\nu} &= g_{\alpha a} R^a_{\beta\mu\nu} \\
                K &= R^{\alpha\beta\mu\nu} R_{\alpha\beta\mu\nu} \\
                R_{\beta\nu} &= R^\alpha_{\beta\alpha\nu} \\
                R &= g^{\beta\nu} R_{\beta\nu} \\
                G_{\beta\nu} &= R_{\beta\nu} - \frac{1}{2}g_{\beta\nu}R
            \end{align}
        """)
        self.assertEqual(0,
            simplify(GammaUDD[0][0][1] - GammaUDD[0][1][0])
        )
        self.assertEqual(str(GammaUDD[0][0][1]),
            '-G*M/(r**2*(2*G*M/r - 1))'
        )
        self.assertEqual(str(GammaUDD[1][0][0]),
            'G*M*(-2*G*M/r + 1)/r**2'
        )
        self.assertEqual(str(GammaUDD[1][1][1]),
            '-G*M/(r**2*(-2*G*M/r + 1))'
        )
        self.assertEqual(str(GammaUDD[1][3][3]),
            '-r*(-2*G*M/r + 1)*sin(theta)**2'
        )
        self.assertEqual(0,
            simplify(GammaUDD[2][1][2] - GammaUDD[2][2][1])
        )
        self.assertEqual(str(GammaUDD[2][1][2]),
            '1/r'
        )
        self.assertEqual(str(GammaUDD[2][3][3]),
            '-sin(theta)*cos(theta)'
        )
        self.assertEqual(0,
            simplify(GammaUDD[2][1][3] - GammaUDD[2][3][1])
        )
        self.assertEqual(str(GammaUDD[3][1][3]),
            '1/r'
        )
        self.assertEqual(0,
            simplify(GammaUDD[3][2][3] - GammaUDD[3][3][2])
        )
        self.assertEqual(str(GammaUDD[3][2][3]),
            'cos(theta)/sin(theta)'
        )
        self.assertEqual(str(simplify(K)),
            '48*G**2*M**2/r**6'
        )
        self.assertEqual(simplify(R), 0)
        self.assertEqual(True,
            all(component == 0 for component in (row for row in simplify(Matrix(GDD))))
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
    suite.addTest(TestParser('test_assignment_6'))
    suite.addTest(TestParser('test_example_1'))
    suite.addTest(TestParser('test_example_2'))
    suite.addTest(TestParser('test_example_3'))
    suite.addTest(TestParser('test_example_4_1'))
    suite.addTest(TestParser('test_example_4_2'))
    suite.addTest(TestParser('test_example_5'))
    suite.addTest(TestParser('test_example_6'))
    result = unittest.TextTestRunner().run(suite)
    sys.exit(not result.wasSuccessful())
