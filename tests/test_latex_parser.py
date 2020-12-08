""" Unit Testing for LaTeX Parser """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

# pylint: disable = import-error, protected-access
# import sys; sys.path.append('..')
from latex_parser import Tensor, Parser, parse_expr, parse
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
        expr = r'e^{{\ln x}} + \sin(\sin^{-1} y) - \tanh(xy)'
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
        function = Function('Tensor')(Symbol('T'))
        tensor = Tensor(function, 4)
        self.assertEqual(
            Parser._generate_covdrv(tensor, function, [('beta', 'D')]),
            r'\nabla_\beta T = \partial_\beta (T)'
        )
        function = Function('Tensor')(Symbol('TUU'), Symbol('mu'), Symbol('nu'))
        tensor = Tensor(function, 4)
        self.assertEqual(
            Parser._generate_covdrv(tensor, function, [('beta', 'D')]),
            r'\nabla_\beta T^{\mu \nu} = \partial_\beta (T^{\mu \nu}) + \Gamma^\mu_{a \beta} (T^{a \nu}) + \Gamma^\nu_{a \beta} (T^{\mu a})'
        )
        function = Function('Tensor')(Symbol('TUD'), Symbol('mu'), Symbol('nu'))
        tensor = Tensor(function, 4)
        self.assertEqual(
            Parser._generate_covdrv(tensor, function, [('beta', 'D')]),
            r'\nabla_\beta T^\mu_\nu = \partial_\beta (T^\mu_\nu) + \Gamma^\mu_{a \beta} (T^a_\nu) - \Gamma^a_{\nu \beta} (T^\mu_a)'
        )
        function = Function('Tensor')(Symbol('TDD'), Symbol('mu'), Symbol('nu'))
        tensor = Tensor(function, 4)
        self.assertEqual(
            Parser._generate_covdrv(tensor, function, [('beta', 'D')]),
            r'\nabla_\beta T_{\mu \nu} = \partial_\beta (T_{\mu \nu}) - \Gamma^a_{\mu \beta} (T_{a \nu}) - \Gamma^a_{\nu \beta} (T_{\mu a})'
        )

    def test_expression_5(self):
        function = Function('Tensor')(Symbol('vU'), Symbol('mu'))
        tensor = Tensor(function, 4)
        self.assertEqual(
            Parser._generate_covdrv(tensor, function, [('a', 'D'), ('b', 'D')]),
            r'\nabla_a \nabla_b v^\mu = \partial_a (\partial_b (v^\mu) + \Gamma^\mu_{c b} (v^c)) + \Gamma^\mu_{c a} (\partial_b (v^c) + \Gamma^c_{d b} (v^d)) - \Gamma^c_{b a} (\partial_c (v^\mu) + \Gamma^\mu_{b c} (v^b))'
        )

    def test_expression_6(self):
        function = Function('Tensor')(Symbol('g'))
        tensor = Tensor(function, dimension=3, weight=2)
        self.assertEqual(
            Parser._generate_liedrv(tensor, function, 'beta'),
            r'\mathcal{L}_\text{beta} g = \text{beta}^a \partial_a g + (2)(\partial_a \text{beta}^a) g'
        )
        function = Function('Tensor')(Symbol('gUU'), Symbol('i'), Symbol('j'))
        tensor = Tensor(function, 3)
        self.assertEqual(
            Parser._generate_liedrv(tensor, function, 'beta'),
            r'\mathcal{L}_\text{beta} g^{i j} = \text{beta}^a \partial_a g^{i j} - (\partial_a \text{beta}^i) g^{a j} - (\partial_a \text{beta}^j) g^{i a}'
        )
        function = Function('Tensor')(Symbol('gUD'), Symbol('i'), Symbol('j'))
        tensor = Tensor(function, 3)
        self.assertEqual(
            Parser._generate_liedrv(tensor, function, 'beta'),
            r'\mathcal{L}_\text{beta} g^i_j = \text{beta}^a \partial_a g^i_j - (\partial_a \text{beta}^i) g^a_j + (\partial_j \text{beta}^a) g^i_a'
        )
        function = Function('Tensor')(Symbol('gDD'), Symbol('i'), Symbol('j'))
        tensor = Tensor(function, 3)
        self.assertEqual(
            Parser._generate_liedrv(tensor, function, 'beta'),
            r'\mathcal{L}_\text{beta} g_{i j} = \text{beta}^a \partial_a g_{i j} + (\partial_i \text{beta}^a) g_{a j} + (\partial_j \text{beta}^a) g_{i a}'
        )

    def test_assignment_1(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % define vU (2D), wU (2D)
                T^{ab}_c = \vphantom{numeric} \partial_c (v^a w^b)
            """)),
            {'vU', 'wU', 'vU_dD', 'wU_dD', 'TUUD'}
        )
        self.assertEqual(str(TUUD[0][0][0]),
            'vU0*wU_dD00 + vU_dD00*wU0'
        )

    def test_assignment_2(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % define vU (2D), const w
                T^a_c = \vphantom{numeric} \partial_c (v^a w)
            """)),
            {'vU', 'vU_dD', 'TUD'}
        )
        self.assertEqual(str(TUD),
            '[[vU_dD00*w, vU_dD01*w], [vU_dD10*w, vU_dD11*w]]'
        )

    def test_assignment_3(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % define metric gDD (4D), vU (4D)
                T^{\mu\nu} = \vphantom{numeric} \nabla^\nu v^\mu
            """)),
            {'gUU', 'gDD', 'vU', 'vU_dD', 'gDD_dD', 'GammaUDD', 'vU_cdD', 'vU_cdU', 'TUU'}
        )

    def test_assignment_4(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % define basis [x, y]
                % define symbolic uD (2D), wD (2D)
                u_0 = x^{{2}} + 2x \\
                u_1 = y\sqrt{x} \\
                v_a = u_a + w_a \\
                T_{ab} = \partial^2_x v_0\, \vphantom{numeric} \partial_b v_a
            """)),
            {'uD', 'wD', 'vD', 'vD_dD', 'wD_dD', 'TDD'}
        )
        self.assertEqual(str(TDD),
            '[[2*wD_dD00 + 4*x + 4, 2*wD_dD01], [2*wD_dD10 + y/sqrt(x), 2*wD_dD11 + 2*sqrt(x)]]'
        )

    def test_assignment_5(self):
        Parser.clear_namespace()
        parse(r"""
            % define basis [\theta, \phi]
            % define const r, deltaDD (2D)
            % parse g_{\mu\nu} = \delta_{\mu\nu}
            \begin{align*}
                g_{0 0} &= r^{{2}} \\
                g_{1 1} &= r^{{2}} \sin^2(\theta)
            \end{align*}
            % assign metric gDD
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
            set(parse(r"""
                % define metric gDD (4D)
                \gamma_{ij} = g_{ij}
            """)),
            {'gUU', 'gDD', 'gammaDD'}
        )
        self.assertEqual(str(gammaDD),
            '[[gDD11, gDD12, gDD13], [gDD12, gDD22, gDD23], [gDD13, gDD23, gDD33]]'
        )
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % define range [i-j] = 0:1
                % define nosym TUU (3D), vD (2D)
                w^\mu = T^{\mu i} v_i
            """)),
            {'TUU', 'vD', 'wU'}
        )
        self.assertEqual(str(wU),
            '[TUU01*vD0 + TUU02*vD1, TUU11*vD0 + TUU12*vD1, TUU21*vD0 + TUU22*vD1]'
        )

    def test_example_1(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % define nosym hUD (4D)
                h = h^\mu{}_\mu
            """)),
            {'hUD', 'h'}
        )
        self.assertEqual(str(h),
            'hUD00 + hUD11 + hUD22 + hUD33'
        )

    def test_example_2(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % define metric gUU (3D), vD (3D)
                v^\mu = g^{\mu\nu} v_\nu
            """)),
            {'gDD', 'gUU', 'vD', 'vU'}
        )
        self.assertEqual(str(vU),
            '[gUU00*vD0 + gUU01*vD1 + gUU02*vD2, gUU01*vD0 + gUU11*vD1 + gUU12*vD2, gUU02*vD0 + gUU12*vD1 + gUU22*vD2]'
        )

    def test_example_3(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % define epsilonDDD (3D)
                % define vU (3D), wU (3D)
                u_i = \epsilon_{ijk} v^j w^k
            """)),
            {'epsilonDDD', 'vU', 'wU', 'uD'}
        )
        self.assertEqual(str(uD),
            '[vU1*wU2 - vU2*wU1, -vU0*wU2 + vU2*wU0, vU0*wU1 - vU1*wU0]'
        )

    def test_example_4_1(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % define anti01 FUU (4D), metric gDD (4D), const k
                J^\mu = (4\pi k)^{-1} \vphantom{numeric} \nabla_\nu F^{\mu\nu}
            """)),
            {'FUU', 'gUU', 'gDD', 'FUU_dD', 'gDD_dD', 'GammaUDD', 'FUU_cdD', 'JU'}
        )

    def test_example_4_2(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % define anti01 FUU (4D), metric ghatDD (4D), const k
                J^\mu = (4\pi k)^{-1} \vphantom{numeric} \hat{\nabla}_\nu F^{\mu\nu}
            """)),
            {'FUU', 'ghatUU', 'ghatDD', 'FUU_dD', 'ghatDD_dD', 'GammahatUDD', 'FUU_cdhatD', 'JU'}
        )

    def test_example_5_1(self):
        Parser.clear_namespace()
        parse(r"""
            % define deltaDD (4D)
            % define const G, const M
            % parse g_{\mu\nu} = \delta_{\mu\nu}
            \begin{align}
                g_{0 0} &= -\left(1 - \frac{2GM}{r}\right) \\
                g_{1 1} &=  \left(1 - \frac{2GM}{r}\right)^{-1} \\
                g_{2 2} &= r^{{2}} \\
                g_{3 3} &= r^{{2}} \sin^2\theta
            \end{align}
            % assign metric gDD
        """)
        self.assertEqual(str(gDD[0][0]),
            '2*G*M/r - 1'
        )
        self.assertEqual(str(gDD[1][1]),
            '1/(-2*G*M/r + 1)'
        )
        self.assertEqual(str(gDD[2][2]),
            'r**2'
        )
        self.assertEqual(str(gDD[3][3]),
            'r**2*sin(theta)**2'
        )
        self.assertEqual(str(gdet),
            'r**4*(2*G*M/r - 1)*sin(theta)**2/(-2*G*M/r + 1)'
        )

    def test_example_5_2(self):
        parse(r"""
            % define basis [t, r, \theta, \phi]
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

    def test_example_6_1(self):
        parse(r"""
            % define basis [r, \theta, \phi]
            \begin{align}
                \gamma_{ij} &= g_{ij} \\
                % assign metric gammaDD
                \beta_i &= g_{0 i} \\
                \alpha &= \sqrt{\gamma^{ij}\beta_i\beta_j - g_{0 0}} \\
                K_{ij} &= \frac{1}{2\alpha}\left(\nabla_i \beta_j + \nabla_j \beta_i\right) \\
                K &= \gamma^{ij} K_{ij}
            \end{align}
        """)
        self.assertEqual(True,
            all(component == 0 for component in (row for row in simplify(Matrix(KDD))))
        )

    def test_example_6_2(self):
        parse(r"""
            \begin{align}
                K^{ij} &= \gamma^{ik} \gamma^{jl} K_{kl} \\
                R_{ij} &= \partial_k \Gamma^k_{ij} - \partial_j \Gamma^k_{ik}
                    + \Gamma^k_{ij}\Gamma^l_{kl} - \Gamma^l_{ik}\Gamma^k_{lj} \\
                R &= \gamma^{ij} R_{ij} \\
                E &= \frac{1}{16\pi}\left(R + K^{{2}} - K_{ij}K^{ij}\right) \\
                p_i &= \frac{1}{8\pi}\left(D_j \gamma^{jk} K_{ki} - D_i K\right)
            \end{align}
        """)
        self.assertEqual(simplify(E), 0)
        self.assertEqual(True,
            all(component == 0 for component in pD)
        )

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(TestParser('test_expression_1'))
    suite.addTest(TestParser('test_expression_2'))
    suite.addTest(TestParser('test_expression_3'))
    suite.addTest(TestParser('test_expression_4'))
    suite.addTest(TestParser('test_expression_5'))
    suite.addTest(TestParser('test_expression_6'))
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
    suite.addTest(TestParser('test_example_5_1'))
    suite.addTest(TestParser('test_example_5_2'))
    suite.addTest(TestParser('test_example_6_1'))
    suite.addTest(TestParser('test_example_6_2'))
    result = unittest.TextTestRunner().run(suite)
    sys.exit(not result.wasSuccessful())
