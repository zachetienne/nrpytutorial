""" Unit Testing for LaTeX Parser """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

# pylint: disable = import-error, protected-access, exec-used
import sys; sys.path.append('..')
from latex_parser import Tensor, Parser, parse_expr, parse
from sympy import Function, Symbol, Matrix, simplify
from UnitTesting.assert_equal import assert_equal
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
        self.assertEqual(
            Parser._generate_covdrv(function, 'beta'),
            r'\nabla_\beta T = \partial_\beta T'
        )
        function = Function('Tensor')(Symbol('TUU'), Symbol('mu'), Symbol('nu'))
        self.assertEqual(
            Parser._generate_covdrv(function, 'beta'),
            r'\nabla_\beta T^{\mu \nu} = \partial_\beta T^{\mu \nu} + \text{Gamma}^\mu_{a \beta} (T^{a \nu}) + \text{Gamma}^\nu_{a \beta} (T^{\mu a})'
        )
        function = Function('Tensor')(Symbol('TUD'), Symbol('mu'), Symbol('nu'))
        self.assertEqual(
            Parser._generate_covdrv(function, 'beta'),
            r'\nabla_\beta T^\mu_\nu = \partial_\beta T^\mu_\nu + \text{Gamma}^\mu_{a \beta} (T^a_\nu) - \text{Gamma}^a_{\nu \beta} (T^\mu_a)'
        )
        function = Function('Tensor')(Symbol('TDD'), Symbol('mu'), Symbol('nu'))
        self.assertEqual(
            Parser._generate_covdrv(function, 'beta'),
            r'\nabla_\beta T_{\mu \nu} = \partial_\beta T_{\mu \nu} - \text{Gamma}^a_{\mu \beta} (T_{a \nu}) - \text{Gamma}^a_{\nu \beta} (T_{\mu a})'
        )

    def test_expression_5(self):
        Parser.clear_namespace()
        parse(r"""
            % vardef -numeric -metric 'gDD' (4D)
            % vardef -numeric 'vU' (4D)
            T^\mu_b = \nabla_b v^\mu
        """)
        function = Parser._namespace['vU_cdD'].equation[0]
        self.assertEqual(
            Parser._generate_covdrv(function, 'a'),
            r'\nabla_a \nabla_b v^\mu = \partial_a \nabla_b v^\mu + \text{Gamma}^\mu_{c a} (\nabla_b v^c) - \text{Gamma}^c_{b a} (\nabla_c v^\mu)'
        )

    def test_expression_6(self):
        function = Function('Tensor')(Symbol('g'))
        self.assertEqual(
            Parser._generate_liedrv(function, 'beta', 2),
            r'\mathcal{L}_\text{beta} g = \text{beta}^a \partial_a g + (2)(\partial_a \text{beta}^a) g'
        )
        function = Function('Tensor')(Symbol('gUU'), Symbol('i'), Symbol('j'))
        self.assertEqual(
            Parser._generate_liedrv(function, 'beta'),
            r'\mathcal{L}_\text{beta} g^{i j} = \text{beta}^a \partial_a g^{i j} - (\partial_a \text{beta}^i) g^{a j} - (\partial_a \text{beta}^j) g^{i a}'
        )
        function = Function('Tensor')(Symbol('gUD'), Symbol('i'), Symbol('j'))
        self.assertEqual(
            Parser._generate_liedrv(function, 'beta'),
            r'\mathcal{L}_\text{beta} g^i_j = \text{beta}^a \partial_a g^i_j - (\partial_a \text{beta}^i) g^a_j + (\partial_j \text{beta}^a) g^i_a'
        )
        function = Function('Tensor')(Symbol('gDD'), Symbol('i'), Symbol('j'))
        self.assertEqual(
            Parser._generate_liedrv(function, 'beta'),
            r'\mathcal{L}_\text{beta} g_{i j} = \text{beta}^a \partial_a g_{i j} + (\partial_i \text{beta}^a) g_{a j} + (\partial_j \text{beta}^a) g_{i a}'
        )

    def test_srepl_macro(self):
        Parser.clear_namespace()
        parse(r"""
            % srepl -global "<1>'" -> "\text{<1>prime}"
            % srepl -global "\text{<1..>}_<2>" -> "\text{(<1..>)<2>}"
            % srepl -global "<1>_{<2>}" -> "<1>_<2>", "<1>_<2>" -> "\text{<1>_<2>}"
            % srepl -global "\text{(<1..>)<2>}" -> "\text{<1..>_<2>}"
            % srepl -global "<1>^{<2>}" -> "<1>^<2>", "<1>^<2>" -> "<1>^{{<2>}}"
        """)
        expr = r"x_n^4 + x'_n \exp(x_n y_n^2)"
        self.assertEqual(
            str(parse_expr(expr)),
            "x_n**4 + xprime_n*exp(x_n*y_n**2)"
        )
        Parser.clear_namespace()
        parse(r""" % srepl -global "<1>'^{<2..>}" -> "\text{<1>prime}" """)
        expr = r"v'^{label}"
        self.assertEqual(
            str(parse_expr(expr)),
            "vprime"
        )

    def test_assignment_1(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % vardef -numeric 'vU' (2D), 'wU' (2D)
                % keydef index [a-z] (2D)
                T^{ab}_c = \partial_c (v^a w^b)
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
                % vardef -const 'w'
                % vardef -numeric 'vU' (2D)
                % keydef index [a-z] (2D)
                T^a_c = \partial_c (v^a w)
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
                % vardef -numeric -metric 'gDD' (4D)
                % vardef -numeric 'vU' (4D)
                % keydef index [a-z] (4D)
                T^{ab} = \nabla^b v^a
            """)),
            {'gUU', 'gdet', 'epsilonUUUU', 'gDD', 'vU', 'vU_dD', 'gDD_dD', 'GammaUDD', 'vU_cdD', 'vU_cdU', 'TUU'}
        )

    def test_assignment_4(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % keydef basis [x, y]
                % vardef 'uD' (2D), 'wD' (2D)
                % keydef index [a-z] (2D)
                u_0 = x^{{2}} + 2x \\
                u_1 = y\sqrt{x} \\
                v_a = u_a + w_a \\
                % assign -numeric 'wD', 'vD'
                T_{ab} = \partial^2_x v_0 (\partial_b v_a)
            """)),
            {'uD', 'wD', 'vD', 'vD_dD', 'wD_dD', 'TDD'}
        )
        self.assertEqual(str(TDD),
            '[[2*wD_dD00 + 4*x + 4, 2*wD_dD01], [2*wD_dD10 + y/sqrt(x), 2*wD_dD11 + 2*sqrt(x)]]'
        )

    def test_assignment_5(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % keydef basis [x, y]
                % vardef 'uD' (2D), 'wD' (2D)
                % assign -symbolic <H> 'uD'
                % keydef index [a-z] (2D)
                u_0 = x^{{2}} + 2x \\
                u_1 = y\sqrt{x} \\
                v_a = u_a + w_a \\
                T_{bc} = \partial^2_x v_0 (\vphantom{numeric} \partial_c v_b)
            """)),
            {'uD', 'wD', 'vD', 'vD_dD', 'wD_dD', 'TDD'}
        )
        self.assertEqual(str(TDD),
            '[[2*wD_dD00 + 4*x + 4, 2*wD_dD01], [2*wD_dD10 + y/sqrt(x), 2*wD_dD11 + 2*sqrt(x)]]'
        )

    def test_assignment_6(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                    % vardef 'vD' (2D), 'uD' (2D), 'wD' (2D)
                    % keydef index [a-z] (2D)
                    T_{abc} = \vphantom{numeric} ((v_a + u_a)_{,b} - w_{a,b})_{,c}
            """)),
            {'vD', 'uD', 'wD', 'TDDD', 'uD_dD', 'vD_dD', 'wD_dD', 'wD_dDD', 'uD_dDD', 'vD_dDD'}
        )
        self.assertEqual(str(TDDD[0][0][0]),
            'uD_dDD000 + vD_dDD000 - wD_dDD000'
        )

    def test_assignment_7(self):
        Parser.clear_namespace()
        parse(r"""
            % keydef basis [\theta, \phi]
            % vardef -const 'r'
            % vardef 'deltaDD' (2D)
            % keydef index [a-z] (2D)
            % parse g_{\mu\nu} = \delta_{\mu\nu}
            \begin{align*}
                g_{0 0} &= r^{{2}} \\
                g_{1 1} &= r^{{2}} \sin^2(\theta)
            \end{align*}
            % assign -metric 'gDD'
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
        assert_equal(GammaUDD[1][0][1] - GammaUDD[1][1][0], 0, suppress_message=True)
        self.assertEqual(str(GammaUDD[1][0][1]),
            'cos(theta)/sin(theta)'
        )
        assert_equal(RDDDD[0][1][0][1] - (-RDDDD[0][1][1][0]) + (-RDDDD[1][0][0][1]) - RDDDD[1][0][1][0], 0, suppress_message=True)
        self.assertEqual(str(RDDDD[0][1][0][1]),
            'r**2*sin(theta)**2'
        )
        assert_equal(RDD[0][0], 1, suppress_message=True)
        self.assertEqual(str(RDD[1][1]),
            'sin(theta)**2'
        )
        assert_equal(RDD[0][1] - RDD[1][0], 0, suppress_message=True)
        assert_equal(RDD[0][1], 0, suppress_message=True)
        self.assertEqual(str(R),
            '2/r**2'
        )

    def test_assignment_8(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % vardef -metric 'gDD' (4D)
                \gamma_{ij} = g_{ij}
            """)),
            {'gUU', 'gdet', 'epsilonUUUU', 'gDD', 'gammaDD'}
        )
        self.assertEqual(str(gammaDD),
            '[[gDD11, gDD12, gDD13], [gDD12, gDD22, gDD23], [gDD13, gDD23, gDD33]]'
        )
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % vardef -nosym 'TUU' (3D)
                % vardef 'vD' (2D)
                % keydef index i (2D)
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
                % vardef -nosym 'hUD' (4D)
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
                % vardef -metric 'gUU' (3D)
                % vardef 'vD' (3D)
                v^\mu = g^{\mu\nu} v_\nu
            """)),
            {'gDD', 'gdet', 'epsilonDDD', 'gUU', 'vD', 'vU'}
        )
        self.assertEqual(str(vU),
            '[gUU00*vD0 + gUU01*vD1 + gUU02*vD2, gUU01*vD0 + gUU11*vD1 + gUU12*vD2, gUU02*vD0 + gUU12*vD1 + gUU22*vD2]'
        )

    def test_example_3(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % vardef 'epsilonDDD' (3D)
                % vardef 'vU' (3D), 'wU' (3D)
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
                % vardef -anti01 'FUU' (4D)
                % vardef -metric 'gDD' (4D)
                % vardef -const 'k'
                J^\mu = (4\pi k)^{-1} \vphantom{numeric} \nabla_\nu F^{\mu\nu}
            """)),
            {'FUU', 'gUU', 'gdet', 'epsilonUUUU', 'gDD', 'FUU_dD', 'gDD_dD', 'GammaUDD', 'FUU_cdD', 'JU'}
        )

    def test_example_4_2(self):
        Parser.clear_namespace()
        self.assertEqual(
            set(parse(r"""
                % vardef -anti01 'FUU' (4D)
                % vardef -metric 'ghatDD' (4D)
                % vardef -const 'k'
                J^\mu = (4\pi k)^{-1} \vphantom{numeric} \hat{\nabla}_\nu F^{\mu\nu}
            """)),
            {'FUU', 'ghatUU', 'ghatdet', 'epsilonUUUU',  'ghatDD', 'FUU_dD', 'ghatDD_dD', 'GammahatUDD', 'FUU_cdhatD', 'JU'}
        )

    def test_example_5_1(self):
        Parser.clear_namespace()
        parse(r"""
            % vardef 'deltaDD' (4D)
            % vardef -const 'G', 'M'
            % parse g_{\mu\nu} = \delta_{\mu\nu}
            \begin{align}
                g_{0 0} &= -\left(1 - \frac{2GM}{r}\right) \\
                g_{1 1} &=  \left(1 - \frac{2GM}{r}\right)^{-1} \\
                g_{2 2} &= r^{{2}} \\
                g_{3 3} &= r^{{2}} \sin^2\theta
            \end{align}
            % assign -metric 'gDD'
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
            % keydef basis [t, r, \theta, \phi]
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
        assert_equal(GammaUDD[0][0][1] - GammaUDD[0][1][0], 0, suppress_message=True)
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
        assert_equal(GammaUDD[2][1][2] - GammaUDD[2][2][1], 0, suppress_message=True)
        self.assertEqual(str(GammaUDD[2][1][2]),
            '1/r'
        )
        self.assertEqual(str(GammaUDD[2][3][3]),
            '-sin(theta)*cos(theta)'
        )
        assert_equal(GammaUDD[2][1][3] - GammaUDD[2][3][1], 0, suppress_message=True)
        self.assertEqual(str(GammaUDD[3][1][3]),
            '1/r'
        )
        assert_equal(GammaUDD[3][2][3] - GammaUDD[3][3][2], 0, suppress_message=True)
        self.assertEqual(str(GammaUDD[3][2][3]),
            'cos(theta)/sin(theta)'
        )
        self.assertEqual(str(simplify(K)),
            '48*G**2*M**2/r**6'
        )
        assert_equal(R, 0, suppress_message=True)
        for i in range(3):
            for j in range(3):
                assert_equal(GDD[i][j], 0, suppress_message=True)

    @staticmethod
    def test_example_6_1():
        parse(r"""
            % keydef basis [r, \theta, \phi]
            \begin{align}
                \gamma_{ij} &= g_{ij} \\
                % assign -metric 'gammaDD'
                \beta_i &= g_{0 i} \\
                \alpha &= \sqrt{\gamma^{ij}\beta_i\beta_j - g_{0 0}} \\
                K_{ij} &= \frac{1}{2\alpha}\left(\nabla_i \beta_j + \nabla_j \beta_i\right) \\
                K &= \gamma^{ij} K_{ij}
            \end{align}
        """)
        for i in range(3):
            for j in range(3):
                assert_equal(KDD[i][j], 0, suppress_message=True)

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
        # assert_equal(E, 0, suppress_message=True)
        self.assertEqual(simplify(E), 0)
        for i in range(3):
            assert_equal(pD[i], 0, suppress_message=True)

    @staticmethod
    def test_metric_inverse():
        for DIM in range(2, 5):
            Parser.clear_namespace()
            parse(r"""
                % vardef -metric 'gDD' ({DIM}D)
                \delta^a_c = g^{{ab}} g_{{bc}}
            """.format(DIM=DIM))
            for i in range(DIM):
                for j in range(DIM):
                    value = 1 if i == j else 0
                    assert_equal(deltaUD[i][j], value, suppress_message=True)
        for DIM in range(2, 5):
            Parser.clear_namespace()
            parse(r"""
                % vardef -metric 'gUU' ({DIM}D)
                \delta^a_c = g^{{ab}} g_{{bc}}
            """.format(DIM=DIM))
            for i in range(DIM):
                for j in range(DIM):
                    value = 1 if i == j else 0
                    assert_equal(deltaUD[i][j], value, suppress_message=True)

    @staticmethod
    def test_example_BSSN():
        import NRPy_param_funcs as par, reference_metric as rfm
        import BSSN.BSSN_RHSs as Brhs, BSSN.BSSN_quantities as Bq
        import BSSN.BSSN_gauge_RHSs as Bgrhs
        Parser.clear_namespace()
        parse(r"""
            \begin{align}
                % keydef basis [x, y, z]
                % ignore "\\%", "\qquad"

                % vardef 'deltaDD'
                % parse \hat{\gamma}_{ij} = \delta_{ij}
                % assign -symbolic <H> -metric 'gammahatDD'
                % vardef -numeric -sym01 'hDD'
                % parse \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}
                % assign -numeric -metric 'gammabarDD'

                % vardef -numeric 'vetU'
                % srepl "\beta" -> "\text{vet}"
                %% upwind pattern inside Lie derivative expansion
                % srepl -global "\text{vet}^<1> \partial_<1>" -> "\text{vet}^<1> \vphantom{upwind} \partial_<1>"
                %% substitute tensor identity (see appropriate BSSN notebook)
                % srepl -global "\bar{D}_k \text{vet}^k" -> "(\partial_k \text{vet}^k + \frac{\partial_k \text{gammahatdet} \text{vet}^k}{2 \text{gammahatdet}})"

                % vardef -numeric 'alpha'
                % vardef -numeric -sym01 'aDD'
                % srepl "\bar{A}" -> "\text{a}"
                % parse \bar{A}^i_j = \bar{\gamma}^{ik} \bar{A}_{kj}
                % srepl "\partial_t \bar{\gamma}" -> "\text{h_rhs}"
                \partial_t \bar{\gamma}_{ij} &= \mathcal{L}_\beta \bar{\gamma}_{ij}
                    + \frac{2}{3} \bar{\gamma}_{ij} \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right)
                    - 2 \alpha \bar{A}_{ij} \\

                % vardef -numeric 'cf', 'trK'
                % srepl "K" -> "\text{trK}"
                %% replace 'phi' with conformal factor cf = W = e^{{-2\phi}}
                % srepl "e^{-4\phi}" -> "\text{cf}^{{2}}"
                % srepl "\partial_t \phi = <1..> \\" -> "\text{cf_rhs} = -2 \text{cf} (<1..>) \\"
                % srepl -global "\partial_<1> \phi"       -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"
                % srepl -global "\partial_<1> \text{phi}" -> "\partial_<1> \text{cf} \frac{-1}{2 \text{cf}}"
                \partial_t \phi &= \mathcal{L}_\beta \phi
                    + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

                % parse \bar{A}^{ij} = \bar{\gamma}^{jk} \bar{A}^i_k
                % srepl "\partial_t \text{trK}" -> "\text{trK_rhs}"
                \partial_t K &= \mathcal{L}_\beta K
                    + \frac{1}{3} \alpha K^{{2}}
                    + \alpha \bar{A}_{ij} \bar{A}^{ij}
                    - e^{-4\phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi\right) \\

                % vardef -numeric 'lambdaU'
                % srepl "\bar{\Lambda}" -> "\text{lambda}"
                % parse \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
                % parse \Delta_{ijk}  = \bar{\gamma}_{il} \Delta^l_{jk}
                % parse \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
                % srepl "\partial_t \text{lambda}" -> "\text{Lambdabar_rhs}"
                \partial_t \bar{\Lambda}^i &= \mathcal{L}_\beta \bar{\Lambda}^i + \bar{\gamma}^{jk} \hat{D}_j \hat{D}_k \beta^i
                    + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
                    &\qquad- 2 \bar{A}^{ij} \left(\partial_j \alpha - 6 \alpha \partial_j \phi\right)
                    + 2 \alpha \bar{A}^{jk} \Delta^i_{jk} - \frac{4}{3} \alpha \bar{\gamma}^{ij} \partial_j K \\

                % vardef -numeric -sym01 'RbarDD'
                X_{ij} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi
                    + 2 \bar{D}_i \alpha \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi
                    - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{ij} \\
                \hat{X}_{ij} &= X_{ij} - \frac{1}{3} \bar{\gamma}_{ij} \bar{\gamma}^{kl} X_{kl}
                % srepl "\partial_t \text{a}" -> "\text{a_rhs}"
                \partial_t \bar{A}_{ij} &= \mathcal{L}_\beta \bar{A}_{ij}
                    - \frac{2}{3} \bar{A}_{ij} \bar{D}_k \beta^k
                    - 2 \alpha \bar{A}_{ik} \bar{A}^k_j
                    + \alpha \bar{A}_{ij} K
                    + e^{-4\phi} \hat{X}_{ij} \\

                % srepl "\partial_t \alpha" -> "\text{alpha_rhs}"
                \partial_t \alpha &= \mathcal{L}_\beta \alpha - 2 \alpha K \\

                % vardef -numeric 'eta', 'betU'
                % srepl "B" -> "\text{bet}"
                % srepl "\partial_t \text{vet}" -> "\text{vet_rhs}"
                \partial_t \beta^i &= \left[\beta^j \vphantom{upwind} \bar{D}_j \beta^i\right] + B^{i} \\
                % srepl "\partial_t \text{bet}" -> "\text{bet_rhs}"
                \partial_t B^i &= \left[\beta^j \vphantom{upwind} \bar{D}_j B^i\right] + \frac{3}{4} \left(\partial_t \bar{\Lambda}^{i} - \left[\beta^j \vphantom{upwind} \bar{D}_j \bar{\Lambda}^{i}\right]\right) - \eta B^{i} \\

                % parse \bar{R} = \bar{\gamma}^{ij} \bar{R}_{ij}
                % srepl "\bar{D}^2" -> "\bar{D}^i \bar{D}_i", "\mathcal{<1>}" -> "<1>"
                \mathcal{H} &= \frac{2}{3} K^{{2}} - \bar{A}_{ij} \bar{A}^{ij} + e^{-4\phi} \left(\bar{R} - 8 \bar{D}^i \phi \bar{D}_i \phi - 8 \bar{D}^2 \phi\right) \\
                \mathcal{M}^i &= e^{-4\phi} \left(\hat{D}_j \bar{A}^{ij} + 6 \bar{A}^{ij}\partial_j \phi - \frac{2}{3} \bar{\gamma}^{ij} \partial_j K + \bar{A}^{jk} \Delta^i_{jk} + \bar{A}^{ik} \Delta^j_{jk}\right) \\

                \bar{R}_{ij} &= -\frac{1}{2} \bar{\gamma}^{kl} \hat{D}_k \hat{D}_l \bar{\gamma}_{ij}
                    + \frac{1}{2} \left(\bar{\gamma}_{ki} \hat{D}_j \bar{\Lambda}^k + \bar{\gamma}_{kj} \hat{D}_i \bar{\Lambda}^k\right)
                    + \frac{1}{2} \Delta^k \left(\Delta_{ijk} + \Delta_{jik}\right) \\%
                    &\qquad+ \bar{\gamma}^{kl} \left(\Delta^m_{ki} \Delta_{jml} + \Delta^m_{kj} \Delta_{iml} + \Delta^m_{ik} \Delta_{mjl}\right)
            \end{align}
        """)
        par.set_parval_from_str('reference_metric::CoordSystem', 'Cartesian')
        par.set_parval_from_str('BSSN.BSSN_quantities::LeaveRicciSymbolic', 'True')
        rfm.reference_metric(); Brhs.BSSN_RHSs(); Bgrhs.BSSN_gauge_RHSs()
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
                      'RbarDD': RbarDD},
                     {'h_rhsDD': Brhs.h_rhsDD,
                      'cf_rhs': Brhs.cf_rhs,
                      'trK_rhs': Brhs.trK_rhs,
                      'Lambdabar_rhsU': Brhs.Lambdabar_rhsU,
                      'a_rhsDD': Brhs.a_rhsDD,
                      'alpha_rhs': Bgrhs.alpha_rhs,
                      'vet_rhsU': Bgrhs.vet_rhsU,
                      'bet_rhsU': Bgrhs.bet_rhsU,
                      'RbarDD': Bq.RbarDD},
                    suppress_message=True)

if __name__ == '__main__':
    suite = unittest.TestSuite()
    for i in range(6):
        suite.addTest(TestParser('test_expression_' + str(i + 1)))
    suite.addTest(TestParser('test_srepl_macro'))
    for i in range(8):
        suite.addTest(TestParser('test_assignment_' + str(i + 1)))
    for i in range(6):
        if i > 2:
            suite.addTest(TestParser('test_example_' + str(i + 1) + '_1'))
            suite.addTest(TestParser('test_example_' + str(i + 1) + '_2'))
        else:
            suite.addTest(TestParser('test_example_' + str(i + 1)))
    suite.addTest(TestParser('test_metric_inverse'))
    suite.addTest(TestParser('test_example_BSSN'))
    result = unittest.TextTestRunner().run(suite)
    sys.exit(not result.wasSuccessful())
