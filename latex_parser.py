""" Convert LaTeX Sentence to SymPy Expression """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

# pylint: disable=attribute-defined-outside-init,exec-used
from sympy import Function, Derivative, Symbol, Integer, Rational, Float, Pow
from sympy import sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
from sympy import pi, exp, log, sqrt, expand, diff
from collections import OrderedDict
from functional import chain, uniquify
from expr_tree import ExprTree
from inspect import currentframe
import indexedexp as ixp
import re, warnings

sympy_env = {'sin': sin, 'cos': cos, 'tan': tan, 'sinh': sinh, 'cosh': cosh, 'tanh': tanh,
    'asin': asin, 'acos': acos, 'atan': atan, 'asinh': asinh, 'acosh': acosh, 'atanh': atanh,
    'pi': pi, 'exp': exp, 'log': log, 'sqrt': sqrt}

class Lexer:
    """ LaTeX Lexer

        The following class will tokenize a LaTeX sentence for parsing.
    """

    def __init__(self, namespace=None):
        if namespace is None: namespace = {}
        greek_symbol = '|'.join(letter for letter in (r'\\[aA]lpha', r'\\[bB]eta', r'\\[gG]amma', r'\\[dD]elta',
            r'\\[eE]psilon', r'\\[zZ]eta', r'\\[eE]ta', r'\\[tT]heta', r'\\[iI]ota', r'\\[kK]appa', r'\\[lL]ambda',
            r'\\[mM]u', r'\\[nN]u', r'\\[xX]i', r'\\[oO]mikron', r'\\[pP]i', r'\\[Rr]ho', r'\\[sS]igma', r'\\[tT]au',
            r'\\[uU]psilon', r'\\[pP]hi', r'\\[cC]hi', r'\\[pP]si', r'\\[oO]mega'))
        symmetry = r'nosym|(?:sym|anti)[0-9]+(?:_(?:sym|anti)[0-9]+)*|metric|permutation|kronecker'
        # define a regex pattern for every token, create a named capture group for
        # every pattern, join together the resulting pattern list using a pipe symbol
        # for regex alternation, and compile the generated regular expression
        self.regex = re.compile('|'.join(['(?P<%s>%s)' % pattern for pattern in
            [ ('SPACE_DELIM',    r'(?:\s|\\,|\{\})+|\&'),
              ('RATIONAL',       r'[0-9]+\/[1-9]+|\\frac{[0-9]+}{[1-9]+}'),
              ('DECIMAL',        r'[0-9]+\.[0-9]+'),
              ('INTEGER',        r'\-?[0-9]+'),
              ('PI',             r'\\pi'),
              ('EULER',          r'e'),
              ('PLUS',           r'\+'),
              ('MINUS',          r'\-'),
              ('DIVIDE',         r'\/'),
              ('EQUAL',          r'\='),
              ('CARET',          r'\^'),
              ('UNDERSCORE',     r'\_'),
              ('COLON',          r'\:'),
              ('COMMA',          r'\,'),
              ('COMMENT',        r'\%'),
              ('LEFT_PAREN',     r'\('),
              ('RIGHT_PAREN',    r'\)'),
              ('LEFT_BRACE',     r'\{'),
              ('RIGHT_BRACE',    r'\}'),
              ('LEFT_BRACKET',   r'\['),
              ('RIGHT_BRACKET',  r'\]'),
              ('BIGL_DELIM',     r'\\[bB]igl'),
              ('BIGR_DELIM',     r'\\[bB]igr'),
              ('LEFT_DELIM',     r'\\left'),
              ('RIGHT_DELIM',    r'\\right'),
              ('LINE_BREAK',     r'(?:\;|\\\\|\\cr)'),
              ('BEGIN_ALIGN',    r'\\begin{align\*?}'),
              ('END_ALIGN',      r'\\end{align\*?}'),
              ('SQRT_CMD',       r'\\sqrt'),
              ('FRAC_CMD',       r'\\frac'),
              ('TRIG_CMD',       r'\\sinh|\\cosh|\\tanh|\\sin|\\cos|\\tan'),
              ('NLOG_CMD',       r'\\ln|\\log'),
              ('PARDRV_OP',      r'\\partial'),
              ('COVDRV_OP',      r'\\nabla'),
              ('DIACRITIC',      r'\\hat|\\tilde|\\bar'),
              ('SYMMETRY',       symmetry),
              ('MATHOP',         r'\\mathop'),
              ('SYMBOL',         r'[a-zA-Z]|' + greek_symbol),
              ('COMMAND',        r'\\[a-z]+')]]))

    def initialize(self, sentence, position=0):
        """ Initialize Lexer

            :arg: sentence (raw string)
            :arg: position [default: 0]
        """
        self.sentence = sentence
        self.token    = None
        self.lexeme   = None
        self.marker   = None
        self.index    = position

    def tokenize(self):
        """ Tokenize Sentence

            :return: token iterator
        """
        while self.index < len(self.sentence):
            token = self.regex.match(self.sentence, self.index)
            if token is None:
                raise ParseError('unexpected \'%s\' at position %d' %
                    (self.sentence[self.index], self.index), self.sentence, self.index)
            self.index = token.end()
            if token.lastgroup not in ('SPACE_DELIM', 'BIGL_DELIM', 'BIGR_DELIM', 'LEFT_DELIM', 'RIGHT_DELIM'):
                self.lexeme = token.group()
                yield token.lastgroup

    def lex(self):
        """ Retrieve Next Token

            :return: next token
        """
        try:
            self.token = next(self.tokenize())
        except StopIteration:
            self.token  = None
            self.lexeme = ''
        return self.token

    def mark(self):
        """ Mark Iterator Position

            :return: previous position
        """
        self.marker = self.index - len(self.lexeme)
        return self.marker

    def reset(self):
        """ Reset Token Iterator """
        if not self.sentence:
            raise RuntimeError('cannot reset uninitialized lexer')
        self.initialize(self.sentence, self.marker)
        self.lex()

class Parser:
    """ LaTeX Parser

        The following class will parse a tokenized LaTeX sentence.

        LaTeX Extended BNF Grammar:
        <ROOT>          -> <STRUCTURE> { <LINE_BREAK> <STRUCTURE> }*
        <STRUCTURE>     -> <CONFIG> | <ENVIRONMENT> | <ASSIGNMENT>
        <ENVIRONMENT>   -> <BEGIN_ALIGN> <ASSIGNMENT> { <LINE_BREAK> <ASSIGNMENT> }* <END_ALIGN>
        <ASSIGNMENT>    -> ( <TENSOR> | <COVDRV> ) = <EXPRESSION>
        <EXPRESSION>    -> <TERM> { ( '+' | '-' ) <TERM> }*
        <TERM>          -> <FACTOR> { [ '/' ] <FACTOR> }*
        <FACTOR>        -> ( <BASE> | <EULER> ) { '^' <EXPONENT> }*
        <BASE>          -> [ '-' ] ( <ATOM> | '(' <EXPRESSION> ')' )
        <EXPONENT>      -> <BASE> | '{' <BASE> '}' | '{{' <BASE> '}}'
        <ATOM>          -> <NUMBER> | <TENSOR> | <COMMAND> | <OPERATOR>
        <NUMBER>        -> <RATIONAL> | <DECIMAL> | <INTEGER> | <PI>
        <COMMAND>       -> <SQRT> | <FRAC> | <NLOG> | <TRIG>
        <OPERATOR>      -> <PARDRV> | <COVDRV>
        <SQRT>          -> <SQRT_CMD> [ '[' <INTEGER> ']' ] '{' <EXPRESSION> '}'
        <FRAC>          -> <FRAC_CMD> '{' <EXPRESSION> '}' '{' <EXPRESSION> '}'
        <NLOG>          -> <NLOG_CMD> [ '_' ( <INTEGER> | { <INTEGER> } ) ] ( <SYMBOL> | <INTEGER> | '(' <EXPRESSION> ')' )
        <TRIG>          -> <TRIG_CMD> [ '^' ( <INTEGER> | { <INTEGER> } ) ] ( <SYMBOL> | <INTEGER> | '(' <EXPRESSION> ')' )
        <PARDRV>        -> { <PARDRV_OP> '_' <SYMBOL> }+ ( <TENSOR> | '(' <EXPRESSION> ')' )
        <COVDRV>        -> { ( <COVDRV_OP> | <DIACRITIC> '{' <COVDRV_OP> '}' ) ( '^' | '_' ) <SYMBOL> }+ ( <TENSOR> | '(' <EXPRESSION> ')' )
        <CONFIG>        -> '%' { <VARIABLE> }+ '[' <INTEGER> ']' ':' <SYMMETRY> { ',' { <VARIABLE> }+ '[' <INTEGER> ']' ':' <SYMMETRY> }*
        <VARIABLE>      -> <SYMBOL> | <DIACRITIC> '{' <SYMBOL> '}' | <MATHOP> '{' <SYMBOL> { <SYMBOL> | <INTEGER> | <UNDERSCORE> }* '}'
        <TENSOR>        -> <VARIABLE> [ ( '_' <LOWER_INDEX> ) | ( '^' <UPPER_INDEX> [ '_' <LOWER_INDEX> ] ) ]
        <LOWER_INDEX>   -> <SYMBOL> | <INTEGER> | '{' { <SYMBOL> | <INTEGER> }* [ ',' { <SYMBOL> }+ ] '}'
        <UPPER_INDEX>   -> <SYMBOL> | <INTEGER> | '{' { <SYMBOL> | <INTEGER> }+ '}'
    """

    namespace = OrderedDict()

    def __init__(self):
        self.lexer     = Lexer()
        self.namespace = OrderedDict()
        self.unevaled  = OrderedDict()
        self.dimension = None

    def parse(self, sentence, expression=False):
        """ Parse LaTeX Sentence

            :arg:    latex sentence (raw string)
            :arg:    expression mode [default: disabled]
            :return: namespace or expression
        """
        stack = []; i_1 = i_2 = i_3 = 0
        for i, lexeme in enumerate(sentence):
            # convert derivative to operator notation for parenthetical expression
            if   lexeme == '(': stack.append(i)
            elif lexeme == ')': i_1, i_2 = stack.pop(), i + 1
            elif lexeme == ',' and sentence[i - 1] == '{':
                i_3 = sentence.find('}', i) + 1
                subexpr, indices = sentence[i_1:i_2], sentence[i_2:i_3][3:-1]
                operator = ''.join('\\partial_' + index for index in indices)
                sentence = sentence.replace(sentence[i_1:i_3], operator + ' ' + subexpr)
            elif lexeme == ';' and sentence[i - 1] == '{':
                i_3 = sentence.find('}', i) + 1
                subexpr, indices = sentence[i_1:i_2], sentence[i_2:i_3][3:-1]
                operator = ''.join('\\nabla_' + index for index in indices)
                sentence = sentence.replace(sentence[i_1:i_3], operator + ' ' + subexpr)
        self.lexer.initialize(sentence)
        self.lexer.lex()
        if expression:
            return self._expression()
        root = self._root()
        for key in chain(*root):
            if key in sympy_env:
                sympy_env.pop(key)
        return root

    # <ROOT> -> <STRUCTURE> { <LINE_BREAK> <STRUCTURE> }*
    def _root(self):
        self._structure()
        while self.accept('LINE_BREAK'):
            self._structure()
        return self.namespace, self.unevaled

    # <STRUCTURE> -> <CONFIG> | <ENVIRONMENT> | <ASSIGNMENT>
    def _structure(self):
        if self.peek('COMMENT'):
            self._config()
        elif self.peek('BEGIN_ALIGN'):
            self._environment()
        else: self._assignment()

    # <ENVIRONMENT> -> <BEGIN_ALIGN> <ASSIGNMENT> { <LINE_BREAK> <ASSIGNMENT> }* <END_ALIGN>
    def _environment(self):
        self.expect('BEGIN_ALIGN')
        self._assignment()
        while self.accept('LINE_BREAK'):
            self._assignment()
        self.expect('END_ALIGN')

    # <ASSIGNMENT> -> ( <TENSOR> | <COVDRV> ) = <EXPRESSION>
    def _assignment(self):
        variable = self._covdrv('LHS') if self.peek('COVDRV_OP') \
              else self._tensor()
        indexed  = variable.func == Function('Tensor')
        variable = Tensor.array_format(variable) if indexed else str(variable)
        self.expect('EQUAL')
        expr = self._expression()
        # distribute over every parenthetical expression
        tree = ExprTree(expand(expr))
        if not indexed:
            for subtree in tree.preorder():
                subexpr = subtree.expr
                if subexpr.func == Function('Tensor'):
                    indexed = True
        if indexed:
            # iterate through every tensor function
            expr = str(tree.root.expr)
            for subtree in tree.preorder():
                subexpr = subtree.expr
                if subexpr.func == Function('Tensor'):
                    # replace function notation with array notation
                    func = Tensor.array_format(subexpr)
                    expr = expr.replace(str(subexpr), func)
            # perform implied summation on tensorial equation
            expanded = _summation((variable, expr), self.dimension)
            # update unevalued with expanded equation
            self.unevaled.update(expanded)
        else: self.namespace.update({variable: expr})

    # <EXPRESSION> -> <TERM> { ( '+' | '-' ) <TERM> }*
    def _expression(self):
        expr = self._term()
        while self.peek('PLUS') or self.peek('MINUS'):
            if self.accept('PLUS'):
                expr += self._term()
            elif self.accept('MINUS'):
                expr -= self._term()
        return expr

    # <TERM> -> <FACTOR> { [ '/' ] <FACTOR> }*
    def _term(self):
        expr = self._factor()
        while any(self.peek(token) for token in ('LEFT_PAREN', 'PARDRV_OP', 'COVDRV_OP',
                'SYMBOL', 'RATIONAL', 'DECIMAL', 'INTEGER', 'PI', 'EULER', 'DIACRITIC',
                'DIVIDE', 'COMMAND', 'SQRT_CMD', 'FRAC_CMD', 'TRIG_CMD', 'NLOG_CMD')):
            if self.accept('DIVIDE'):
                expr /= self._factor()
            else: expr *= self._factor()
        return expr

    # <FACTOR> -> ( <BASE> | <EULER> ) { '^' <EXPONENT> }*
    def _factor(self):
        stack = ['e'] if self.accept('EULER') else [self._base()]
        while self.accept('CARET'):
            stack.append(self._exponent())
        expr = stack.pop()
        for subexpr in reversed(stack):
            expr = exp(expr) if subexpr == 'e' else subexpr ** expr
        return expr

    # <BASE> -> [ '-' ] ( <ATOM> | '(' <EXPRESSION> ')' )
    def _base(self):
        sign = -1 if self.accept('MINUS') else 1
        if self.accept('LEFT_PAREN'):
            expr = sign * self._expression()
            self.expect('RIGHT_PAREN')
            return expr
        return sign * self._atom()

    # <EXPONENT> -> <BASE> | '{' <BASE> '}' | '{{' <BASE> '}}'
    def _exponent(self):
        if self.accept('LEFT_BRACE'):
            if self.accept('LEFT_BRACE'):
                base = self._base()
                self.expect('RIGHT_BRACE')
            else: base = self._base()
            self.expect('RIGHT_BRACE')
            return base
        return self._base()

    # <ATOM> -> <NUMBER> | <TENSOR> | <COMMAND> | <OPERATOR>
    def _atom(self):
        if any(self.peek(token) for token in
                ('RATIONAL', 'DECIMAL', 'INTEGER', 'PI')):
            return self._number()
        if any(self.peek(token) for token in
                ('SYMBOL', 'DIACRITIC', 'MATHOP')):
            function = self._tensor()
            if function.func == Function('Tensor'):
                symbol = str(function.args[0])
                # reserved keyword for christoffel symbol
                if 'Gamma' in symbol and symbol not in self.namespace:
                    self._instantiate_christoffel(function)
            return function
        if any(self.peek(token) for token in
                ('COMMAND', 'SQRT_CMD', 'FRAC_CMD', 'NLOG_CMD', 'TRIG_CMD')):
            return self._command()
        if any(self.peek(token) for token in
                ('PARDRV_OP', 'COVDRV_OP')):
            return self._operator()
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    # <NUMBER> -> <RATIONAL> | <DECIMAL> | <INTEGER> | <PI>
    def _number(self):
        number = self.lexer.lexeme
        if self.accept('RATIONAL'):
            rational = re.match(r'([1-9][0-9]*)\/([1-9][0-9]*)', number)
            if not rational:
                rational = re.match(r'\\frac{([0-9]+)}{([1-9]+)}', number)
            return Rational(rational.group(1), rational.group(2))
        if self.accept('DECIMAL'):
            return Float(number)
        if self.accept('INTEGER'):
            return Integer(number)
        if self.accept('PI'):
            return pi
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    # <COMMAND> -> <SQRT> | <FRAC> | <NLOG> | <TRIG>
    def _command(self):
        command = self.lexer.lexeme
        if self.peek('SQRT_CMD'):
            return self._sqrt()
        if self.peek('FRAC_CMD'):
            return self._frac()
        if self.peek('NLOG_CMD'):
            return self._nlog()
        if self.peek('TRIG_CMD'):
            return self._trig()
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unsupported command \'%s\' at position %d' %
            (command, position), sentence, position)

    # <OPERATOR> -> <PARDRV> | <COVDRV>
    def _operator(self):
        operator = self.lexer.lexeme
        if self.peek('PARDRV_OP'):
            return self._pardrv()
        if self.peek('COVDRV_OP'):
            return self._covdrv('RHS')
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unsupported operator \'%s\' at position %d' %
            (operator, position), sentence, position)

    # <SQRT> -> <SQRT_CMD> [ '[' <INTEGER> ']' ] '{' <EXPRESSION> '}'
    def _sqrt(self):
        self.expect('SQRT_CMD')
        if self.accept('LEFT_BRACKET'):
            integer = self.lexer.lexeme
            self.expect('INTEGER')
            root = Rational(1, integer)
            self.expect('RIGHT_BRACKET')
        else: root = Rational(1, 2)
        self.expect('LEFT_BRACE')
        expr = self._expression()
        self.expect('RIGHT_BRACE')
        if root == Rational(1, 2):
            return sqrt(expr)
        return Pow(expr, root)

    # <FRAC> -> <FRAC_CMD> '{' <EXPRESSION> '}' '{' <EXPRESSION> '}'
    def _frac(self):
        self.expect('FRAC_CMD')
        self.expect('LEFT_BRACE')
        numerator = self._expression()
        self.expect('RIGHT_BRACE')
        self.expect('LEFT_BRACE')
        denominator = self._expression()
        self.expect('RIGHT_BRACE')
        return numerator / denominator

    # <NLOG> -> <NLOG_CMD> [ '_' <INTEGER> | { <INTEGER> } ] ( <SYMBOL> | <INTEGER> | '(' <EXPRESSION> ')' )
    def _nlog(self):
        func = self.lexer.lexeme[1:]
        self.expect('NLOG_CMD')
        if func == 'log':
            if self.accept('UNDERSCORE'):
                if self.accept('LEFT_BRACE'):
                    base = self.lexer.lexeme
                    self.expect('INTEGER')
                    self.expect('RIGHT_BRACE')
                else:
                    base = self.lexer.lexeme
                    self.expect('INTEGER')
                base = int(base)
            else: base = 10
        if self.peek('SYMBOL'):
            expr = self.lexer.lexeme
            self.expect('SYMBOL')
        elif self.peek('INTEGER'):
            expr = self.lexer.lexeme
            self.expect('INTEGER')
        elif self.accept('LEFT_PAREN'):
            expr = self._expression()
            self.expect('RIGHT_PAREN')
        else:
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('unexpected \'%s\' at position %d' %
                (sentence[position], position), sentence, position)
        if func == 'ln': return log(expr)
        return log(expr, base)

    # <TRIG> -> <TRIG_CMD> [ '^' <INTEGER> | { <INTEGER> } ] ( <SYMBOL> | <INTEGER> | '(' <EXPRESSION> ')' )
    def _trig(self):
        func = self.lexer.lexeme[1:]
        self.expect('TRIG_CMD')
        if self.accept('CARET'):
            if self.accept('LEFT_BRACE'):
                exponent = self.lexer.lexeme
                self.expect('INTEGER')
                self.expect('RIGHT_BRACE')
            else:
                exponent = self.lexer.lexeme
                self.expect('INTEGER')
            exponent = int(exponent)
        else: exponent = 1
        if   func == 'cosh': trig = acosh if exponent == -1 else cosh
        elif func == 'sinh': trig = asinh if exponent == -1 else sinh
        elif func == 'tanh': trig = atanh if exponent == -1 else tanh
        elif func == 'cos':  trig = acos  if exponent == -1 else cos
        elif func == 'sin':  trig = asin  if exponent == -1 else sin
        elif func == 'tan':  trig = atan  if exponent == -1 else tan
        if self.peek('SYMBOL'):
            expr = self.lexer.lexeme
            self.expect('SYMBOL')
        elif self.peek('INTEGER'):
            expr = self.lexer.lexeme
            self.expect('INTEGER')
        elif self.accept('LEFT_PAREN'):
            expr = self._expression()
            self.expect('RIGHT_PAREN')
        else:
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('unexpected \'%s\' at position %d' %
                (sentence[position], position), sentence, position)
        if exponent == -1: return trig(expr)
        return trig(expr) ** exponent

    # <PARDRV> -> { <PARDRV_OP> '_' <SYMBOL> }+ ( <TENSOR> | '(' <EXPRESSION> ')' )
    def _pardrv(self):
        indices = []
        while self.accept('PARDRV_OP'):
            self.expect('UNDERSCORE')
            index = self._strip(self.lexer.lexeme)
            self.expect('SYMBOL')
            indices.append(Symbol(index))
        if self.accept('LEFT_PAREN'):
            tree = ExprTree(self._expression())
            self.expect('RIGHT_PAREN')
            for subtree in tree.preorder():
                subexpr = subtree.expr
                if subexpr.func == Function('Tensor'):
                    # insert temporary symbol '_x' for symbolic differentiation
                    subtree.expr = Function('Tensor')(Symbol('_x'), *subexpr.args)
                    tree.build(subtree)
            expr = tree.reconstruct()
            # differentiate the expression, including product rule expansion
            tree = ExprTree(diff(expr, Symbol('_x')))
            for subtree in tree.preorder():
                subexpr = subtree.expr
                if subexpr.func == Derivative:
                    # instantiate each derivative and update namespace
                    order, function = len(indices), subexpr.args[0]
                    indices = list(function.args[2:]) + indices
                    tensor = Tensor(function, self.dimension)
                    subtree.expr = self._instantiate_derivative(tensor, order, indices)
                    tree.build(subtree)
                elif subexpr.func == Function('Tensor'):
                    # remove temporary symbol '_x' from tensor function
                    subtree.expr = Function('Tensor')(*subexpr.args[1:])
                    tree.build(subtree)
            return tree.reconstruct()
        # instantiate partial derivative and update namespace
        order, function = len(indices), self._tensor()
        indices = list(function.args[1:]) + indices
        tensor = Tensor(function, self.dimension)
        return self._instantiate_derivative(tensor, order, indices)

    # <COVDRV> -> { <COVDRV_OP> ( '^' | '_' ) <SYMBOL> }+ ( <TENSOR> | '(' <EXPRESSION> ')' )
    def _covdrv(self, location):
        indices, config, equation = [], [], ['', ' = ', '', '']
        while self.accept('COVDRV_OP'):
            equation[0] += '\\nabla'
            equation[3] += '\\nabla'
            if self.accept('CARET'):
                index = self.lexer.lexeme
                equation[0] += '^' + index
                bound_index = next(x for x in (chr(97 + n) for n in range(26)) if x != index)
                equation[2] += 'g^{%s %s} ' % (index, bound_index)
                config.append(equation[2] + '[%d]: metric' % self.dimension)
                equation[3] += '_' + bound_index
                index = self._strip(index)
                self.expect('SYMBOL')
                indices.append((Symbol(index), 'U'))
            elif self.accept('UNDERSCORE'):
                index = self.lexer.lexeme
                equation[0] += '_' + index
                equation[3] += '_' + index
                index = self._strip(index)
                self.expect('SYMBOL')
                indices.append((Symbol(index), 'D'))
            else:
                sentence, position = self.lexer.sentence, self.lexer.mark()
                raise ParseError('unexpected \'%s\' at position %d' %
                    (sentence[position], position), sentence, position)
        marker_1 = self.lexer.mark()
        order, function = len(indices), self._tensor()
        marker_2 = self.lexer.index
        equation[0] += ' ' + self.lexer.sentence[marker_1:marker_2]
        equation[3] += ' ' + self.lexer.sentence[marker_1:marker_2]
        tensor = Tensor(function, self.dimension)
        if location == 'RHS':
            if equation[2]:
                sentence, position = self.lexer.sentence, self.lexer.mark()
                if config: config = '% ' + ';\n'.join(config) + ';\n'
                self.parse(config + ''.join(equation))
                self.lexer.initialize(sentence, position)
                self.lexer.lex()
            else:
                sentence, position = self.lexer.sentence, self.lexer.mark()
                self.parse(self._generate_covdrv(tensor, [index[0] for index in indices]))
                self.lexer.initialize(sentence, position)
                self.lexer.lex()
        return self._instantiate_derivative(tensor, order, list(function.args[1:]) +
            [index[0] for index in indices], [index[1] for index in indices])

    # <CONFIG> -> '%' { <VARIABLE> }+ '[' <INTEGER> ']' ':' <SYMMETRY> { ',' { <VARIABLE> }+ '[' <INTEGER> ']' ':' <SYMMETRY> }*
    def _config(self):
        self.expect('COMMENT')
        while True:
            variable = []
            if self.peek('EULER'):
                self.lexer.token = 'SYMBOL'
            while any(self.peek(token) for token in ('SYMBOL', 'DIACRITIC', 'MATHOP')):
                variable.append(self._variable())
                if self.peek('EULER'):
                    self.lexer.token = 'SYMBOL'
            variable = ''.join(variable)
            self.expect('LEFT_BRACKET')
            dimension = self.lexer.lexeme
            self.expect('INTEGER')
            dimension = int(dimension)
            if self.dimension and self.dimension != dimension:
                raise TensorError('inconsistent tensor dimension')
            self.dimension = dimension
            self.expect('RIGHT_BRACKET')
            self.expect('COLON')
            symmetry = self.lexer.lexeme
            self.expect('SYMMETRY')
            self._instantiate_tensor(Tensor(variable, dimension), symmetry, invertible=(symmetry == 'metric'),
                permutation=(symmetry == 'permutation'), kronecker=(symmetry == 'kronecker'))
            if not self.accept('COMMA'): break

    # <TENSOR> -> <VARIABLE> [ ( '_' <LOWER_INDEX> ) | ( '^' <UPPER_INDEX> [ '_' <LOWER_INDEX> ] ) ]
    def _tensor(self):
        variable = self._strip(self._variable())
        symbol, indices = list(variable), []
        if self.accept('UNDERSCORE'):
            index, order = self._lower_index()
            indices.extend(index)
            symbol.extend((len(index) - order) * ['D'])
            function = Function('Tensor')(Symbol(''.join(symbol)), *indices)
            if order > 0:
                tensor = Tensor(function, self.dimension)
                return self._instantiate_derivative(tensor, order, indices)
            return function
        self.lexer.mark()
        if self.accept('CARET'):
            if self.accept('LEFT_BRACE'):
                if self.accept('LEFT_BRACE'):
                    self.lexer.reset()
                    symbol = ''.join(symbol)
                    sympy_env[symbol] = Symbol(symbol)
                    return sympy_env[symbol]
                self.lexer.reset()
                self.lexer.lex()
            index = self._upper_index()
            indices.extend(index)
            symbol.extend(len(index) * ['U'])
            if self.accept('UNDERSCORE'):
                index, order = self._lower_index()
                indices.extend(index)
                symbol.extend((len(index) - order) * ['D'])
                function = Function('Tensor')(Symbol(''.join(symbol)), *indices)
                if order > 0:
                    tensor = Tensor(function, self.dimension)
                    return self._instantiate_derivative(tensor, order, indices)
            return Function('Tensor')(Symbol(''.join(symbol)), *indices)
        symbol = ''.join(symbol)
        sympy_env[symbol] = Symbol(symbol)
        return sympy_env[symbol]

    # <LOWER_INDEX> -> <SYMBOL> | '{' { <SYMBOL> }* [ ',' { <SYMBOL> }+ ] '}'
    def _lower_index(self):
        indices, order = [], 0
        def append_index():
            index = self._strip(self.lexer.lexeme)
            self.expect('SYMBOL')
            indices.append(Symbol(index))
        if self.peek('SYMBOL'):
            append_index()
            return indices, order
        if self.accept('LEFT_BRACE'):
            while self.peek('SYMBOL'):
                append_index()
            if self.accept('COMMA'):
                while self.peek('SYMBOL'):
                    order += 1
                    append_index()
            self.expect('RIGHT_BRACE')
            return indices, order
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    # <UPPER_INDEX> -> <SYMBOL> | '{' { <SYMBOL> }+ '}'
    def _upper_index(self):
        indices = []
        def append_index():
            index = self._strip(self.lexer.lexeme)
            self.lexer.lex()
            indices.append(Symbol(index))
        if self.peek('SYMBOL') or self.peek('INTEGER'):
            append_index()
            return indices
        if self.accept('LEFT_BRACE'):
            while self.peek('SYMBOL') or self.peek('INTEGER'):
                append_index()
            self.expect('RIGHT_BRACE')
            return indices
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    # <VARIABLE> -> <SYMBOL> | <DIACRITIC> '{' <SYMBOL> '}' | <MATHOP> '{' <SYMBOL> { <SYMBOL> | <INTEGER> | <UNDERSCORE> }* '}'
    def _variable(self):
        variable = self.lexer.lexeme
        if self.accept('SYMBOL'):
            return variable
        if self.accept('DIACRITIC'):
            self.expect('LEFT_BRACE')
            symbol = self.lexer.lexeme
            self.expect('SYMBOL')
            self.expect('RIGHT_BRACE')
            return symbol + variable[1:]
        if self.peek('MATHOP'):
            self.expect('MATHOP')
            self.expect('LEFT_BRACE')
            symbol = [self.lexer.lexeme]
            if not (self.accept('SYMBOL') or self.accept('EULER')):
                sentence, position = self.lexer.sentence, self.lexer.mark()
                raise ParseError('unexpected \'%s\' at position %d' %
                    (sentence[position], position), sentence, position)
            while any(self.peek(token) for token in
                    ('SYMBOL', 'EULER', 'INTEGER', 'UNDERSCORE')):
                symbol.extend([self.lexer.lexeme])
                self.lexer.lex()
            self.expect('RIGHT_BRACE')
            return ''.join(symbol)
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    def _instantiate_tensor(self, tensor, symmetry, invertible=False, permutation=False, kronecker=False):
        def sgn(sequence):
            """ Permutation Signature (Parity)"""
            cycle_length = 0
            for n, i in enumerate(sequence[:-1]):
                for j in sequence[(n + 1):]:
                    if i == j: return 0
                    cycle_length += i > j
            return (-1)**cycle_length
        symbol, rank, dimension = tensor.symbol, tensor.rank, tensor.dimension
        symmetry = 'sym01' if symmetry == 'metric' \
              else None    if symmetry == 'nosym' \
              else symmetry
        if symbol in self.namespace:
            # avoid repeated instantiation
            array = self.namespace[symbol]
        elif symbol in Parser.namespace:
            # pylint: disable=unused-argument
            def formatwarning(message, category, filename=None, lineno=None, file=None, line=None):
                return '%s: %s\n' % (category.__name__, message)
            warnings.formatwarning = formatwarning
            # throw warning whenever duplicate namespace variable
            warnings.warn(symbol, OverrideWarning)
        if permutation:
            # instantiate permutation or Levi-Civita symbol using parity
            index  = [chr(105 + n) for n in range(rank)]
            prefix = '[' * rank + 'sgn([' + ', '.join(index) + '])'
            suffix = ''.join(' for %s in range(%d)]' % (index[rank - i], dimension) for i in range(1, rank + 1))
            array  = eval(prefix + suffix, {'sgn': sgn})
        elif kronecker:
            if rank != 2:
                raise TensorError('invalid rank for kronecker delta')
            array = ixp.declare_indexedexp(rank=rank, dimension=dimension)
            for i in range(dimension): array[i][i] = 1
        else:
            array = ixp.declare_indexedexp(rank, symbol, symmetry, dimension)
        if invertible:
            inverse_symbol = symbol.replace('U', 'D') if 'U' in symbol else symbol.replace('D', 'U')
            if dimension == 2:
                inverse, determinant = ixp.symm_matrix_inverter2x2(array)
            elif dimension == 3:
                inverse, determinant = ixp.symm_matrix_inverter3x3(array)
            elif dimension == 4:
                inverse, determinant = ixp.symm_matrix_inverter4x4(array)
            self.namespace[inverse_symbol] = inverse
            if symbol[-2:] == 'DD':
                self.namespace[symbol[:-2] + 'det'] = determinant
            else:
                if dimension == 2:
                    _, determinant = ixp.symm_matrix_inverter2x2(inverse)
                elif dimension == 3:
                    _, determinant = ixp.symm_matrix_inverter3x3(inverse)
                elif dimension == 4:
                    _, determinant = ixp.symm_matrix_inverter4x4(inverse)
                self.namespace[symbol[:-2] + 'det'] = determinant
        self.namespace[symbol] = array

    def _instantiate_derivative(self, tensor, order, indices, covariant=None):
        prefix, suffix = '_cd' if covariant else '_d', ''.join(covariant) if covariant else order * 'D'
        symbol = tensor.symbol + ('' if prefix in tensor.symbol else prefix) + suffix
        function, rank = Function('Tensor')(Symbol(symbol), *indices), tensor.rank
        if not covariant:
            if order == 2:
                if symmetry and symmetry != 'nosym':
                    symmetry = tensor.symmetry + '_sym%d%d' % (rank, rank + order - 1)
                else: symmetry = 'sym%d%d' % (rank, rank + order - 1)
            self._instantiate_tensor(Tensor(function, self.dimension), 'nosym')
        return function

    def _instantiate_christoffel(self, function):
        sentence, position = self.lexer.sentence, self.lexer.mark()
        if self.dimension is None:
            raise ParseError('cannot instantiate from inference', sentence, position)
        self.parse(self._generate_christoffel(function, self.dimension))
        self.lexer.initialize(sentence, position)
        self.lexer.lex()

    @staticmethod
    def _generate_christoffel(function, dimension):
        symbol, indices = '\\' + str(function.args[0])[:-3], function.args[1:]
        diacritic = 'bar'   if 'bar'   in symbol \
               else 'hat'   if 'hat'   in symbol \
               else 'tilde' if 'tilde' in symbol \
               else None
        metric = '\\%s{g}' % diacritic if diacritic else 'g'
        if diacritic: symbol = '\\%s{%s}' % (diacritic, symbol[:-len(diacritic)])
        indices = [('\\' if len(str(index)) > 1 else '') + str(index) for index in indices]
        bound_index = next(x for x in (chr(97 + n) for n in range(26)) if x not in indices)
        return (('% {metric}UU [{dim}]: metric;\n{symbol}^{i1}_{{{i2}{i3}}} = \\frac{{1}}{{2}} {metric}^{{{i1} {bound_index}}}({metric}_{{{i3} {bound_index},{i2}}} + {metric}_{{{bound_index} {i2},{i3}}} - {metric}_{{{i2} {i3},{bound_index}}})')
                .format(i1 = indices[0], i2 = indices[1], i3 = indices[2], symbol = symbol, metric = metric, bound_index = bound_index, dim = dimension))

    @staticmethod
    def _generate_covdrv(tensor, indices):
        order, LHS, RHS = len(indices), str(tensor), ''
        while order > 0:
            indices = [str(index[0]) for index in tensor.indexing] + \
                [str(index) for index in indices]
            for i, index in enumerate(indices):
                if index in indices[:i]:
                    alphabet = (chr(97 + n) for n in range(26))
                    indices[i] = next(x for x in alphabet if x not in indices)
            diff_index = indices[len(indices) - order]
            if len(str(diff_index)) > 1:
                diff_index = '\\' + str(diff_index)
            LHS = '\\nabla_%s %s' % (diff_index, LHS)
            RHS += '\\partial_%s %s' % (diff_index, str(tensor))
            for index, position in tensor.indexing:
                alphabet = (chr(97 + n) for n in range(26))
                bound_index  = next(x for x in alphabet if x not in indices)
                tensor_latex = Tensor.latex_format(tensor.symbol,
                    [(bound_index, index_[1]) if index_[0] == index else index_ for index_ in tensor.indexing])
                if len(str(index)) > 1:
                    index = '\\' + str(index)
                if position == 'U':
                    RHS += ' + \\Gamma^%s_{%s %s} %s' % (index, bound_index, diff_index, tensor_latex)
                else:
                    RHS += ' - \\Gamma^%s_{%s %s} %s' % (bound_index, index, diff_index, tensor_latex)
            order -= 1
        return LHS + ' = ' + RHS

    @staticmethod
    def _strip(symbol):
        return symbol[1:] if symbol[0] == '\\' else symbol

    def peek(self, token):
        return self.lexer.token == token

    def accept(self, token):
        if self.peek(token):
            self.lexer.lex()
            return True
        return False

    def expect(self, token):
        if not self.accept(token):
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('expected token %s at position %d' %
                (token, position), sentence, position)

class ParseError(Exception):
    """ Invalid LaTeX Sentence """

    def __init__(self, message, sentence, position):
        super().__init__('%s\n%s^\n' % (sentence, (12 + position) * ' ') + message)

class Tensor:
    """ Tensor Structure """

    def __init__(self, function, dimension):
        if isinstance(function, Function('Tensor')):
            self.function  = function
            self.dimension = dimension
            self.symbol    = str(function.args[0])
            self.rank      = len(function.args) - 1
            self.indexing  = list(zip(function.args[1:], re.findall(r'[UD]', self.symbol)))
        else:
            self.function  = None
            self.dimension = dimension
            self.symbol    = function
            self.rank      = len(re.findall(r'[UD]', self.symbol))
            self.indexing  = None

    @staticmethod
    def array_format(function):
        """ Tensor Notation for Array Formatting """
        fields = [str(arg) for arg in function.args]
        suffix = ''.join(['[' + n + ']' for n in fields[1:]])
        return fields[0] + suffix

    @staticmethod
    def latex_format(name, indexing):
        """ Tensor Notation for LaTeX Formatting """
        latex = [re.split('[UD]', name)[0], [], []]
        U_count, D_count = 0, 0
        for index, position in indexing:
            index = str(index)
            if len(index) > 1:
                index = '\\' + index
            if position == 'U':
                latex[1].append(index)
                U_count += 1
            else:
                latex[2].append(index)
                D_count += 1
        latex[1] = ' '.join(latex[1])
        latex[2] = ' '.join(latex[2])
        if U_count > 0:
            if U_count > 1:
                latex[1] = '^{' + latex[1] + '}'
            else: latex[1] = '^' + latex[1]
        if D_count > 0:
            if D_count > 1:
                latex[2] = '_{' + latex[2] + '}'
            else: latex[2] = '_' + latex[2]
        return ''.join(latex)

    def __str__(self):
        if self.indexing is None: return self.symbol
        return self.latex_format(self.symbol, self.indexing)

def _summation(equation, dimension):
    var, expr = equation
    # count every index on LHS to determine the rank
    rank = len(re.findall(r'\[[^\]]+\]', var))
    # construct a tuple list of every LHS free index
    LHS, RHS = zip(re.findall(r'\[([a-zA-Z]+)\]', var), re.findall(r'[UD]', var)), []
    # iterate over every subexpression containing a product
    for product in re.split(r'\s[\+\-]\s', expr):
        # extract every index present in the subexpression
        idx_list = re.findall(r'\[([a-zA-Z]+)\]', product)
        # extract every index position (ex: U or D)
        pos_list = re.findall(r'[UD]', product)
        free_index, bound_index = [], []
        # iterate over every unique index in the subexpression
        for idx in sorted(uniquify((idx_list))):
            count = U = D = 0; index_tuple = []
            # count index occurrence and position occurrence
            for idx_, pos_ in zip(idx_list, pos_list):
                if idx_ == idx:
                    index_tuple.append((idx_, pos_))
                    if pos_ == 'U': U += 1
                    if pos_ == 'D': D += 1
                    count += 1
            # identify every bound index on the RHS
            if count > 1:
                if count != 2 or U != D:
                    # raise exception upon violation of the following rule:
                    # a bound index must appear exactly once as a superscript
                    # and exactly once as a subscript in any single term
                    raise TensorError('illegal bound index')
                bound_index.append(idx)
            # identify every free index on the RHS
            else: free_index.extend(index_tuple)
        RHS.append(sorted(uniquify(free_index)))
        summation = product
        # generate implied summation over every bound index
        for idx in bound_index:
            summation = 'sum([%s for %s in range(%d)])' % \
                (summation, idx, dimension)
        expr = expr.replace(product, summation)
    if LHS:
        LHS = sorted(uniquify(LHS))
        for i in range(len(RHS)):
            if LHS != RHS[i]:
                # raise exception upon violation of the following rule:
                # a free index must appear in every term with the same
                # position and cannot be summed over in any term
                raise TensorError('unbalanced free index')
        # generate tensor instantiation with implied summation
        for idx, _ in LHS:
            expr = '[%s for %s in range(%d)]' % (expr, idx, dimension)
    if rank == len(re.findall(r'\[[^0-9\]]+\]', var)):
        return {var.split('[')[0]: expr}
    return {re.sub(r'\[[^0-9\]]+\]', '[:]', var): expr}

class TensorError(Exception):
    """ Invalid Tensor Indexing or Dimension """
class OverrideWarning(UserWarning):
    """ Overridden Namespace Variable """

def parse_expr(sentence):
    """ Convert LaTeX Sentence to SymPy Expression (Expression Mode)

        :arg:    latex sentence (raw string)
        :return: expression
    """
    return Parser().parse(sentence, expression=True)

def parse(sentence, evaluate=True):
    """ Convert LaTeX Sentence to SymPy Expression

        :arg:    latex sentence (raw string)
        :arg:    evaluate mode [default: enabled]
        :return: namespace
    """
    namespace, unevaled = Parser().parse(sentence)
    if evaluate:
        # update namespace with SymPy environment
        namespace.update(sympy_env)
        for key in unevaled:
            if isinstance(unevaled[key], str):
                # evaluate each implied summation and update namespace
                exec('%s = %s' % (key, unevaled[key]), namespace)
        # remove SymPy environment from namespace
        for key in sympy_env: namespace.pop(key)
        # update (static) class namespace for global persistance
        Parser.namespace.update(namespace)
    else: namespace.update(unevaled)
    # inject namespace into the previous stack frame
    frame = currentframe().f_back
    frame.f_globals.update(namespace)
    return list(namespace.keys())
