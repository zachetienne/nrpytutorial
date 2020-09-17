""" Convert LaTeX Sentence to SymPy Expression """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

# pylint: disable=attribute-defined-outside-init
from sympy import Function, Derivative, Symbol, Integer, Rational, Float, Pow
from sympy import sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
from sympy import pi, exp, log, sqrt, expand, diff, var
from collections import OrderedDict
from expr_tree import ExprTree
from inspect import currentframe
from warnings import warn
import indexedexp as ixp
import re

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
        symmetry = r'nosym|(?:sym|anti)[0-9]+(?:_(?:sym|anti)[0-9]+)*|metric|permutation'
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
              ('COVDRV_OP',      r'\\nabla|D'),
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
        <NLOG>          -> <NLOG_CMD> [ '_' <INTEGER> | { <INTEGER> } ] ( <SYMBOL> | <INTEGER> | '(' <EXPRESSION> ')' )
        <TRIG>          -> <TRIG_CMD> [ '^' <INTEGER> | { <INTEGER> } ] ( <SYMBOL> | <INTEGER> | '(' <EXPRESSION> ')' )
        <PARDRV>        -> { <PARDRV_OP> '_' <SYMBOL> }+ ( <TENSOR> | '(' <EXPRESSION> ')' )
        <COVDRV>        -> { <COVDRV_OP> ( '^' | '_' ) <SYMBOL> }+ ( <TENSOR> | '(' <EXPRESSION> ')' )
        <CONFIG>        -> '%' <VARIABLE> '[' <INTEGER> ']' [ ':' <SYMMETRY> ] { ',' <VARIABLE> '[' <INTEGER> ']' [ ':' <SYMMETRY> ] }*
        <TENSOR>        -> <VARIABLE> [ ( '_' <LOWER_INDEX> ) | ( '^' <UPPER_INDEX> [ '_' <LOWER_INDEX> ] ) ]
        <VARIABLE>      -> <SYMBOL> | <STRING> | <DIACRITIC> '{' <SYMBOL> '}'
        <STRING>        -> <MATHOP> '{' <SYMBOL> { <SYMBOL> | <INTEGER> | <UNDERSCORE> }* '}'
        <LOWER_INDEX>   -> <SYMBOL> | '{' { <SYMBOL> }* [ ',' { <SYMBOL> }+ ] '}'
        <UPPER_INDEX>   -> <SYMBOL> | '{' { <SYMBOL> }+ '}'
    """

    namespace = OrderedDict()

    def __init__(self):
        self.lexer     = Lexer()
        self.namespace = OrderedDict()
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
        return self._root()

    # <ROOT> -> <STRUCTURE> { <LINE_BREAK> <STRUCTURE> }*
    def _root(self):
        self._structure()
        while self.accept('LINE_BREAK'):
            self._structure()
        return self.namespace

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
        if variable.func == Function('Tensor'):
            variable = Tensor.notation(variable)
        else: variable = str(variable)
        self.expect('EQUAL')
        expr = self._expression()
        # distribute over every parenthetical expression
        tree = ExprTree(expand(expr))
        tensorial = False
        for subtree in tree.preorder():
            subexpr = subtree.expr
            if subexpr.func == Function('Tensor'):
                tensorial = True
        if tensorial:
            # iterate through every tensor function
            expr = str(tree.root.expr)
            for subtree in tree.preorder():
                subexpr = subtree.expr
                if subexpr.func == Function('Tensor'):
                    # replace function notation with array notation
                    func = Tensor.notation(subexpr)
                    expr = expr.replace(str(subexpr), func)
            # perform implied summation on tensorial equation
            expanded = _summation((variable, expr), self.dimension)
            # update namespace and lexer with expanded equation
            self.namespace.update(expanded)
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
                ('SYMBOL', 'DIACRITIC')):
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
                    tensor = self.namespace[str(function.args[1])]
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
        tensor = self.namespace[str(function.args[0])]
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
        tensor = self.namespace[str(function.args[0])]
        if location == 'RHS':
            if equation[2]:
                sentence, position = self.lexer.sentence, self.lexer.mark()
                if config: config = '% ' + ';\n'.join(config) + ';\n'
                self.parse(config + ''.join(equation))
                self.lexer.initialize(sentence, position)
                self.lexer.lex()
            else:
                sentence, position = self.lexer.sentence, self.lexer.mark()
                self.parse(self._generate_gradient(tensor, [index[0] for index in indices]))
                self.lexer.initialize(sentence, position)
                self.lexer.lex()
        return self._instantiate_derivative(tensor, order, list(function.args[1:]) +
            [index[0] for index in indices], [index[1] for index in indices])

    # <CONFIG> -> '%' <TENSOR> '[' <INTEGER> ']' [ ':' <SYMMETRY> ] { ',' <TENSOR> '[' <INTEGER> ']' [ ':' <SYMMETRY> ] }*
    def _config(self):
        self.expect('COMMENT')
        while True:
            function = self._tensor()
            self.expect('LEFT_BRACKET')
            dimension = self.lexer.lexeme
            self.expect('INTEGER')
            dimension = int(dimension)
            if self.dimension and self.dimension != dimension:
                raise TensorError('inconsistent tensor dimension')
            self.dimension = dimension
            self.expect('RIGHT_BRACKET')
            if self.accept('COLON'):
                symmetry = self.lexer.lexeme
                self.expect('SYMMETRY')
            else: symmetry = None
            tensor = Tensor(function, dimension, symmetry,
                invertible=(symmetry == 'metric'), permutation=(symmetry == 'permutation'))
            self.namespace[tensor.symbol] = tensor
            if tensor.inverse and tensor.rank == 2:
                inverse = tensor.symbol.replace('U', 'D') if 'U' in tensor.symbol else tensor.symbol.replace('D', 'U')
                self.namespace[inverse] = self.namespace[tensor.symbol].inverse
            if not self.accept('COMMA'): break

    # <TENSOR> -> <VARIABLE> [ ( '_' <LOWER_INDEX> ) | ( '^' <UPPER_INDEX> [ '_' <LOWER_INDEX> ] ) ]
    def _tensor(self):
        variable = self._strip(self._variable())
        symbol, indices = list(variable), []
        if self.accept('UNDERSCORE'):
            index, order = self._lower_index()
            indices.extend(index)
            symbol.extend((len(index) - order) * ['D'])
            if order > 0:
                tensor = self.namespace[''.join(symbol)]
                return self._instantiate_derivative(tensor, order, indices)
            return Function('Tensor')(Symbol(''.join(symbol)), *indices)
        self.lexer.mark()
        if self.accept('CARET'):
            if self.accept('LEFT_BRACE'):
                if self.accept('LEFT_BRACE'):
                    self.lexer.reset()
                    return var(''.join(symbol))
                self.lexer.reset()
                self.lexer.lex()
            index = self._upper_index()
            indices.extend(index)
            symbol.extend(len(index) * ['U'])
            if self.accept('UNDERSCORE'):
                index, order = self._lower_index()
                indices.extend(index)
                symbol.extend((len(index) - order) * ['D'])
                if order > 0:
                    tensor = self.namespace[''.join(symbol)]
                    return self._instantiate_derivative(tensor, order, indices)
            return Function('Tensor')(Symbol(''.join(symbol)), *indices)
        return var(''.join(symbol))

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
            self.expect('SYMBOL')
            indices.append(Symbol(index))
        if self.peek('SYMBOL'):
            append_index()
            return indices
        if self.accept('LEFT_BRACE'):
            while self.peek('SYMBOL'):
                append_index()
            self.expect('RIGHT_BRACE')
            return indices
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    # <VARIABLE> -> <SYMBOL> | <STRING> | <DIACRITIC> '{' <SYMBOL> '}'
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
            return self._string()
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    # <STRING> -> <MATHOP> '{' <SYMBOL> { <SYMBOL> | <INTEGER> | <UNDERSCORE> }* '}'
    def _string(self):
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

    def _instantiate_derivative(self, tensor, order, indices, covariant=None):
        prefix, suffix = '_cd' if covariant else '_d', ''.join(covariant) if covariant else order * 'D'
        symbol = tensor.symbol + ('' if prefix in tensor.symbol else prefix) + suffix
        rank, dimension, symmetry = tensor.rank, tensor.dimension, tensor.symmetry
        function = Function('Tensor')(Symbol(symbol), *indices)
        if not covariant:
            if order == 2:
                if symmetry and symmetry != 'nosym':
                    symmetry = tensor.symmetry + '_sym%d%d' % (rank, rank + order - 1)
                else: symmetry = 'sym%d%d' % (rank, rank + order - 1)
            self.namespace[symbol] = Tensor(function, dimension, symmetry)
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
        metric = '\\tilde{\\gamma}' if diacritic == 'tilde' \
            else '\\%s{g}' % diacritic if diacritic \
            else 'g'
        if diacritic: symbol = '\\%s{%s}' % (diacritic, symbol[:-len(diacritic)])
        indices = [('\\' if len(str(index)) > 1 else '') + str(index) for index in indices]
        bound_index = next(x for x in (chr(97 + n) for n in range(26)) if x not in indices)
        return (('% {metric}^{{{i1}{bound_index}}} [{dim}]: metric, {metric}_{{{i3}{bound_index}}} [{dim}]: metric, {metric}_{{{bound_index} {i2}}} [{dim}]: metric, {metric}_{{{i2}{i3}}} [{dim}]: metric;\n'
                '{symbol}^{i1}_{{{i2}{i3}}} = \\frac{{1}}{{2}} {metric}^{{{i1} {bound_index}}}({metric}_{{{i3} {bound_index},{i2}}} + {metric}_{{{bound_index} {i2},{i3}}} - {metric}_{{{i2} {i3},{bound_index}}})')
                .format(i1 = indices[0], i2 = indices[1], i3 = indices[2], symbol = symbol, metric = metric, bound_index = bound_index, dim = dimension))

    @staticmethod
    def _generate_gradient(tensor, indices):
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

    def __init__(self, function, dimension, symmetry=None, permutation=False, invertible=False):
        self.function = function
        self.symbol = str(function.args[0])
        self.rank = len(function.args) - 1
        self.dimension = dimension
        self.symmetry  = 'sym01' if symmetry == 'metric' \
                    else 'nosym' if symmetry == 'permutation' \
                    else symmetry
        self.indexing = list(zip(function.args[1:], re.findall(r'[UD]', self.symbol)))
        # avoid repeated instantiation of the same derivative
        if '_d' in self.symbol and self.symbol in Parser.namespace:
            self.struct = Parser.namespace[self.symbol]
        elif permutation:
            # instantiate permutation or Levi-Civita symbol using parity
            index = [chr(105 + n) for n in range(self.rank)]
            prefix = '[' * self.rank + 'sgn([' + ', '.join(index) + '])'
            suffix = ''.join(' for %s in range(%d)]' % (index[self.rank - i], dimension) for i in range(1, self.rank + 1))
            self.struct = eval(prefix + suffix, {'sgn': self._sgn})
        else: self.struct = ixp.declare_indexedexp(self.rank, self.symbol, self.symmetry, self.dimension)
        if invertible and self.rank > 0 and self.symmetry == 'sym01':
            if   self.dimension == 2:
                self.inverse = ixp.symm_matrix_inverter2x2(self.struct)[0]
            elif self.dimension == 3:
                self.inverse = ixp.symm_matrix_inverter3x3(self.struct)[0]
            elif self.dimension == 4:
                self.inverse = ixp.symm_matrix_inverter4x4(self.struct)[0]
        else: self.inverse = None

    @staticmethod
    def _sgn(sequence):
        """ Permutation Signature (Parity)"""
        cycle_length = 0
        for n, i in enumerate(sequence[:-1]):
            for j in sequence[(n + 1):]:
                if i == j: return 0
                cycle_length += i > j
        return (-1)**cycle_length

    @staticmethod
    def notation(function):
        """ Tensor Notation for Array Indexing """
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

    def __repr__(self):
        return '%s -> (dimension: %d, symmetry: %s)' % \
            (self.function, self.dimension, self.symmetry)

    def __str__(self):
        return self.latex_format(self.symbol, self.indexing)

def _summation(equation, dimension):
    var, expr = equation
    # construct a tuple list of every LHS free index
    LHS, RHS = set(zip(re.findall(r'\[([a-zA-Z]+)\]', var),
        re.findall(r'[UD]', var))), []
    # iterate over every subexpression containing a product
    for product in re.split(r'\s[\+\-]\s', expr):
        # extract every index present in the subexpression
        idx_list = re.findall(r'\[([a-zA-Z]+)\]', product)
        # extract every index position (ex: U or D)
        pos_list = re.findall(r'[UD]', product)
        free_index, bound_index = [], []
        # iterate over every unique index in the subexpression
        for idx in set(idx_list):
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
        RHS.append(set(free_index))
        summation = product
        # generate implied summation over every bound index
        for idx in bound_index:
            summation = 'sum([%s for %s in range(%d)])' % \
                (summation, idx, dimension)
        expr = expr.replace(product, summation)
    if LHS:
        for i in range(len(RHS)):
            if LHS != RHS[i]:
                # raise exception upon violation of the following rule:
                # a free index must appear in every term with the same
                # position and cannot be summed over in any term
                raise TensorError('unbalanced free index')
        # generate tensor instantiation with implied summation
        for idx, _ in LHS:
            expr = '[%s for %s in range(%d)]' % (expr, idx, dimension)
    return {var.split('[')[0]: expr}

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

def parse(sentence, debug=False):
    """ Convert LaTeX Sentence to SymPy Expression

        :arg:    latex sentence (raw string)
        :arg:    debug mode [default: disabled]
        :return: namespace
    """
    namespace = Parser().parse(sentence)
    if not debug:
        for var in namespace:
            # throw warning whenever duplicate namespace variable
            if var in Parser.namespace:
                warn('\'' + var + '\'', OverrideWarning, stacklevel=2)
            # extract array field from Tensor wrapper class
            if isinstance(namespace[var], Tensor):
                namespace[var] = namespace[var].struct
        # update global scope with tensor namespace
        globals().update(namespace)
        for var in namespace:
            if isinstance(namespace[var], str):
                # evaluate each implied summation
                namespace[var] = eval(namespace[var])
                globals()[var] = namespace[var]
        # update static namespace for persistance
        Parser.namespace.update(namespace)
    # inject namespace into the previous stack frame
    frame = currentframe().f_back
    frame.f_globals.update(namespace)
    return list(namespace.keys())

if __name__ == "__main__":
    import doctest
    doctest.testmod()
