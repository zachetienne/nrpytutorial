""" Convert LaTeX Sentence to SymPy Expression """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

from sympy import Function, Derivative, Symbol, Integer, Rational, Float, Pow, Add
from sympy import sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
from sympy import pi, exp, log, sqrt, expand, diff
from functional import chain, uniquify
from inspect import currentframe
from expr_tree import ExprTree
import indexedexp as ixp
import re, sys, warnings

# pylint: disable = attribute-defined-outside-init, protected-access, exec-used
sympy_env = (('sin', sin), ('cos', cos), ('tan', tan), ('sinh', sinh), ('cosh', cosh), ('tanh', tanh),
    ('asin', asin), ('acos', acos), ('atan', atan), ('asinh', asinh), ('acosh', acosh), ('atanh', atanh),
    ('pi', pi), ('exp', exp), ('log', log), ('sqrt', sqrt), ('diff', diff))

class Lexer:
    """ LaTeX Lexer

        The following class will tokenize a LaTeX sentence for parsing.
    """

    def __init__(self):
        greek_pattern = '|'.join(letter for letter in (r'\\[aA]lpha', r'\\[bB]eta', r'\\[gG]amma', r'\\[dD]elta',
            r'\\[eE]psilon', r'\\[zZ]eta', r'\\[eE]ta', r'\\[tT]heta', r'\\[iI]ota', r'\\[kK]appa', r'\\[lL]ambda',
            r'\\[mM]u', r'\\[nN]u', r'\\[xX]i', r'\\[oO]mikron', r'\\[pP]i', r'\\[Rr]ho', r'\\[sS]igma', r'\\[tT]au',
            r'\\[uU]psilon', r'\\[pP]hi', r'\\[cC]hi', r'\\[pP]si', r'\\[oO]mega'))
        symmetry = r'nosym|(?:sym|anti)[0-9]+(?:_(?:sym|anti)[0-9]+)*'
        # define a regex pattern for every token, create a named capture group for
        # every pattern, join together the resulting pattern list using a pipe symbol
        # for regex alternation, and compile the generated regular expression
        self.regex = re.compile('|'.join(['(?P<%s>%s)' % pattern for pattern in
            [ ('SPACE_DELIM',    r'(?:\s|\\,|\{\})+|\&'),
              ('DIMENSION',      r'[2-4]D'),
              ('RATIONAL',       r'[0-9]+\/[1-9]+|\\frac{[0-9]+}{[1-9]+}'),
              ('DECIMAL',        r'[0-9]+\.[0-9]+'),
              ('INTEGER',        r'\-?[0-9]+'),
              ('NABLA',          r'\\nabla'),
              ('PI',             r'\\pi'),
              ('EULER',          r'e'),
              ('PLUS',           r'\+'),
              ('MINUS',          r'\-'),
              ('DIVIDE',         r'\/'),
              ('EQUAL',          r'\='),
              ('CARET',          r'\^'),
              ('COMMA',          r'\,'),
              ('COLON',          r'\:'),
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
              ('PARTIAL',        r'\\partial'),
              ('SQRT_CMD',       r'\\sqrt'),
              ('FRAC_CMD',       r'\\frac'),
              ('TRIG_CMD',       r'\\sinh|\\cosh|\\tanh|\\sin|\\cos|\\tan'),
              ('NLOG_CMD',       r'\\ln|\\log'),
              ('VPHANTOM',       r'\\vphantom'),
              ('DEFINE_MACRO',   r'define'),
              ('ASSIGN_MACRO',   r'assign'),
              ('PARSE_MACRO',    r'parse'),
              ('INDEX_KWRD',     r'index'),
              ('BASIS_KWRD',     r'basis'),
              ('DERIV_KWRD',     r'deriv'),
              ('DERIV_TYPE',     r'symbolic|_d'),
              ('UNDERSCORE',     r'\_'),
              ('DIACRITIC',      r'\\hat|\\tilde|\\bar'),
              ('SYMMETRY',       r'const|metric|permutation|kronecker|' + symmetry),
              ('MATHOP',         r'\\mathop'),
              ('LETTER',         r'[a-zA-Z]|' + greek_pattern),
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
        <STRUCTURE>     -> <CONFIG> | <ASSIGNMENT> | <ENVIRONMENT>
        <CONFIG>        -> '%' ( <DEFINE> | <ASSIGN> | <PARSE> )
        <PARSE>         -> <PARSE_MACRO> <ASSIGNMENT>
        <ENVIRONMENT>   -> <BEGIN_ALIGN> ( <CONFIG> | <ASSIGNMENT> ) { <LINE_BREAK> ( <CONFIG> | <ASSIGNMENT> ) }* <END_ALIGN>
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
        <NLOG>          -> <NLOG_CMD> [ '_' ( <INTEGER> | { <INTEGER> } ) ] ( <TENSOR> | <INTEGER> | '(' <EXPRESSION> ')' )
        <TRIG>          -> <TRIG_CMD> [ '^' ( <INTEGER> | { <INTEGER> } ) ] ( <TENSOR> | <INTEGER> | '(' <EXPRESSION> ')' )
        <PARDRV>        -> [ <VPHANTOM> '{' <DERIV_TYPE> '}' ] { <PARTIAL> [ '^' <INTEGER> ] '_' <LETTER> }+ ( <TENSOR> | '(' <EXPRESSION> ')' )
        <COVDRV>        -> [ <VPHANTOM> '{' <DERIV_TYPE> '}' ] { ( <NABLA> | <DIACRITIC> '{' <NABLA> '}' ) ( '^' | '_' ) <LETTER> }+ ( <TENSOR> | '(' <EXPRESSION> ')' )
        <DEFINE>        -> <DEFINE_MACRO> ( <GLOBAL> | [ <SYMMETRY> ] { <SYMBOL> }+ [ '(' <DIMENSION> ')' ] ) { ',' ( <GLOBAL> | [ <SYMMETRY> ] { <SYMBOL> }+ [ '(' <DIMENSION> ')' ] ) }*
        <GLOBAL>        -> <BASIS_KWRD> <BASIS> | <DERIV_KWRD> <DERIV_TYPE> | <INDEX_KWRD> <RANGE>
        <BASIS>         -> <BASIS_KWRD> <LEFT_BRACKET> <LETTER> [ ',' <LETTER> ]* <RIGHT_BRACKET>
        <RANGE>         -> ( <LETTER> | '[' <LETTER> '-' <LETTER> ']' ) '=' <INTEGER> ':' <INTEGER>
        <SYMBOL>        -> <LETTER> | <DIACRITIC> '{' <LETTER> '}' | <MATHOP> '{' <LETTER> { <LETTER> | <INTEGER> | <UNDERSCORE> }* '}'
        <TENSOR>        -> <SYMBOL> [ ( '_' <LOWER_INDEX> ) | ( '^' <UPPER_INDEX> [ '_' <LOWER_INDEX> ] ) ]
        <LOWER_INDEX>   -> <LETTER> | <INTEGER> | '{' { <LETTER> | <INTEGER> }* [ ',' { <LETTER> }+ ] '}'
        <UPPER_INDEX>   -> <LETTER> | <INTEGER> | '{' { <LETTER> | <INTEGER> }+ '}'
    """

    _namespace = {}
    continue_parsing = True

    def __init__(self, debug=False):
        self.lexer = Lexer()
        if 'basis' not in self._namespace:
            self._namespace['basis'] = []
        if 'deriv' not in self._namespace:
            self._namespace['deriv'] = 'symbolic'
        if 'index' not in self._namespace:
            self._namespace['index'] = {chr(105 + n): (0, 3) for n in range(4)}
        def excepthook(exception_type, exception, traceback):
            if not debug:
                print('%s: %s' % (exception_type.__name__, exception))
            else: sys.__excepthook__(exception_type, exception, traceback)
        sys.excepthook = excepthook

    def parse(self, sentence, expression=False):
        """ Parse LaTeX Sentence

            :arg:    latex sentence (raw string)
            :arg:    expression mode [default: disabled]
            :return: namespace or expression
        """
        stack = []; i_1 = i_2 = i_3 = 0
        for i, lexeme in enumerate(sentence):
            if   lexeme == '(': stack.append(i)
            elif lexeme == ')': i_1, i_2 = stack.pop(), i + 1
            elif lexeme == ',' and sentence[i - 1] == '{':
                i_3 = sentence.find('}', i) + 1
                subexpr, indexing = sentence[i_1:i_2], sentence[i_2:i_3][3:-1]
                operator = ''.join('\\partial_' + index for index in indexing)
                sentence = sentence.replace(sentence[i_1:i_3], operator + ' ' + subexpr)
            elif lexeme == ';' and sentence[i - 1] == '{':
                i_3 = sentence.find('}', i) + 1
                subexpr, indexing = sentence[i_1:i_2], sentence[i_2:i_3][3:-1]
                operator = ''.join('\\nabla_' + index for index in indexing)
                sentence = sentence.replace(sentence[i_1:i_3], operator + ' ' + subexpr)
        self.lexer.initialize(sentence)
        self.lexer.lex()
        if expression:
            tree = ExprTree(self._expression())
            for subtree in tree.preorder():
                subexpr, rank = subtree.expr, len(subtree.expr.args)
                if subexpr.func == Function('Tensor') and rank == 1:
                    subtree.expr = subexpr.args[0]
                    del subtree.children[:]
            return tree.reconstruct()
        self._root()
        return self._namespace

    # <ROOT> -> <STRUCTURE> { <LINE_BREAK> <STRUCTURE> }*
    def _root(self):
        self._structure()
        while self.accept('LINE_BREAK') or self.peek('COMMENT'):
            self._structure()

    # <STRUCTURE> -> <CONFIG> | <ASSIGNMENT> | <ENVIRONMENT>
    def _structure(self):
        if self.peek('COMMENT'):
            self._config()
        elif self.peek('BEGIN_ALIGN'):
            self._environment()
        else: self._assignment()

    # <CONFIG> -> '%' ( <DEFINE> | <ASSIGN> | <PARSE> )
    def _config(self):
        self.expect('COMMENT')
        if self.peek('PARSE_MACRO'):
            self._parse()
        elif self.peek('DEFINE_MACRO'):
            self._define()
        elif self.peek('ASSIGN_MACRO'):
            self._assign()
        else:
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('unsupported macro at position %d' %
                position, sentence, position)

    # <PARSE> -> <PARSE_MACRO> <ASSIGNMENT>
    def _parse(self):
        self.expect('PARSE_MACRO')
        self._assignment()

    # <ENVIRONMENT> -> <BEGIN_ALIGN> ( <CONFIG> | <ASSIGNMENT> ) { <LINE_BREAK> ( <CONFIG> | <ASSIGNMENT> ) }* <END_ALIGN>
    def _environment(self):
        self.expect('BEGIN_ALIGN')
        if self.peek('COMMENT'):
            self._config()
        else: self._assignment()
        while self.accept('LINE_BREAK'):
            if self.peek('COMMENT'):
                self._config()
            else: self._assignment()
        self.expect('END_ALIGN')

    # <ASSIGNMENT> -> ( <TENSOR> | <COVDRV> ) = <EXPRESSION>
    def _assignment(self):
        covdrv = self.peek('NABLA')
        self.lexer.mark()
        if self.accept('DIACRITIC'):
            self.expect('LEFT_BRACE')
            covdrv = self.peek('NABLA')
            self.lexer.reset()
        LHS = self._covdrv('LHS') if covdrv else self._tensor()
        indexed = LHS.func == Function('Tensor') and len(LHS.args) > 1
        self.expect('EQUAL')
        tree = ExprTree(self._expression())
        if not indexed:
            for subtree in tree.preorder():
                subexpr, rank = subtree.expr, len(subtree.expr.args)
                if subexpr.func == Function('Tensor') and rank > 1:
                    indexed = True
        if indexed:
            function, RHS = LHS, expand(tree.root.expr)
            # perform implied summation on indexed expression
            (LHS, RHS), dimension = self._summation(LHS, RHS)
            global_env = dict(sympy_env)
            global_env.update(self._namespace)
            for key in global_env:
                if isinstance(global_env[key], Tensor):
                    global_env[key] = global_env[key].structure
                if isinstance(global_env[key], Function('Constant')):
                    global_env[key] = global_env[key].args[0]
            # evaluate each implied summation and update namespace
            exec('%s = %s' % (LHS, RHS), global_env)
            symbol = LHS.split('[')[0]
            tensor = Tensor(function, dimension, global_env[symbol])
            self._namespace.update({symbol: tensor})
        else:
            # replace each tensor function with symbol whenever scalar
            for subtree in tree.preorder():
                subexpr, rank = subtree.expr, len(subtree.expr.args) - 1
                if subexpr.func in (Function('Tensor'), Function('Constant')):
                    if rank == 0:
                        subtree.expr = subexpr.args[0]
                        del subtree.children[:]
            self._namespace.update({str(LHS.args[0]): tree.reconstruct()})

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
        while any(self.peek(token) for token in ('LEFT_PAREN', 'PARTIAL',
                'LETTER', 'RATIONAL', 'DECIMAL', 'INTEGER', 'NABLA', 'PI', 'EULER', 'DIACRITIC',
                'DIVIDE', 'COMMAND', 'SQRT_CMD', 'FRAC_CMD', 'TRIG_CMD', 'NLOG_CMD', 'VPHANTOM')):
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
        self.lexer.mark()
        if self.accept('DIACRITIC'):
            self.expect('LEFT_BRACE')
            if self.peek('NABLA'):
                self.lexer.reset()
                return self._operator()
            self.lexer.reset()
        if any(self.peek(token) for token in
                ('RATIONAL', 'DECIMAL', 'INTEGER', 'PI')):
            return self._number()
        if any(self.peek(token) for token in
                ('LETTER', 'DIACRITIC', 'MATHOP')):
            sentence, position = self.lexer.sentence, self.lexer.mark()
            tensor = Tensor(self._tensor(), None)
            symbol, indexing = tensor.symbol, tensor.indexing
            if symbol[:5] == 'Gamma': # reserved keyword for christoffel symbol
                metric = 'gamma' if symbol[5:-3] == 'tilde' else 'g'
                if (metric + symbol[5:-3] + 'DD') not in self._namespace:
                    raise ParseError('cannot generate christoffel symbol without defined metric \'%s\'' %
                        (metric + symbol[5:-3]), sentence, position)
                sentence, position = self.lexer.sentence, self.lexer.mark()
                self.parse(self._generate_christoffel(tensor.function))
                self.lexer.initialize(sentence, position)
                self.lexer.lex()
            if symbol not in self._namespace:
                if indexing:
                    raise ParseError('cannot index undefined tensor \'%s\' at position %d' %
                        (symbol, position), sentence, position)
                self._define_tensor(Tensor(tensor.function, 0))
            return tensor.function
        if any(self.peek(token) for token in
                ('COMMAND', 'SQRT_CMD', 'FRAC_CMD', 'NLOG_CMD', 'TRIG_CMD')):
            return self._command()
        if any(self.peek(token) for token in
                ('VPHANTOM', 'PARTIAL', 'NABLA')):
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
        self.lexer.mark()
        if self.accept('VPHANTOM'):
            self.lexer.lex()
            self.lexer.lex()
            self.lexer.lex()
            operator = self.lexer.lexeme
            if self.peek('PARTIAL'):
                self.lexer.reset()
                return self._pardrv()
            if self.peek('NABLA') or self.peek('DIACRITIC'):
                self.lexer.reset()
                return self._covdrv('RHS')
            self.lexer.reset()
        if self.peek('PARTIAL'):
            return self._pardrv()
        if self.peek('NABLA') or self.peek('DIACRITIC'):
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

    # <NLOG> -> <NLOG_CMD> [ '_' ( <INTEGER> | { <INTEGER> } ) ] ( <TENSOR> | <INTEGER> | '(' <EXPRESSION> ')' )
    def _nlog(self):
        func = self._strip(self.lexer.lexeme)
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
        if self.peek('LETTER'):
            sentence, position = self.lexer.sentence, self.lexer.mark()
            tensor = Tensor(self._tensor(), None)
            symbol, indexing = tensor.symbol, tensor.indexing
            if symbol not in self._namespace:
                if indexing:
                    raise ParseError('cannot index undefined tensor \'%s\' at position %d' %
                        (symbol, position), sentence, position)
                self._define_tensor(Tensor(tensor.function, 0))
            expr = tensor.function
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

    # <TRIG> -> <TRIG_CMD> [ '^' ( <INTEGER> | { <INTEGER> } ) ] ( <TENSOR> | <INTEGER> | '(' <EXPRESSION> ')' )
    def _trig(self):
        func = self._strip(self.lexer.lexeme)
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
        if self.peek('LETTER'):
            sentence, position = self.lexer.sentence, self.lexer.mark()
            tensor = Tensor(self._tensor(), None)
            symbol, indexing = tensor.symbol, tensor.indexing
            if symbol not in self._namespace:
                if indexing:
                    raise ParseError('cannot index undefined tensor \'%s\' at position %d' %
                        (symbol, position), sentence, position)
                self._define_tensor(Tensor(tensor.function, 0))
            expr = tensor.function
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

    # <PARDRV> -> [ <VPHANTOM> '{' <DERIV_TYPE> '}' ] { <PARTIAL> [ '^' <INTEGER> ] '_' <LETTER> }+ ( <TENSOR> | '(' <EXPRESSION> ')' )
    def _pardrv(self):
        deriv_type = self._namespace['deriv']
        if self.accept('VPHANTOM'):
            self.expect('LEFT_BRACE')
            deriv_type = self.lexer.lexeme
            self.expect('DERIV_TYPE')
            self.expect('RIGHT_BRACE')
        indexing, order = [], 1
        while self.accept('PARTIAL'):
            if self.accept('CARET'):
                order = self.lexer.lexeme
                self.expect('INTEGER')
            self.expect('UNDERSCORE')
            index = self._strip(self.lexer.lexeme)
            self.expect('LETTER')
            indexing.extend(int(order)*[Symbol(index)])
        if all(index in self._namespace['basis'] for index in indexing):
            deriv_type = 'symbolic'
        if self.accept('LEFT_PAREN'):
            if deriv_type == 'symbolic':
                derivative = Derivative(self._expression(), *indexing)
                self.expect('RIGHT_PAREN')
                return derivative
            if deriv_type == '_d':
                tree = ExprTree(self._expression())
                self.expect('RIGHT_PAREN')
                for subtree in tree.preorder():
                    subexpr = subtree.expr
                    if subexpr.func == Function('Tensor'):
                        # insert temporary symbol '_x' for symbolic differentiation
                        subtree.expr = Function('_Tensor')(subexpr, Symbol('_x'))
                        del subtree.children[:]
                expr = tree.reconstruct()
                # differentiate the expression, including product rule expansion
                tree = ExprTree(diff(expr, Symbol('_x')))
                for subtree in tree.preorder():
                    subexpr = subtree.expr
                    if subexpr.func == Derivative:
                        # remove temporary symbol '_x' from tensor function
                        symbol, order = str(subexpr.args[0].args[0].args[0]), len(indexing)
                        _indexing = list(subexpr.args[0].args[0].args[1:]) + indexing
                        if symbol not in self._namespace:
                            raise ParseError('cannot differentiate undefined tensor \'%s\'' %
                                symbol, self.lexer.sentence)
                        tensor = self._namespace[symbol]
                        # instantiate partial derivative and update namespace
                        symbol = symbol + ('' if '_d' in symbol else '_d') + order * 'D'
                        function, rank = Function('Tensor')(symbol, *_indexing), len(_indexing) - order
                        symmetry = 'sym%d%d' % (rank, rank + order - 1) if order == 2 else 'nosym'
                        self._define_tensor(Tensor(function, tensor.dimension), symmetry)
                        subtree.expr = function
                        del subtree.children[:]
                    elif subexpr.func == Function('_Tensor'):
                        # remove temporary symbol '_x' from tensor function
                        subtree.expr = subexpr.args[0]
                        del subtree.children[:]
                return tree.reconstruct()
        if deriv_type == 'symbolic':
            return Derivative(self._tensor(), *indexing)
        sentence, position = self.lexer.sentence, self.lexer.mark()
        function = self._tensor()
        symbol, order = str(function.args[0]), len(indexing)
        indexing = list(function.args[1:]) + indexing
        if symbol not in self._namespace:
            raise ParseError('cannot differentiate undefined tensor \'%s\' at position %d' %
                (symbol, position), sentence, position)
        tensor = self._namespace[symbol]
        # instantiate partial derivative and update namespace
        symbol = symbol + ('' if '_d' in symbol else '_d') + order * 'D'
        function, rank = Function('Tensor')(symbol, *indexing), len(indexing) - order
        symmetry = 'sym%d%d' % (rank, rank + order - 1) if order == 2 else 'nosym'
        self._define_tensor(Tensor(function, tensor.dimension), symmetry)
        return function

    # <COVDRV> -> [ <VPHANTOM> '{' <DERIV_TYPE> '}' ] { ( <NABLA> | <DIACRITIC> '{' <NABLA> '}' ) ( '^' | '_' ) <LETTER> }+ ( <TENSOR> | '(' <EXPRESSION> ')' )
    def _covdrv(self, location):
        deriv_type = self._namespace['deriv']
        if self.accept('VPHANTOM'):
            self.expect('LEFT_BRACE')
            deriv_type = self.lexer.lexeme
            self.expect('DERIV_TYPE')
            self.expect('RIGHT_BRACE')
        indexing, equation, diacritic = [], ['', ' = ', '', ''], ''
        sentence, position = self.lexer.sentence, self.lexer.mark()
        while self.peek('NABLA') or self.peek('DIACRITIC'):
            lexeme = self._strip(self.lexer.lexeme)
            operator = '\\nabla'
            if self.accept('DIACRITIC'):
                diacritic = lexeme
                operator = '\\%s{\\nabla}' % diacritic
                self.expect('LEFT_BRACE')
                self.expect('NABLA')
                self.expect('RIGHT_BRACE')
            else: self.expect('NABLA')
            metric = 'gamma' if diacritic == 'tilde' else 'g'
            if (metric + diacritic + 'DD') not in self._namespace:
                raise ParseError('cannot generate covariant derivative without defined metric \'%s\'' %
                    (metric + diacritic), sentence, position)
            equation[0] += operator
            if deriv_type != 'symbolic':
                operator = '\\vphantom{%s} \\nabla' % deriv_type
            equation[3] += operator
            if self.accept('CARET'):
                index = self.lexer.lexeme
                equation[0] += '^' + index + ' '
                bound_index = next(x for x in (chr(97 + n) for n in range(26)) if x != index)
                metric = '\\%s{g}' % diacritic if diacritic in ('bar', 'hat') \
                    else '\\%s{\\gamma}' % diacritic if diacritic == 'tilde' \
                    else 'g'
                equation[2] += '%s^{%s %s} ' % (metric, index, bound_index)
                equation[3] += '_' + bound_index + ' '
                index = self._strip(index)
                self.expect('LETTER')
                indexing.append((Symbol(index), 'U'))
            elif self.accept('UNDERSCORE'):
                index = self.lexer.lexeme
                equation[0] += '_' + index + ' '
                equation[3] += '_' + index + ' '
                index = self._strip(index)
                self.expect('LETTER')
                indexing.append((Symbol(index), 'D'))
            else:
                sentence, position = self.lexer.sentence, self.lexer.mark()
                raise ParseError('unexpected \'%s\' at position %d' %
                    (sentence[position], position), sentence, position)
        # instantiate covariant derivative and update namespace
        marker_1 = self.lexer.mark()
        function = self._tensor()
        marker_2 = self.lexer.index
        equation[0] += self.lexer.sentence[marker_1:marker_2].strip()
        equation[3] += self.lexer.sentence[marker_1:marker_2].strip()
        if location == 'RHS':
            global_deriv = self._namespace['deriv']
            if equation[2]:
                sentence, position = self.lexer.sentence, self.lexer.mark()
                self.parse(''.join(equation))
                self.lexer.initialize(sentence, position)
                self.lexer.lex()
            else:
                sentence, position = self.lexer.sentence, self.lexer.mark()
                symbol = str(function.args[0])
                if symbol not in self._namespace:
                    sentence, position = self.lexer.sentence, self.lexer.mark()
                    raise ParseError('cannot differentiate undefined tensor \'%s\'' %
                        (symbol), sentence, position)
                tensor = Tensor(function, self._namespace[symbol].dimension)
                self.parse(self._generate_covdrv(tensor, indexing, deriv_type, diacritic))
                self.lexer.initialize(sentence, position)
                self.lexer.lex()
            self._namespace['deriv'] = global_deriv
        symbol, suffix = str(function.args[0]), ''.join([index[1] for index in indexing])
        symbol = symbol + ('' if '_cd' in symbol else '_cd' + diacritic) + suffix
        indexing = list(function.args[1:]) + [index[0] for index in indexing]
        return Function('Tensor')(symbol, *indexing)

    # <DEFINE> -> <DEFINE_MACRO> ( <GLOBAL> | [ <SYMMETRY> ] { <SYMBOL> }+ [ '(' <INTEGER> ')' ] ) { ',' ( <GLOBAL> | [ <SYMMETRY> ] { <SYMBOL> }+ [ '(' <INTEGER> ')' ] ) }*
    def _define(self):
        self.expect('DEFINE_MACRO')
        while True:
            if any(self.peek(kwrd) for kwrd in ('BASIS_KWRD', 'DERIV_KWRD', 'INDEX_KWRD')):
                self._global()
            else:
                symmetry = self.lexer.lexeme
                if self.peek('SYMMETRY'):
                    self.lexer.lex()
                else: symmetry = None
                symbol = []
                if self.peek('EULER'):
                    self.lexer.token = 'LETTER'
                while any(self.peek(token) for token in ('LETTER', 'DIACRITIC', 'MATHOP')):
                    symbol.append(self._symbol())
                    if self.peek('EULER'):
                        self.lexer.token = 'LETTER'
                symbol = ''.join(symbol)
                if self.accept('LEFT_PAREN'):
                    dimension = self.lexer.lexeme[:-1]
                    self.expect('DIMENSION')
                    dimension = int(dimension)
                    self.expect('RIGHT_PAREN')
                else: dimension = None
                if symmetry != 'const' and not dimension:
                    raise TensorError('dimension only omittable for constant')
                if symmetry == 'const':
                    self._namespace[symbol] = Function('Constant')(Symbol(symbol))
                else:
                    tensor = Tensor(Function('Tensor')(symbol), dimension)
                    self._define_tensor(tensor, symmetry, invertible=(symmetry == 'metric'),
                        permutation=(symmetry == 'permutation'), kronecker=(symmetry == 'kronecker'))
            if not self.accept('COMMA'): break

    # <GLOBAL> -> <BASIS_KWRD> <BASIS> | <DERIV_KWRD> <DERIV_TYPE> | <INDEX_KWRD> <RANGE>
    def _global(self):
        if self.accept('BASIS_KWRD'):
            self._basis()
        elif self.accept('DERIV_KWRD'):
            deriv_type = self.lexer.lexeme
            self.expect('DERIV_TYPE')
            self._namespace['deriv'] = deriv_type
        elif self.accept('INDEX_KWRD'):
            self._range()
        else:
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('unexpected keyword at position %d' %
                position, sentence, position)

    # <BASIS> -> <BASIS_KWRD> <LEFT_BRACKET> <LETTER> [ ',' <LETTER> ]* <RIGHT_BRACKET>
    def _basis(self):
        self.expect('LEFT_BRACKET')
        self._namespace['basis'].clear()
        while True:
            symbol = self._strip(self.lexer.lexeme)
            self.expect('LETTER')
            if symbol in self._namespace['basis']:
                sentence, position = self.lexer.sentence, self.lexer.mark()
                raise ParseError('duplicate basis symbol \'%s\' at position %d' %
                    (sentence[position], position), sentence, position)
            self._namespace['basis'].append(Symbol(symbol))
            if not self.accept('COMMA'): break
        self.expect('RIGHT_BRACKET')

    # <RANGE> -> ( <LETTER> | '[' <LETTER> '-' <LETTER> ']' ) '=' <INTEGER> ':' <INTEGER>
    def _range(self):
        if self.accept('LEFT_BRACKET'):
            index_1 = self._strip(self.lexer.lexeme)
            self.expect('LETTER')
            self.expect('MINUS')
            index_2 = self._strip(self.lexer.lexeme)
            self.expect('LETTER')
            index = [chr(i) for i in range(ord(index_1), ord(index_2) + 1)]
            self.expect('RIGHT_BRACKET')
        else:
            index = [self._strip(self.lexer.lexeme)]
            self.expect('LETTER')
        self.expect('EQUAL')
        lower = self.lexer.lexeme
        self.expect('INTEGER')
        self.expect('COLON')
        upper = self.lexer.lexeme
        self.expect('INTEGER')
        self._namespace['index'].update({i: (int(lower), int(upper) + 1) for i in index})

    # <ASSIGN> -> <ASSIGN_MACRO> [ <SYMMETRY> ] { <SYMBOL> }+
    def _assign(self):
        self.expect('ASSIGN_MACRO')
        symmetry = self.lexer.lexeme
        if self.peek('SYMMETRY'):
            self.lexer.lex()
        symbol = []
        if self.peek('EULER'):
            self.lexer.token = 'LETTER'
        while any(self.peek(token) for token in ('LETTER', 'DIACRITIC', 'MATHOP')):
            symbol.append(self._symbol())
            if self.peek('EULER'):
                self.lexer.token = 'LETTER'
        symbol = ''.join(symbol)
        if symbol not in self._namespace:
            raise TensorError('cannot update undefined tensor \'%s\'' % symbol)
        tensor = self._namespace[symbol]
        structure, dimension = tensor.structure, tensor.dimension
        if symmetry == 'metric':
            inverse_symbol = symbol.replace('U', 'D') if 'U' in symbol else symbol.replace('D', 'U')
            if dimension == 2:
                inverse, determinant = ixp.symm_matrix_inverter2x2(structure)
            elif dimension == 3:
                inverse, determinant = ixp.symm_matrix_inverter3x3(structure)
            elif dimension == 4:
                inverse, determinant = ixp.symm_matrix_inverter4x4(structure)
            _symbol = symbol.replace('U', 'D') if 'U' in symbol else symbol.replace('D', 'U')
            function = Function('Tensor')(_symbol, *tensor.function.args[1:])
            self._namespace[_symbol] = Tensor(function, dimension, inverse)
            _symbol = symbol[:-2] + 'det'
            function = Function('Tensor')(_symbol)
            self._namespace[_symbol] = Tensor(function, 0, determinant \
                if symbol[-2:] == 'DD' else (determinant)**(-1))
            self._namespace.update(self._namespace)

    # <TENSOR> -> <SYMBOL> [ ( '_' <LOWER_INDEX> ) | ( '^' <UPPER_INDEX> [ '_' <LOWER_INDEX> ] ) ]
    def _tensor(self):
        indexing = []
        symbol = list(self._strip(self._symbol()))
        if self.accept('UNDERSCORE'):
            index, order, _ = self._lower_index()
            indexing.extend(index)
            symbol.extend((len(index) - order) * ['D'])
            if order > 0:
                if symbol not in self._namespace:
                    sentence, position = self.lexer.sentence, self.lexer.mark()
                    raise ParseError('cannot differentiate undefined tensor \'%s\'' %
                        (symbol), sentence, position)
                tensor = self._namespace[symbol]
                symbol = ''.join(symbol) + ('' if '_d' in ''.join(symbol) else '_d') + order * 'D'
                function, rank = Function('Tensor')(symbol, *indexing), len(indexing) - order
                symmetry = 'sym%d%d' % (rank, rank + order - 1) if order == 2 else 'nosym'
                self._define_tensor(Tensor(function, tensor.dimension), symmetry)
                return function
        self.lexer.mark()
        if self.accept('CARET'):
            if self.accept('LEFT_BRACE'):
                if self.accept('LEFT_BRACE'):
                    self.lexer.reset()
                    symbol = ''.join(symbol)
                    function = Function('Tensor')(symbol)
                    if symbol in self._namespace:
                        if isinstance(self._namespace[symbol], Function('Constant')):
                            return self._namespace[symbol]
                    return function
                self.lexer.reset()
                self.lexer.lex()
            index = self._upper_index()
            indexing.extend(index)
            symbol.extend(len(index) * ['U'])
            if self.accept('UNDERSCORE'):
                index, order, _ = self._lower_index()
                indexing.extend(index)
                symbol.extend((len(index) - order) * ['D'])
                if order > 0:
                    if symbol not in self._namespace:
                        sentence, position = self.lexer.sentence, self.lexer.mark()
                        raise ParseError('cannot differentiate undefined tensor \'%s\'' %
                            symbol, sentence, position)
                    tensor = self._namespace[symbol]
                    symbol = ''.join(symbol) + ('' if '_d' in ''.join(symbol) else '_d') + order * 'D'
                    function, rank = Function('Tensor')(symbol, *indexing), len(indexing) - order
                    symmetry = 'sym%d%d' % (rank, rank + order - 1) if order == 2 else 'nosym'
                    self._define_tensor(Tensor(function, tensor.dimension), symmetry)
                    return function
        scalar, symbol = len(symbol) == 1, ''.join(symbol)
        if scalar and symbol in self._namespace:
            if isinstance(self._namespace[symbol], Function('Constant')):
                return self._namespace[symbol]
        return Function('Tensor')(symbol, *indexing)

    # <LOWER_INDEX> -> <LETTER> | <INTEGER> | '{' { <LETTER> | <INTEGER> }* [ ',' { <LETTER> }+ ] '}'
    def _lower_index(self):
        indexing, covariant = [], False
        def append_index():
            index = self._strip(self.lexer.lexeme)
            self.lexer.lex()
            indexing.append(Symbol(index))
        order = 0
        if self.peek('LETTER') or self.peek('INTEGER'):
            append_index()
            return indexing, order, covariant
        if self.accept('LEFT_BRACE'):
            while self.peek('LETTER') or self.peek('INTEGER'):
                append_index()
            # TODO: ADD SEMICOLON NOTATION FOR COVARIANT DERIVATIVE
            if self.accept('COMMA'):
                while self.peek('LETTER'):
                    order += 1
                    append_index()
            elif self.accept('SEMICOLON'):
                covariant = True
                while self.peek('LETTER'):
                    order += 1
                    append_index()
            self.expect('RIGHT_BRACE')
            return indexing, order, covariant
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    # <UPPER_INDEX> -> <LETTER> | <INTEGER> | '{' { <LETTER> | <INTEGER> }+ '}'
    def _upper_index(self):
        indexing = []
        def append_index():
            index = self._strip(self.lexer.lexeme)
            self.lexer.lex()
            indexing.append(Symbol(index))
        if self.peek('LETTER') or self.peek('INTEGER'):
            append_index()
            return indexing
        if self.accept('LEFT_BRACE'):
            while self.peek('LETTER') or self.peek('INTEGER'):
                append_index()
            self.expect('RIGHT_BRACE')
            return indexing
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    # <SYMBOL> -> <LETTER> | <DIACRITIC> '{' <LETTER> '}' | <MATHOP> '{' <LETTER> { <LETTER> | <INTEGER> | <UNDERSCORE> }* '}'
    def _symbol(self):
        lexeme = self.lexer.lexeme
        if self.accept('LETTER'):
            return lexeme
        if self.accept('DIACRITIC'):
            self.expect('LEFT_BRACE')
            symbol = self.lexer.lexeme + lexeme[1:]
            self.expect('LETTER')
            self.expect('RIGHT_BRACE')
            return symbol
        if self.peek('MATHOP'):
            self.expect('MATHOP')
            self.expect('LEFT_BRACE')
            symbol = [self.lexer.lexeme]
            if not (self.accept('LETTER') or self.accept('EULER')):
                sentence, position = self.lexer.sentence, self.lexer.mark()
                raise ParseError('unexpected \'%s\' at position %d' %
                    (sentence[position], position), sentence, position)
            while any(self.peek(token) for token in
                    ('LETTER', 'EULER', 'INTEGER', 'UNDERSCORE')):
                symbol.extend([self.lexer.lexeme])
                self.lexer.lex()
            self.expect('RIGHT_BRACE')
            return ''.join(symbol).replace('\\', '')
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    def _define_tensor(self, tensor, symmetry=None, invertible=False, permutation=False, kronecker=False):
        def sgn(sequence):
            """ Permutation Signature (Parity)"""
            cycle_length = 0
            for n, i in enumerate(sequence[:-1]):
                for j in sequence[(n + 1):]:
                    if i == j: return 0
                    cycle_length += i > j
            return (-1)**cycle_length
        symbol, rank = tensor.symbol, tensor.rank
        if symbol not in self._namespace:
            symmetry = 'sym01' if symmetry == 'metric' \
                  else None    if symmetry == 'nosym' \
                  else symmetry
            dimension = tensor.dimension
            if permutation:
                # instantiate permutation (Levi-Civita) symbol using parity
                index  = [chr(105 + n) for n in range(rank)]
                prefix = '[' * rank + 'sgn([' + ', '.join(index) + '])'
                suffix = ''.join(' for %s in range(%d)]' % (index[rank - i], dimension) for i in range(1, rank + 1))
                tensor.structure  = eval(prefix + suffix, {'sgn': sgn})
            elif kronecker:
                if rank != 2:
                    raise TensorError('cannot instantiate kronecker delta of rank ' + str(rank))
                tensor.structure = ixp.declare_indexedexp(rank=rank, dimension=dimension)
                for i in range(dimension): tensor.structure[i][i] = 1
            elif rank == 0:
                tensor.structure = Symbol(symbol)
            else:
                tensor.structure = ixp.declare_indexedexp(rank, symbol, symmetry, dimension)
            if invertible:
                if rank != 2:
                    raise TensorError('cannot invert tensor of rank ' + str(rank))
                if dimension == 2:
                    inverse, determinant = ixp.symm_matrix_inverter2x2(tensor.structure)
                elif dimension == 3:
                    inverse, determinant = ixp.symm_matrix_inverter3x3(tensor.structure)
                elif dimension == 4:
                    inverse, determinant = ixp.symm_matrix_inverter4x4(tensor.structure)
                _symbol = symbol.replace('U', 'D') if 'U' in symbol else symbol.replace('D', 'U')
                function = Function('Tensor')(_symbol, *tensor.function.args[1:])
                self._namespace[_symbol] = Tensor(function, dimension, inverse)
                _symbol = symbol[:-2] + 'det'
                function = Function('Tensor')(_symbol)
                self._namespace[_symbol] = Tensor(function, 0, determinant \
                    if symbol[-2:] == 'DD' else (determinant)**(-1))
            if symbol in self._namespace:
                # pylint: disable=unused-argument
                def formatwarning(message, category, filename=None, lineno=None, file=None, line=None):
                    return '%s: %s\n' % (category.__name__, message)
                warnings.formatwarning = formatwarning
                # throw warning whenever duplicate namespace variable
                warnings.warn(symbol, OverrideWarning)
            self._namespace[symbol] = tensor

    def _summation(self, LHS, RHS):
        def replace_function(sentence, subexpr, idx_map):
            # replace every tensor function with array notation
            tree = ExprTree(subexpr)
            for subtree in tree.preorder():
                subexpr = subtree.expr
                if subexpr.func == Function('Tensor'):
                    symbol = str(subexpr.args[0])
                    dimension = self._namespace[symbol].dimension
                    tensor = Tensor(subexpr, dimension)
                    for index in subexpr.args[1:]:
                        if str(index) in self._namespace['index']:
                            lower, upper = self._namespace['index'][str(index)]
                        else: lower, upper = (0, dimension)
                        if str(index) in idx_map and (lower, upper) != idx_map[str(index)]:
                            raise ParseError('inconsistent indexing range for index \'%s\'' %
                                index, self.lexer.sentence)
                        idx_map[str(index)] = (lower, upper)
                    sentence = sentence.replace(str(subexpr), tensor.array_format())
                elif subexpr.func == Function('Constant'):
                    symbol = str(subexpr.args[0])
                    sentence = sentence.replace(str(subexpr), symbol)
            return sentence
        def separate_indexing(subexpr):
            # extract every index present in the subexpression
            idx_list = re.findall(r'\[([a-zA-Z]+)\]', subexpr)
            # extract every index position (ex: U or D)
            pos_list = re.findall(r'[UD]', subexpr)
            if len(idx_list) != len(pos_list):
                pos_list.extend((len(idx_list) - len(pos_list)) * ['D'])
            free_index, bound_index = [], []
            # iterate over every unique index in the subexpression
            for idx in uniquify((idx_list)):
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
            return uniquify(free_index), bound_index
        iterable = RHS.args if RHS.func == Add else [RHS]
        LHS, RHS = Tensor(LHS, None).array_format(), str(RHS)
        # count every index on LHS to determine the rank
        rank = len(re.findall(r'\[[^\]]+\]', LHS))
        # construct a tuple list of every LHS free index
        free_index_LHS, _ = separate_indexing(LHS)
        # construct a tuple list of every RHS free index
        free_index_RHS = []
        for element in iterable:
            original, idx_map = str(element), {}
            if original[0] == '-':
                original = original[1:]
            modified = original
            tree = ExprTree(element)
            for subtree in tree.preorder():
                subexpr = subtree.expr
                if subexpr.func == Derivative:
                    argument = subexpr.args[0]
                    derivative = 'diff(' + str(argument)
                    argument = replace_function(str(argument), argument, idx_map)
                    free_index, _ = separate_indexing(argument)
                    for idx, _ in reversed(free_index):
                        lower, upper = idx_map[idx]
                        dimension = upper - lower
                    if not free_index: dimension = 0
                    for index, order in subexpr.args[1:]:
                        if str(index) in self._namespace['index']:
                            lower, upper = self._namespace['index'][str(index)]
                        else: lower, upper = (0, dimension)
                        if str(index) in idx_map and (lower, upper) != idx_map[str(index)]:
                            raise ParseError('inconsistent indexing range for index \'%s\'' %
                                index, self.lexer.sentence)
                        idx_map[str(index)] = (lower, upper)
                        if index not in self._namespace['basis']:
                            if not self._namespace['basis']:
                                message = 'cannot differentiate symbolically without specifying a basis'
                                raise ParseError(message, self.lexer.sentence)
                            derivative += ', (basis[%s], %s)' % (index, order)
                        else: derivative += ', (%s, %s)' % (index, order)
                    derivative += ')'
                    modified = modified.replace(str(subexpr), derivative)
            modified = replace_function(modified, element, idx_map)
            free_index, bound_index = separate_indexing(modified)
            free_index_RHS.append(free_index)
            # generate implied summation over every bound index
            for idx in bound_index:
                lower, upper = idx_map[idx]
                modified = 'sum(%s for %s in range(%d, %d))' % (modified, idx, lower, upper)
            RHS = RHS.replace(original, modified)
        for i in range(len(free_index_RHS)):
            if sorted(free_index_LHS) != sorted(free_index_RHS[i]):
                # raise exception upon violation of the following rule:
                # a free index must appear in every term with the same
                # position and cannot be summed over in any term
                raise TensorError('unbalanced free index')
        # generate tensor instantiation with implied summation
        for idx, _ in reversed(free_index_LHS):
            lower, upper = idx_map[idx]
            RHS = '[%s for %s in range(%d, %d)]' % (RHS, idx, lower, upper)
            LHS_dimension = upper - lower
        if not free_index_LHS: LHS_dimension = 0
        # shift tensor indexing forward whenever dimension > upper bound
        for subtree in tree.preorder():
            subexpr = subtree.expr
            if subexpr.func == Function('Tensor'):
                symbol = str(subexpr.args[0])
                dimension = self._namespace[symbol].dimension
                tensor = Tensor(subexpr, dimension)
                array_format = tensor.array_format()
                for index in subexpr.args[1:]:
                    if str(index) in self._namespace['index']:
                        _, upper = self._namespace['index'][str(index)]
                        if dimension > upper:
                            shift = dimension - upper
                            for i, (idx, pos) in enumerate(tensor.indexing):
                                if str(idx) == str(index):
                                    tensor.indexing[i] = ('%s + %s' % (idx, shift), pos)
                RHS = RHS.replace(array_format, tensor.array_format())
        if rank == len(re.findall(r'\[[^0-9\]]+\]', LHS)):
            return (LHS.split('[')[0], RHS), LHS_dimension
        LHS_dimension = self._namespace[LHS.split('[')[0]].dimension
        return (re.sub(r'\[[^0-9\]]+\]', '[:]', LHS), RHS), LHS_dimension

    @staticmethod
    def _generate_christoffel(function):
        symbol, indexing = '\\' + str(function.args[0])[:-3], function.args[1:]
        diacritic = 'bar'   if 'bar'   in symbol \
               else 'hat'   if 'hat'   in symbol \
               else 'tilde' if 'tilde' in symbol \
               else None
        metric = '\\%s{g}' % diacritic if diacritic in ('bar', 'hat') \
            else '\\%s{\\gamma}' % diacritic if diacritic == 'tilde' \
            else 'g'
        if diacritic: symbol = '\\%s{%s}' % (diacritic, symbol[:-len(diacritic)])
        indexing = [('\\' if len(str(index)) > 1 else '') + str(index) for index in indexing]
        bound_index = next(x for x in (chr(97 + n) for n in range(26)) if x not in indexing)
        return (('{symbol}^{i1}_{{{i2}{i3}}} = \\frac{{1}}{{2}} {metric}^{{{i1} {bound_index}}}(\\partial_{i2} {metric}_{{{i3} {bound_index}}} + \\partial_{i3} {metric}_{{{bound_index} {i2}}} - \\partial_{bound_index} {metric}_{{{i2} {i3}}})')
                .format(i1 = indexing[0], i2 = indexing[1], i3 = indexing[2], symbol = symbol, metric = metric, bound_index = bound_index))

    @staticmethod
    def _generate_covdrv(tensor, deriv_index, deriv_type, diacritic=''):
        indexing = [str(index[0]) for index in chain(tensor.indexing, deriv_index)]
        alphabet, order, LHS = (chr(97 + n) for n in range(26)), len(deriv_index), ''
        for i, index in enumerate(indexing):
            if index in indexing[:i]:
                indexing[i] = next(x for x in alphabet if x not in indexing)
        for diff_index in indexing[-order:]:
            if len(diff_index) > 1:
                diff_index = '\\' + diff_index
            LHS += ('\\%s{\\nabla}' % diacritic if diacritic else '\\nabla') + ('_%s ' % diff_index)
        LHS += tensor.latex_format()
        def generate_RHS(symbol, order, indexing):
            if order == 0:
                _tensor = Tensor(tensor.function, tensor.dimension)
                _tensor.indexing = [(index, position)
                    for index, (_, position) in zip(indexing, tensor.indexing)]
                return _tensor.latex_format()
            diff_index, RHS = indexing[len(indexing) - order], ''
            if len(diff_index) > 1:
                diff_index = '\\' + diff_index
            latex = generate_RHS(symbol, order - 1, indexing)
            RHS += '\\partial_%s (%s)' % (diff_index, latex)
            for index, (_, position) in zip(indexing, tensor.indexing):
                alphabet = (chr(97 + n) for n in range(26))
                bound_index = next(x for x in alphabet if x not in indexing)
                latex = generate_RHS(symbol, order - 1,
                    [bound_index if i == str(index) else i for i in indexing])
                if len(str(index)) > 1:
                    index = '\\' + str(index)
                RHS += ' + ' if position == 'U' else ' - '
                RHS += '\\%s{\\Gamma}' % diacritic if diacritic else '\\Gamma'
                if position == 'U':
                    RHS += '^%s_{%s %s} (%s)' % (index, bound_index, diff_index, latex)
                else:
                    RHS += '^%s_{%s %s} (%s)' % (bound_index, index, diff_index, latex)
            return RHS
        if deriv_type != 'symbolic':
            LHS = '% define deriv ' + deriv_type + ';\n' + LHS
        return LHS + ' = ' + generate_RHS(tensor.symbol, order, indexing)

    @staticmethod
    def ignore_override():
        warnings.filterwarnings('ignore', category=OverrideWarning)

    @staticmethod
    def clear_namespace():
        Parser._namespace.clear()

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

    def __init__(self, message, sentence, position=None):
        if position is not None:
            length = 0
            for i, substring in enumerate(sentence.split('\n')):
                if position - length <= len(substring):
                    sentence = substring.strip()
                    position -= length + i
                    position += len(sentence) - len(substring)
                    break
                length += len(substring)
            super(ParseError, self).__init__('%s\n%s^\n' % (sentence, (12 + position) * ' ') + message)
        else: super(ParseError, self).__init__(message)

class Tensor:
    """ Tensor Structure """

    def __init__(self, function, dimension, structure=None):
        self.function  = function
        self.dimension = dimension
        self.structure = structure
        self.symbol    = str(function.args[0])
        location = re.findall(r'[UD]', self.symbol)
        self.rank      = len(location)
        self.indexing  = list(zip(function.args[1:], location))

    def array_format(self):
        """ Tensor Notation for Array Formatting """
        if not self.indexing:
            return self.symbol
        return self.symbol + ''.join(['[' + str(index) + ']' for index, _ in self.indexing])

    def latex_format(self):
        """ Tensor Notation for LaTeX Formatting """
        latex = [re.split('[UD]', self.symbol)[0], [], []]
        if len(latex[0]) > 1: latex[0] = '\\' + latex[0]
        U_count, D_count = 0, 0
        for index, position in self.indexing:
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
        if self.rank == 0 and self.dimension == 0:
            return 'Scalar(%s)' % self.symbol
        return 'Tensor(%s, %dD)' % (self.symbol, self.dimension)

    __str__ = __repr__

class TensorError(Exception):
    """ Invalid Tensor Indexing or Dimension """
class OverrideWarning(UserWarning):
    """ Overridden Namespace Variable """

def parse_expr(sentence, verbose=False):
    """ Convert LaTeX Sentence to SymPy Expression (Expression Mode)

        :arg: latex sentence (raw string)
        :arg: verbose mode [default: disabled]
        :return: expression
    """
    return Parser(verbose).parse(sentence, expression=True)

def parse(sentence, verbose=False):
    """ Convert LaTeX Sentence to SymPy Expression

        :arg: latex sentence (raw string)
        :arg: verbose mode [default: disabled]
        :return: namespace
    """
    if not Parser.continue_parsing:
        Parser.clear_namespace()
    _namespace = Parser._namespace.copy()
    namespace = Parser(verbose).parse(sentence)
    kwrd_dict = {}
    for kwrd in ('basis', 'deriv', 'index'):
        kwrd_dict[kwrd] = namespace[kwrd]
        del namespace[kwrd]
    key_diff = tuple(key for key in namespace if key not in _namespace)
    # inject updated namespace into the previous stack frame
    frame = currentframe().f_back
    for key in namespace:
        if isinstance(namespace[key], Tensor):
            frame.f_globals[key] = namespace[key].structure
        elif isinstance(namespace[key], Function('Constant')):
            frame.f_globals[key] = namespace[key].args[0]
        else:
            frame.f_globals[key] = namespace[key]
    namespace.update(kwrd_dict)
    if verbose:
        return tuple(namespace[key] for key in key_diff)
    return key_diff
