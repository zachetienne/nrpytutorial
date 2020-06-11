""" Convert LaTeX Sentence to SymPy Expression """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

# pylint: disable=unused-import
from sympy import Function, Symbol, Integer, Rational, Float, Pow
from sympy import var, sqrt, expand
import re

# pylint: disable=attribute-defined-outside-init
class Lexer:
    """ LaTeX Lexer

        The following class will tokenize a sentence for parsing.
    """

    def __init__(self, namespace):
        greek  = '|'.join([r'\\' + i for i in ('[aA]lpha', '[bB]eta', '[gG]amma', '[dD]elta',
            '[eE]psilon', '[zZ]eta', '[eE]ta', '[tT]heta', '[iI]ota', '[kK]appa', '[lL]ambda',
            '[mM]u', '[nN]u', '[xX]i', '[oO]mikron', '[pP]i', '[Rr]ho', '[sS]igma', '[tT]au',
            '[uU]psilon', '[pP]hi', '[cC]hi', '[pP]si', '[oO]mega')])
        namespace = '|'.join(map(re.escape, ('\\' + i if len(i) > 1 else i
            for i in re.split(r'[\s,]+', namespace)))) if namespace else '(?!)'
        self.regex = re.compile('|'.join(['(?P<%s>%s)' % pattern for pattern in
            [ ('RATIONAL',       r'[0-9]+\/[1-9]+|\\frac{[0-9]+}{[1-9]+}'),
              ('DECIMAL',        r'[0-9]+\.[0-9]+'),
              ('INTEGER',        r'[0-9]+'),
              ('PLUS',           r'\+'),
              ('MINUS',          r'\-'),
              ('DIVIDE',         r'\/'),
              ('EQUAL',          r'\='),
              ('CARET',          r'\^'),
              ('UNDERSCORE',     r'\_'),
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
              ('SPACE_DELIM',    r'(?:\s|\\,)+'),
              ('SQRT_CMD',       r'\\sqrt'),
              ('FRAC_CMD',       r'\\frac'),
              ('TENSOR',         namespace),
              ('SYMBOL',         greek + r'|[a-zA-Z]'),
              ('COMMAND',        r'\\[a-z]+')]]))

    def initialize(self, sentence):
        """ Initialize Lexer

            :arg: sentence (raw string)
        """
        self.sentence = sentence
        self.token    = None
        self.lexeme   = None
        self.index    = 0

    def tokenize(self):
        """ Tokenize Sentence

            :return: token iterator
        """
        while self.index < len(self.sentence):
            token = self.regex.match(self.sentence, self.index)
            if token is None:
                raise ParseError('unexpected \'%s\' at position %d' %
                    (self.sentence[self.index], self.index), self.sentence, self.index)
            if token.lastgroup in ('BIGL_DELIM', 'BIGR_DELIM',
                    'LEFT_DELIM', 'RIGHT_DELIM', 'SPACE_DELIM'):
                self.index = token.end()
            else:
                self.index  = token.end()
                self.lexeme = token.group()
                yield token.lastgroup

    def lex(self):
        """ Retrieve Next Token

            :return: next token in iterator
        """
        try:
            self.token = next(self.tokenize())
        except StopIteration:
            self.token = None
        return self.token

    def reset(self):
        """ Reset Token Iterator """
        if not self.sentence:
            raise RuntimeError('cannot reset uninitialized lexer')
        self.initialize(self.sentence)

# pylint: disable=anomalous-backslash-in-string
class Parser:
    """ LaTeX Parser

        The following class will parse a tokenized sentence.

        LaTeX Grammar:
        <ROOT>          -> <VARIABLE> = <EXPR> | <EXPR>
        <EXPR>          -> [ - ] <TERM> { ( + | - ) <TERM> }
        <TERM>          -> <FACTOR> { [ / ] <FACTOR> }
        <FACTOR>        -> <SUBEXPR> { ^( <SUBEXPR> | {<EXPR>} ) }
        <SUBEXPR>       -> <OPERAND> | (<EXPR>) | [<EXPR>]
        <OPERAND>       -> <VARIABLE> | <NUMBER> | <COMMAND>
        <VARIABLE>      -> <ARRAY> | <SYMBOL> [ _( <SYMBOL> | <INTEGER> ) ]
        <NUMBER>        -> <RATIONAL> | <DECIMAL> | <INTEGER>
        <COMMAND>       -> <SQRT> | <FRAC>
        <SQRT>          -> \ sqrt [ [<INTEGER>] ] {<EXPR>}
        <FRAC>          -> \ frac {<EXPR>} {<EXPR>}
        <ARRAY>         -> <TENSOR> [ _( <SYMBOL> | {{ <SYMBOL> }} ) [ ^( <SYMBOL> | {{ <SYMBOL> }} ) ]
                            | ^( <SYMBOL> | {{ <SYMBOL> }} ) [ _( <SYMBOL> | {{ <SYMBOL> }} ) ] ]
    """

    def __init__(self, namespace):
        self.lexer = Lexer(namespace)
        self.evaluate = namespace is None

    def parse(self, sentence):
        """ Parse Sentence

            :arg:    sentence (raw string)
            :return: symbolic expression
        """
        self.lexer.initialize(sentence)
        self.lexer.lex()
        root = self.__root()
        if self.lexer.token:
            sentence = self.lexer.sentence
            position = self.lexer.index - len(self.lexer.lexeme)
            raise ParseError('unexpected \'%s\' at position %d' %
                (sentence[position], position), sentence, position)
        return root

    def __root(self):
        if self.peek('SYMBOL') or self.peek('TENSOR'):
            if self.peek('SYMBOL') and not self.evaluate:
                self.expect('TENSOR')
            variable = self.__variable()
            if self.accept('EQUAL'):
                expr = self.__expr()
                if self.evaluate:
                    return {variable[8:-2]: eval(expr)}
                variable = str(eval(variable))
                expr = str(expand(eval(expr)))
                args = variable[7:-1].split(', ')
                suffix = ''.join(['[' + n + ']' for n in args[1:]])
                variable = args[0] + suffix
                for match in re.findall(r'Tensor\([^\)]+\)', expr):
                    args = match[7:-1].split(', ')
                    suffix = ''.join(['[' + n + ']' for n in args[1:]])
                    expr = expr.replace(match, args[0] + suffix)
                return {variable: expr}
            self.lexer.reset()
            self.lexer.lex()
        if self.evaluate:
            return eval(self.__expr())
        return self.__expr()

    def __expr(self):
        sign = '-' if self.accept('MINUS') else ''
        expr = sign + self.__term()
        while self.peek('PLUS') or self.peek('MINUS'):
            operator = self.lexer.lexeme
            self.lexer.lex()
            expr += ' %s %s' % (operator, self.__term())
        return expr

    def __term(self):
        expr = self.__factor()
        while any(self.peek(i) for i in ('LEFT_PAREN', 'LEFT_BRACKET',
                'SYMBOL', 'TENSOR', 'RATIONAL', 'DECIMAL', 'INTEGER', 'DIVIDE',
                'SQRT_CMD', 'FRAC_CMD')):
            operator = '/' if self.accept('DIVIDE') else '*'
            expr += '%s%s' % (operator, self.__factor())
        return expr

    def __factor(self):
        expr = self.__subexpr()
        while self.accept('CARET'):
            if self.accept('LEFT_BRACE'):
                expr += '**(' + self.__expr() + ')'
                self.expect('RIGHT_BRACE')
            else: expr += '**' + self.__subexpr()
        return expr

    def __subexpr(self):
        if self.accept('LEFT_PAREN'):
            expr = '(' + self.__expr() + ')'
            self.expect('RIGHT_PAREN')
            return expr
        if self.accept('LEFT_BRACKET'):
            expr = '(' + self.__expr() + ')'
            self.expect('RIGHT_BRACKET')
            return expr
        return self.__operand()

    def __operand(self):
        if self.peek('SYMBOL') or self.peek('TENSOR'):
            return self.__variable()
        if any(self.peek(i) for i in ('SQRT_CMD', 'FRAC_CMD', 'COMMAND')):
            return self.__command()
        return self.__number()

    def __variable(self):
        variable = self.lexer.lexeme
        if variable[0] == '\\':
            variable = variable[1:]
        if self.accept('SYMBOL'):
            if self.accept('UNDERSCORE'):
                if self.peek('SYMBOL') or self.peek('INTEGER'):
                    subscript = self.lexer.lexeme
                    if subscript[0] == '\\':
                        subscript = subscript[1:]
                    variable += '_' + subscript
                    self.lexer.lex(); var(variable)
                    return 'Symbol(\'' + variable + '\')'
                sentence = self.lexer.sentence
                position = self.lexer.index - len(self.lexer.lexeme)
                raise ParseError('unexpected \'%s\' at position %d' %
                    (sentence[position], position), sentence, position)
            var(variable)
            return 'Symbol(\'' + variable + '\')'
        return self.__array()

    def __array(self):
        array = self.lexer.lexeme
        if array[0] == '\\':
            array = array[1:]
        self.lexer.lex()
        array, index = list(array), []
        def add_index(suffix):
            symbol = self.lexer.lexeme
            if symbol[0] == '\\':
                symbol = symbol[1:]
            var(symbol)
            self.lexer.lex()
            array.append(suffix)
            index.append(symbol)
        if self.accept('UNDERSCORE'):
            if self.peek('SYMBOL'):
                add_index('D')
            elif self.accept('LEFT_BRACE'):
                while self.peek('SYMBOL'):
                    add_index('D')
                self.expect('RIGHT_BRACE')
            if self.accept('CARET'):
                if self.peek('SYMBOL'):
                    add_index('U')
                elif self.accept('LEFT_BRACE'):
                    while self.peek('SYMBOL'):
                        add_index('U')
                    self.expect('RIGHT_BRACE')
        elif self.accept('CARET'):
            if self.peek('SYMBOL'):
                add_index('U')
            elif self.accept('LEFT_BRACE'):
                while self.peek('SYMBOL'):
                    add_index('U')
                self.expect('RIGHT_BRACE')
            if self.accept('UNDERSCORE'):
                if self.peek('SYMBOL'):
                    add_index('D')
                elif self.accept('LEFT_BRACE'):
                    while self.peek('SYMBOL'):
                        add_index('D')
                    self.expect('RIGHT_BRACE')
        array = ''.join(array); var(array)
        if not index:
            return 'Function(\'Tensor\')(' + array + ')'
        return 'Function(\'Tensor\')(%s, %s)' % (array, ', '.join(index))

    def __number(self):
        number = self.lexer.lexeme
        if self.accept('RATIONAL'):
            rational = re.match(r'([1-9][0-9]*)\/([1-9][0-9]*)', number)
            if not rational:
                rational = re.match(r'\\frac{([0-9]+)}{([1-9]+)}', number)
            return 'Rational(%s, %s)' % (rational.group(1), rational.group(2))
        if self.accept('DECIMAL'):
            return 'Float(' + number + ')'
        if self.accept('INTEGER'):
            return 'Integer(' + number + ')'
        sentence = self.lexer.sentence
        position = self.lexer.index - len(self.lexer.lexeme)
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    def __command(self):
        command = self.lexer.lexeme
        if self.accept('SQRT_CMD'):
            return self.__sqrt()
        if self.accept('FRAC_CMD'):
            return self.__frac()
        position = self.lexer.index - len(self.lexer.lexeme)
        raise ParseError('unsupported command \'%s\' at position %d' %
            (command, position), self.lexer.sentence, position)

    def __sqrt(self):
        if self.accept('LEFT_BRACKET'):
            root = 'Integer(' + self.lexer.lexeme + ')'
            self.expect('INTEGER')
            self.expect('RIGHT_BRACKET')
        else: root = 2
        self.expect('LEFT_BRACE')
        expr = self.__expr()
        self.expect('RIGHT_BRACE')
        return 'Pow(%s, Rational(1, %s))' % (expr, root)

    def __frac(self):
        self.expect('LEFT_BRACE')
        numerator = self.__expr()
        self.expect('RIGHT_BRACE')
        self.expect('LEFT_BRACE')
        denominator = self.__expr()
        self.expect('RIGHT_BRACE')
        return '(%s)/(%s)' % (numerator, denominator)

    def peek(self, token_type):
        return self.lexer.token == token_type

    def accept(self, token_type):
        if self.peek(token_type):
            self.lexer.lex()
            return True
        return False

    def expect(self, token_type):
        if not self.accept(token_type):
            position = self.lexer.index - len(self.lexer.lexeme)
            raise ParseError('expected token %s at position %d' %
                (token_type, position), self.lexer.sentence, position)

class ParseError(Exception):
    """ Invalid LaTeX Sentence """

    def __init__(self, message, sentence, position):
        super().__init__('%s\n%s^\n' % (sentence, (12 + position) * ' ') + message)

def __summation(equation, dimension):
    for var_ in equation:
        var, expr = var_, equation[var_]
    LHS, RHS = list(zip(re.findall(r'\[([a-zA-Z]+)\]', var),
        re.findall(r'[UD]', var))), []
    loop_count = 0
    for product in re.split(r'\s[\+\-]\s', expr):
        loop_index = (2 * chr(97 + n) for n in range(26))
        idx_lst = re.findall(r'\[([a-zA-Z]+)\]', product)
        pos_lst = re.findall(r'[UD]', product)
        bound_count, index_map = 0, {}
        for idx in set(idx_lst):
            count = U = D = 0
            free_index = []
            for idx_, pos_ in zip(idx_lst, pos_lst):
                if idx_ == idx:
                    free_index.append((idx_, pos_))
                    if pos_ == 'U': U += 1
                    if pos_ == 'D': D += 1
                    count += 1
            if count > 1:
                if count % 2 != 0 or U != D:
                    raise TensorError('illegal bound index')
                lp_idx = next(loop_index)
                if bound_count >= loop_count:
                    index_map[idx] = lp_idx
                else: expr = expr.replace(idx, lp_idx)
                bound_count += 1
            else:
                RHS.extend(free_index)
            if bound_count > loop_count: loop_count = bound_count
        for idx in index_map:
            expr = expr.replace(idx, index_map[idx])
        for idx in index_map:
            expr = 'sum(%s for %s in range(%d))' % \
                (expr, index_map[idx], dimension)
    if LHS:
        if set(LHS) != set(RHS):
            raise TensorError('unbalanced free index')
        for idx, _ in LHS:
            expr = '[%s for %s in range(%d)]' % (expr, idx, dimension)
    return var.split('[')[0], expr

# pylint: disable=missing-class-docstring
class TensorError(Exception): pass

def parse(sentence, namespace=None):
    """ Convert LaTeX Sentence to SymPy Expression

        :arg:    LaTeX Sentence
        :arg:    Tensor Namespace
        :return: SymPy Expression

        >>> from latex_parser import parse
        >>> parse(r'-(x\\frac{2}{3} + \\sqrt[5]{x + 3})')
        -2*x/3 - (x + 3)**(1/5)
        >>> parse(r'x_n = \\sqrt[5]{x + 3}')
        {'x_n': (x + 3)**(1/5)}
        >>> parse(r'v^\\mu = g^{\\mu\\nu}v_\\nu', 'v g')
        {'vU[mu]': 'vD[nu]*gUU[mu][nu]'}

        >>> import indexedexp as ixp
        >>> hUU = ixp.declarerank2('hUU', 'nosym', 2)
        >>> hDD = ixp.declarerank2('hDD', 'nosym', 2)
        >>> namespace = {'h': 0, 'hUU': hUU, 'hDD': hDD}
        >>> parse(r'h = h_{\\mu\\mu}h^{\\mu\\mu}', namespace)
        {'h': hDD00*hUU00 + hDD11*hUU11}
    """
    if namespace is None or isinstance(namespace, str):
        return Parser(namespace).parse(sentence)
    if isinstance(namespace, dict):
        dim_list = [len(namespace[tensor])
            for tensor in namespace if isinstance(namespace[tensor], list)]
        if any(dimension != dim_list[0] for dimension in dim_list):
            raise TensorError('inconsistent tensor dimension')
        equation = Parser(' '.join(re.split(r'[UD]', tensor)[0]
            for tensor in namespace)).parse(sentence)
        var, expr = __summation(equation, dim_list[0])
        globals().update(namespace)
        return {var: eval(expr)}
    raise TypeError('inappropriate type for tensor namespace')

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    