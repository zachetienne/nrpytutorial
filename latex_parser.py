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
        # split namespace string at whitespace or a comma, append a backslash to the
        # front of every multi-character tensor, and join a tensor list together with
        # a pipe symbol for regex alternation
        namespace = '|'.join(map(re.escape, ('\\' + i if len(i) > 1 else i
            for i in re.split(r'[\s,]+', namespace)))) if namespace else '(?!)'
        # define a regex pattern for every token, create a named capture group for
        # every pattern, join a pattern list together with a pipe symbol for regex
        # alternation, and compile the resulting regex expression
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
            # raise exception whenever lexer has a remaining token
            sentence = self.lexer.sentence
            position = self.lexer.index - len(self.lexer.lexeme)
            raise ParseError('unexpected \'%s\' at position %d' %
                (sentence[position], position), sentence, position)
        return root

    # <ROOT> -> <VARIABLE> = <EXPR> | <EXPR>
    def __root(self):
        if self.peek('SYMBOL') or self.peek('TENSOR'):
            # expect a tensor variable on LHS of a tensorial equation
            if self.peek('SYMBOL') and not self.evaluate:
                self.expect('TENSOR')
            variable = self.__variable()
            if self.accept('EQUAL'):
                expr = self.__expr()
                # return a dictionary mapping for a non-tensorial equation
                if self.evaluate:
                    return {variable[8:-2]: eval(expr)}
                # distribute over every parenthetical expression
                expr = str(expand(eval(expr)))
                # remove variable cruft and extract tensor information
                args = variable[19:-1].split(', ')
                # construct tensor indexing for array notation
                suffix = ''.join(['[' + n + ']' for n in args[1:]])
                # construct array notation for a tensor variable
                variable = args[0] + suffix
                # iterate over every tensor in the expression
                for match in re.findall(r'Tensor\([^\)]+\)', expr):
                    # remove cruft and extract tensor information
                    args = match[7:-1].split(', ')
                    # construct tensor indexing for array notation
                    suffix = ''.join(['[' + n + ']' for n in args[1:]])
                    # replace function notation with array notation
                    expr = expr.replace(match, args[0] + suffix)
                    # return a dictionary mapping for a tensorial equation
                return {variable: expr}
            # reset the token iterator for expression parsing
            self.lexer.reset()
            self.lexer.lex()
        if self.evaluate:
            return eval(self.__expr())
        return self.__expr()

    # <EXPR> -> [ - ] <TERM> { ( + | - ) <TERM> }
    def __expr(self):
        sign = '-' if self.accept('MINUS') else ''
        expr = sign + self.__term()
        while self.peek('PLUS') or self.peek('MINUS'):
            operator = self.lexer.lexeme
            self.lexer.lex()
            expr += ' %s %s' % (operator, self.__term())
        return expr

    # <TERM> -> <FACTOR> { [ / ] <FACTOR> }
    def __term(self):
        expr = self.__factor()
        while any(self.peek(i) for i in ('LEFT_PAREN', 'LEFT_BRACKET',
                'SYMBOL', 'TENSOR', 'RATIONAL', 'DECIMAL', 'INTEGER', 'DIVIDE',
                'SQRT_CMD', 'FRAC_CMD')):
            operator = '/' if self.accept('DIVIDE') else '*'
            expr += '%s%s' % (operator, self.__factor())
        return expr

    # <FACTOR> -> <SUBEXPR> { ^( <SUBEXPR> | {<EXPR>} ) }
    def __factor(self):
        expr = self.__subexpr()
        while self.accept('CARET'):
            if self.accept('LEFT_BRACE'):
                expr += '**(' + self.__expr() + ')'
                self.expect('RIGHT_BRACE')
            else: expr += '**' + self.__subexpr()
        return expr

    # <SUBEXPR> -> <OPERAND> | (<EXPR>) | [<EXPR>]
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

    # <OPERAND> -> <VARIABLE> | <NUMBER> | <COMMAND>
    def __operand(self):
        if self.peek('SYMBOL') or self.peek('TENSOR'):
            return self.__variable()
        if any(self.peek(i) for i in ('SQRT_CMD', 'FRAC_CMD', 'COMMAND')):
            return self.__command()
        return self.__number()

    # <VARIABLE> -> <ARRAY> | <SYMBOL> [ _( <SYMBOL> | <INTEGER> ) ]
    def __variable(self):
        variable = self.lexer.lexeme
        # remove backslash from variable whenever present
        if variable[0] == '\\':
            variable = variable[1:]
        if self.accept('SYMBOL'):
            if self.accept('UNDERSCORE'):
                if self.peek('SYMBOL') or self.peek('INTEGER'):
                    subscript = self.lexer.lexeme
                    # remove backslash from subscript whenever present
                    if subscript[0] == '\\':
                        subscript = subscript[1:]
                    variable += '_' + subscript
                    self.lexer.lex()
                    # create variable symbol and insert into global namespace
                    var(variable)
                    return 'Symbol(\'' + variable + '\')'
                # raise exception whenever illegal subscript
                sentence = self.lexer.sentence
                position = self.lexer.index - len(self.lexer.lexeme)
                raise ParseError('unexpected \'%s\' at position %d' %
                    (sentence[position], position), sentence, position)
            # create variable symbol and insert into global namespace
            var(variable)
            return 'Symbol(\'' + variable + '\')'
        return self.__array()

    # <ARRAY> -> <TENSOR> [ _( <SYMBOL> | {{ <SYMBOL> }} ) [ ^( <SYMBOL> | {{ <SYMBOL> }} ) ]
    #             | ^( <SYMBOL> | {{ <SYMBOL> }} ) [ _( <SYMBOL> | {{ <SYMBOL> }} ) ] ]
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

    # <NUMBER> -> <RATIONAL> | <DECIMAL> | <INTEGER>
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

    # <COMMAND> -> <SQRT> | <FRAC>
    def __command(self):
        command = self.lexer.lexeme
        if self.accept('SQRT_CMD'):
            return self.__sqrt()
        if self.accept('FRAC_CMD'):
            return self.__frac()
        position = self.lexer.index - len(self.lexer.lexeme)
        raise ParseError('unsupported command \'%s\' at position %d' %
            (command, position), self.lexer.sentence, position)

    # <SQRT> -> \ sqrt [ [<INTEGER>] ] {<EXPR>}
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

    # <FRAC> -> \ frac {<EXPR>} {<EXPR>}
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
    # construct a tuple list of every LHS free index
    LHS, RHS = set(zip(re.findall(r'\[([a-zA-Z]+)\]', var),
        re.findall(r'[UD]', var))), []
    # iterate over every subexpression containing a product
    for product in re.split(r'\s[\+\-]\s', expr):
        # generate a loop index iterator (ex: aa, bb, ...)
        loop_index = (2 * chr(97 + n) for n in range(26))
        # extract every index present in the subexpression
        idx_lst = re.findall(r'\[([a-zA-Z]+)\]', product)
        # extract every index position (ex: U or D)
        pos_lst = re.findall(r'[UD]', product)
        free_index, bound_index = [], {}
        # iterate over every unique index in the subexpression
        for idx in set(idx_lst):
            count = U = D = 0; index_tuple = []
            # count index occurrence and position occurrence
            for idx_, pos_ in zip(idx_lst, pos_lst):
                if idx_ == idx:
                    index_tuple.append((idx_, pos_))
                    if pos_ == 'U': U += 1
                    if pos_ == 'D': D += 1
                    count += 1
            # assign mapping: bound index -> loop index
            if count > 1:
                if count % 2 != 0 or U != D:
                    # raise exception upon violation of the following rule:
                    # a bound index must appear exactly once as a superscript
                    # and exactly once as a subscript in any single term
                    raise TensorError('illegal bound index')
                bound_index[idx] = next(loop_index)
            # identify every free index on the RHS
            else: free_index.extend(index_tuple)
        RHS.append(set(free_index))
        summation = product
        # replace every bound index with a loop index
        for idx in bound_index:
            summation = summation.replace(idx, bound_index[idx])
        # generate implied summation over every bound index
        for idx in bound_index:
            summation = 'sum([%s for %s in range(%d)])' % \
                (summation, bound_index[idx], dimension)
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
    return var.split('[')[0], expr

# pylint: disable=missing-class-docstring
class TensorError(Exception): pass

def parse(sentence, namespace=None, debug=False):
    """ Convert LaTeX Sentence to SymPy Expression

        :arg:    LaTeX Sentence
        :arg:    Tensor Namespace
        :arg:    Debug Mode
        :return: SymPy Expression

        >>> from latex_parser import parse
        >>> parse(r'-(x\\frac{2}{3} + \\sqrt[5]{x + 3})')
        -2*x/3 - (x + 3)**(1/5)
        >>> parse(r'x_n = \\sqrt[5]{x + 3}')
        {'x_n': (x + 3)**(1/5)}
        >>> parse(r'v^\\mu = g^{\\mu\\nu}v_\\nu', 'v g')
        {'vU[mu]': 'vD[nu]*gUU[mu][nu]'}

        >>> import indexedexp as ixp
        >>> hUU = ixp.declarerank2('hUU', 'sym01', 2)
        >>> hDD = ixp.declarerank2('hDD', 'sym01', 2)
        >>> namespace = {'h': 0, 'hUU': hUU, 'hDD': hDD}
        >>> parse(r'h = h_{\\mu\\mu}h^{\\mu\\mu}', namespace)
        {'h': hDD00*hUU00 + hDD11*hUU11}

        >>> vD  = ixp.declarerank1('vD', DIM=2)
        >>> wD  = ixp.declarerank1('wD', DIM=2)
        >>> gUU = ixp.declarerank2('gUU', 'sym01', DIM=2)
        >>> namespace = {'vD': vD, 'wD': wD, 'gUU': gUU}
        >>> parse(r'v^i = g^{ij}(v_j + w_j)', namespace)
        {'vU': [gUU00*vD0 + gUU00*wD0 + gUU01*vD1 + gUU01*wD1, \
gUU01*vD0 + gUU01*wD0 + gUU11*vD1 + gUU11*wD1]}
    """
    if namespace is None or isinstance(namespace, str):
        return Parser(namespace).parse(sentence)
    if isinstance(namespace, dict):
        # construct a list for every tensor dimension
        dim_list = [len(namespace[tensor])
            for tensor in namespace if isinstance(namespace[tensor], list)]
        if any(dimension != dim_list[0] for dimension in dim_list):
            raise TensorError('inconsistent tensor dimension')
        # extract every tensor from namespace and parse sentence
        equation = Parser(' '.join(re.split(r'[UD]', tensor)[0]
            for tensor in namespace)).parse(sentence)
        # perform implied summation on the tensorial equation
        var, expr = __summation(equation, dim_list[0])
        # update the global namespace with the tensor namespace
        globals().update(namespace)
        if debug: return {var: expr}
        return {var: eval(expr)}
    raise TypeError('inappropriate type for tensor namespace')

if __name__ == "__main__":
    import doctest
    doctest.testmod()
