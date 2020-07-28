""" Convert LaTeX Sentence to SymPy Expression """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

# pylint: disable=attribute-defined-outside-init
from sympy import Function, Derivative, Symbol, Integer, Rational, Float, Pow
from sympy import expand, diff, var
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
        # extract every tensor from namespace, append a backslash to the front of every
        # multi-letter tensor, and join together the resulting tensor list using a pipe
        # symbol for regex alternation
        nameset = set(re.match(r'[^UD]*', tensor).group() for tensor in namespace)
        array_pattern = '|'.join(r'\\' + name if len(name) > 1 else name
            for name in nameset) if nameset else '(?!)'
        greek_pattern = '|'.join(r'\\' + letter for letter in ('[aA]lpha', '[bB]eta', '[gG]amma', '[dD]elta',
            '[eE]psilon', '[zZ]eta', '[eE]ta', '[tT]heta', '[iI]ota', '[kK]appa', '[lL]ambda',
            '[mM]u', '[nN]u', '[xX]i', '[oO]mikron', '[pP]i', '[Rr]ho', '[sS]igma', '[tT]au',
            '[uU]psilon', '[pP]hi', '[cC]hi', '[pP]si', '[oO]mega'))
        # define a regex pattern for every token, create a named capture group for
        # every pattern, join together the resulting pattern list using a pipe symbol
        # for regex alternation, and compile the generated regular expression
        self.regex = re.compile('|'.join(['(?P<%s>%s)' % pattern for pattern in
            [ ('SPACE_DELIM',    r'(?:\s|\\,|\{\})+|\&'),
              ('RATIONAL',       r'[0-9]+\/[1-9]+|\\frac{[0-9]+}{[1-9]+}'),
              ('DECIMAL',        r'[0-9]+\.[0-9]+'),
              ('INTEGER',        r'[0-9]+'),
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
              ('DERV_CMD',       r'\\partial'),
              ('SYMMETRY',       r'nosym|sym[0-9]+(?:_sym[0-9]+)*|metric|permutation'),
              ('ARRAY',          array_pattern),
              ('SYMBOL',         greek_pattern + r'|[a-zA-Z]'),
              ('COMMAND',        r'\\[a-z]+')]]))

    def initialize(self, sentence, position=0):
        """ Initialize Lexer

            :arg: sentence (raw string)
            :arg: position [default: 0]
        """
        self.sentence = sentence
        self.token    = None
        self.lexeme   = None
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

            :return: next token in iterator
        """
        try:
            self.token = next(self.tokenize())
        except StopIteration:
            self.token = None
        return self.token

    def reset(self, position=0):
        """ Reset Token Iterator

            :arg: reset position [default: 0]
        """
        if not self.sentence:
            raise RuntimeError('cannot reset uninitialized lexer')
        self.initialize(self.sentence, position)

    def update(self, namespace):
        """ Update Tensor Pattern

            :arg: tensor namespace
        """
        nameset = set(re.match(r'[^UD]*', tensor).group() for tensor in namespace)
        tensor_pattern = '|'.join(r'\\\\' + name if len(name) > 1 else name
            for name in nameset) if nameset else '(?!)'
        self.regex = re.compile(re.sub(r'<ARRAY>.+?(?=\)\|)',
            '<ARRAY>' + tensor_pattern, self.regex.pattern))

class Parser:
    """ LaTeX Parser

        The following class will parse a tokenized LaTeX sentence.

        LaTeX Extended BNF Grammar:
        <ROOT>          -> <EXPRESSION> | <STRUCTURE> { <LINE_BREAK> <STRUCTURE> }*
        <STRUCTURE>     -> <CONFIG> | <ENVIRONMENT> | <ASSIGNMENT>
        <ENVIRONMENT>   -> <BEGIN_ALIGN> <ASSIGNMENT> { <LINE_BREAK> <ASSIGNMENT> }* <END_ALIGN>
        <ASSIGNMENT>    -> <VARIABLE> = <EXPRESSION>
        <EXPRESSION>    -> <TERM> { ( '+' | '-' ) <TERM> }*
        <TERM>          -> <FACTOR> { [ '/' ] <FACTOR> }*
        <FACTOR>        -> <BASE> { '^' <EXPONENT> }*
        <BASE>          -> [ '-' ] ( <ATOM> | '(' <EXPRESSION> ')' [ '_' '{' ',' { <SYMBOL> }+ '}' ] | '[' <EXPRESSION> ']' )
        <EXPONENT>      -> <BASE> | '{' <BASE> '}'
        <ATOM>          -> <VARIABLE> | <NUMBER> | <COMMAND>
        <VARIABLE>      -> <TENSOR> | <SYMBOL> [ '_' ( <SYMBOL> | <INTEGER> ) ]
        <NUMBER>        -> <RATIONAL> | <DECIMAL> | <INTEGER>
        <COMMAND>       -> <SQRT> | <FRAC> | <DERV>
        <SQRT>          -> '\\sqrt' [ '[' <INTEGER> ']' ] '{' <EXPRESSION> '}'
        <FRAC>          -> '\\frac' '{' <EXPRESSION> '}' '{' <EXPRESSION> '}'
        <DERV>          -> { '\\partial' '_' <SYMBOL> }+ ( <TENSOR> | '(' <EXPRESSION> ')' )
        <CONFIG>        -> '%' <TENSOR> '[' <INTEGER> ']' [ ':' <SYMMETRY> ] { ',' <TENSOR> '[' <INTEGER> ']' [ ':' <SYMMETRY> ] }*
        <TENSOR>        -> <ARRAY> ( '_' <LOWER_INDEX> ) | ( '^' <UPPER_INDEX> [ '_' <LOWER_INDEX> ] ] )
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
            :arg:    expression mode [default: False]
            :return: symbolic expression or namespace tuple
        """
        self.lexer.initialize(sentence)
        self.lexer.lex()
        if expression:
            return self._root(expression)
        self._root(expression)
        for var in self.namespace:
            if var in Parser.namespace:
                warn('\'' + var + '\'', OverrideWarning, stacklevel=3)
            if isinstance(self.namespace[var], Tensor):
                self.namespace[var] = self.namespace[var].array
            Parser.namespace[var] = self.namespace[var]
        return self.namespace

    # <ROOT> -> <EXPRESSION> | <STRUCTURE> { <LINE_BREAK> <STRUCTURE> }*
    def _root(self, expression):
        if expression: return self._expression()
        self._structure()
        while self.accept('LINE_BREAK'):
            self._structure()
        return None

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

    # <ASSIGNMENT> -> <VARIABLE> = <EXPRESSION>
    def _assignment(self):
        tensorial = self._tensorial()
        # expect a tensor on LHS of a tensorial equation
        if tensorial and self.peek('SYMBOL'):
            self.lexer.token = 'ARRAY'
        variable = self._variable()
        self.expect('EQUAL')
        expr = self._expression()
        if self.namespace:
            # distribute over every parenthetical expression
            expr_ = expand(expr); expr = str(expr_)
            variable = ''.join(Tensor.notation(variable))
            # iterate over every tensor in the expression
            for subtree in ExprTree(expr_).preorder():
                subexpr = subtree.expr
                if subexpr.func == Function('Tensor'):
                    # replace function notation with array notation
                    func = ''.join(Tensor.notation(subexpr))
                    expr = expr.replace(str(subexpr), func)
        if tensorial:
            # perform implied summation on a tensorial equation
            self.namespace.update(_summation((str(variable), expr), self.dimension))
            self.lexer.update(self.namespace)
        else: self.namespace[str(variable)] = expr

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
        while any(self.peek(i) for i in ('LEFT_PAREN', 'LEFT_BRACKET',
                'SYMBOL', 'ARRAY', 'RATIONAL', 'DECIMAL', 'INTEGER',
                'DIVIDE', 'SQRT_CMD', 'FRAC_CMD', 'DERV_CMD')):
            if self.accept('DIVIDE'):
                expr /= self._factor()
            else: expr *= self._factor()
        return expr

    # <FACTOR> -> <BASE> { '^' <EXPONENT> }*
    def _factor(self):
        stack = [self._base()]
        while self.accept('CARET'):
            stack.append(self._exponent())
        expr = stack.pop()
        for subexpr in reversed(stack):
            expr = subexpr ** expr
        return expr

    # <BASE> -> [ '-' ] ( <ATOM> | '(' <EXPRESSION> ')' [ '_' '{' ',' { <SYMBOL> }+ '}' ] | '[' <EXPRESSION> ']' )
    def _base(self):
        sign = -1 if self.accept('MINUS') else 1
        if self.accept('LEFT_PAREN'):
            expr = sign * self._expression()
            self.expect('RIGHT_PAREN')
            if self.accept('UNDERSCORE'):
                self.expect('LEFT_BRACE')
                self.expect('COMMA')
                indices = [self.lexer.lexeme]
                self.expect('SYMBOL')
                while self.peek('SYMBOL'):
                    index = self.lexer.lexeme
                    self.expect('SYMBOL')
                    indices.extend(index)
                self.expect('RIGHT_BRACE')
                tree = ExprTree(expr)
                for subtree in tree.preorder():
                    subexpr = subtree.expr
                    if subexpr.func == Function('Tensor'):
                        subtree.expr = Function('Tensor')(Symbol('_x'), *subexpr.args)
                        tree.build(subtree)
                expr = tree.reconstruct()
                tree = ExprTree(diff(expr, Symbol('_x')))
                for subtree in tree.preorder():
                    subexpr = subtree.expr
                    if subexpr.func == Derivative:
                        order, array = len(indices), subexpr.args[1]
                        indexing = list(array.args[1:-1]) + indices
                        tensor = self.namespace[str(array.args[0])]
                        subtree.expr = self._differentiate(tensor, order, *indexing)
                        tree.build(subtree)
                    elif subexpr.func == Function('Tensor'):
                        subtree.expr = Function('Tensor')(*subexpr.args[1:])
                        tree.build(subtree)
                return tree.reconstruct()
            return expr
        if self.accept('LEFT_BRACKET'):
            expr = sign * self._expression()
            self.expect('RIGHT_BRACKET')
            return expr
        return sign * self._atom()

    # <EXPONENT> -> <BASE> | '{' <BASE> '}'
    def _exponent(self):
        if self.accept('LEFT_BRACE'):
            base = self._base()
            self.expect('RIGHT_BRACE')
            return base
        return self._base()

    # <ATOM> -> <VARIABLE> | <NUMBER> | <COMMAND>
    def _atom(self):
        if self.peek('SYMBOL') or self.peek('ARRAY'):
            return self._variable()
        if any(self.peek(i) for i in ('COMMAND', 'SQRT_CMD', 'FRAC_CMD', 'DERV_CMD')):
            return self._command()
        return self._number()

    # <VARIABLE> -> <TENSOR> | <SYMBOL> [ '_' ( <SYMBOL> | <INTEGER> ) ]
    def _variable(self):
        if self.peek('SYMBOL'):
            variable = self.lexer.lexeme
            self.lexer.lex()
            # remove backslash from variable whenever present
            if variable[0] == '\\':
                variable = variable[1:]
            if self.accept('UNDERSCORE'):
                if self.peek('SYMBOL') or self.peek('INTEGER'):
                    subscript = self.lexer.lexeme
                    self.lexer.lex()
                    # remove backslash from subscript whenever present
                    if subscript[0] == '\\':
                        subscript = subscript[1:]
                    variable += '_' + subscript
                    var(variable)
                    return Symbol(variable)
                # raise exception whenever illegal subscript
                sentence = self.lexer.sentence
                position = self.lexer.index - len(self.lexer.lexeme)
                raise ParseError('unexpected \'%s\' at position %d' %
                    (sentence[position], position), sentence, position)
            var(variable)
            return Symbol(variable)
        return self._tensor()

    # <NUMBER> -> <RATIONAL> | <DECIMAL> | <INTEGER>
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
        sentence = self.lexer.sentence
        position = self.lexer.index - len(self.lexer.lexeme)
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    # <COMMAND> -> <SQRT> | <FRAC>
    def _command(self):
        command = self.lexer.lexeme
        if self.peek('SQRT_CMD'):
            return self._sqrt()
        if self.peek('FRAC_CMD'):
            return self._frac()
        if self.peek('DERV_CMD'):
            return self._derivative()
        position = self.lexer.index - len(self.lexer.lexeme)
        raise ParseError('unsupported command \'%s\' at position %d' %
            (command, position), self.lexer.sentence, position)

    # <SQRT> -> '\\sqrt' [ '[' <INTEGER> ']' ] '{' <EXPRESSION> '}'
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
        return Pow(expr, root)

    # <FRAC> -> '\\frac' '{' <EXPRESSION> '}' '{' <EXPRESSION> '}'
    def _frac(self):
        self.expect('FRAC_CMD')
        self.expect('LEFT_BRACE')
        numerator = self._expression()
        self.expect('RIGHT_BRACE')
        self.expect('LEFT_BRACE')
        denominator = self._expression()
        self.expect('RIGHT_BRACE')
        return numerator / denominator

    # <DERIVATIVE> -> { '\\partial' '_' <SYMBOL> }+ ( <TENSOR> | '(' <EXPRESSION> ')' )
    def _derivative(self):
        indices = []
        while True:
            self.expect('DERV_CMD')
            self.expect('UNDERSCORE')
            index = self.lexer.lexeme
            self.expect('SYMBOL')
            indices.extend(index)
            if not self.peek('DERV_CMD'): break
        if self.accept('LEFT_PAREN'):
            tree = ExprTree(self._expression())
            self.expect('RIGHT_PAREN')
            for subtree in tree.preorder():
                subexpr = subtree.expr
                if subexpr.func == Function('Tensor'):
                    subtree.expr = Function('Tensor')(Symbol('_x'), *subexpr.args)
                    tree.build(subtree)
            expr = tree.reconstruct()
            tree = ExprTree(diff(expr, Symbol('_x')))
            for subtree in tree.preorder():
                subexpr = subtree.expr
                if subexpr.func == Derivative:
                    order, array = len(indices), subexpr.args[1]
                    indexing = list(array.args[1:-1]) + indices
                    tensor = self.namespace[str(array.args[0])]
                    subtree.expr = self._differentiate(tensor, order, *indexing)
                    tree.build(subtree)
                elif subexpr.func == Function('Tensor'):
                    subtree.expr = Function('Tensor')(*subexpr.args[1:])
                    tree.build(subtree)
            return tree.reconstruct()
        order, array = len(indices), self._tensor()
        indices = list(array.args[1:]) + indices
        tensor = self.namespace[str(array.args[0])]
        return self._differentiate(tensor, order, *indices)

    # <CONFIG> -> '%' <TENSOR> '[' <INTEGER> ']' [ ':' <SYMMETRY> ] { ',' <TENSOR> '[' <INTEGER> ']' [ ':' <SYMMETRY> ] }*
    def _config(self):
        self.expect('COMMENT')
        while True:
            array = self._tensor()
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
            tensor = Tensor(array, dimension,
                'sym01' if symmetry == 'metric' else 'nosym' if symmetry == 'permutation' else symmetry,
                invertible=(symmetry == 'metric'), permutation=(symmetry == 'permutation'))
            self.namespace[tensor.name] = tensor
            if tensor.inverse and tensor.rank == 2:
                inverse = tensor.name.replace('U', 'D') if 'U' in tensor.name else tensor.name.replace('D', 'U')
                self.namespace[inverse] = self.namespace[tensor.name].inverse
            if not self.accept('COMMA'): break
        self.lexer.update(self.namespace)

    # <TENSOR> -> <ARRAY> ( '_' <LOWER_INDEX> ) | ( '^' <UPPER_INDEX> [ '_' <LOWER_INDEX> ] ] )
    def _tensor(self):
        if not (self.peek('SYMBOL') or self.peek('ARRAY')):
            sentence = self.lexer.sentence
            position = self.lexer.index - len(self.lexer.lexeme)
            raise ParseError('unexpected \'%s\' at position %d' %
                (sentence[position], position), sentence, position)
        array = self.lexer.lexeme
        self.lexer.lex()
        if array[0] == '\\':
            array = array[1:]
        array, indices = list(array), []
        if self.accept('UNDERSCORE'):
            index, order = self._lower_index()
            indices.extend(index)
            array.extend((len(index) - order) * ['D'])
            if order > 0:
                tensor = self.namespace[''.join(array)]
                return self._differentiate(tensor, order, *indices)
        if self.accept('CARET'):
            index = self._upper_index()
            indices.extend(index)
            array.extend(len(index) * ['U'])
            if self.accept('UNDERSCORE'):
                index, order = self._lower_index()
                indices.extend(index)
                array.extend((len(index) - order) * ['D'])
                if order > 0:
                    tensor = self.namespace[''.join(array)]
                    return self._differentiate(tensor, order, *indices)
        return Function('Tensor')(Symbol(''.join(array)), *indices)

    # <LOWER_INDEX> -> <SYMBOL> | '{' { <SYMBOL> }* [ ',' { <SYMBOL> }+ ] '}'
    def _lower_index(self):
        index, order = [], 0
        def append_index():
            symbol = self.lexer.lexeme
            self.lexer.lex()
            if symbol[0] == '\\':
                symbol = symbol[1:]
            index.append(Symbol(symbol))
        if self.peek('SYMBOL'):
            append_index()
        if self.accept('LEFT_BRACE'):
            while self.peek('SYMBOL'):
                append_index()
            if self.accept('COMMA'):
                while self.peek('SYMBOL'):
                    order += 1
                    append_index()
            self.expect('RIGHT_BRACE')
        return index, order

    # <UPPER_INDEX> -> <SYMBOL> | '{' { <SYMBOL> }+ '}'
    def _upper_index(self):
        index = []
        def append_index():
            symbol = self.lexer.lexeme
            self.lexer.lex()
            if symbol[0] == '\\':
                symbol = symbol[1:]
            index.append(Symbol(symbol))
        if self.peek('SYMBOL'):
            append_index()
        if self.accept('LEFT_BRACE'):
            while self.peek('SYMBOL'):
                append_index()
            self.expect('RIGHT_BRACE')
        return index

    def _tensorial(self):
        mark_1 = self.lexer.index - 1
        for token in self.lexer.tokenize():
            if token == 'LINE_BREAK': break
        mark_2 = self.lexer.index
        assignment = self.lexer.sentence[mark_1:mark_2]
        self.lexer.reset(mark_1); self.lexer.lex()
        nameset = set(re.match(r'[^UD]*', tensor).group() for tensor in self.namespace)
        return any(tensor in assignment for tensor in nameset)

    def _differentiate(self, tensor, order, *index):
        name = tensor.name + ('' if '_d' in tensor.name else '_d') + order * 'D'
        rank, dimension, symmetry = tensor.rank, tensor.dimension, tensor.symmetry
        function = Function('Tensor')(Symbol(name), *index)
        if order == 2:
            if symmetry and symmetry != 'nosym':
                symmetry = tensor.symmetry + '_sym%d%d' % (rank, rank + order - 1)
            else: symmetry = 'sym%d%d' % (rank, rank + order - 1)
        self.namespace[name] = Tensor(function, dimension, symmetry)
        return function

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

class Tensor:
    """ Tensor Structure """

    def __init__(self, function, dimension, symmetry=None, permutation=False, invertible=False):
        self.function = function
        self.name = str(function.args[0])
        self.rank = len(function.args) - 1
        self.dimension = dimension
        self.symmetry  = symmetry
        if '_d' in self.name and self.name in Parser.namespace:
            self.array = Parser.namespace[self.name]
        elif permutation:
            index = [chr(105 + n) for n in range(self.rank)]
            prefix = '[' * self.rank + 'sgn([' + ', '.join(index) + '])'
            suffix = ''.join(' for %s in range(%d)]' % (index[self.rank - i], dimension) for i in range(1, self.rank + 1))
            self.array = eval(prefix + suffix, {'sgn': self._sgn})
        else: self.array = ixp.declare_indexedexp(self.rank, self.name, symmetry, dimension)
        if invertible and symmetry == 'sym01':
            if   dimension == 2:
                self.inverse = ixp.symm_matrix_inverter2x2(self.array)[0]
            elif dimension == 3:
                self.inverse = ixp.symm_matrix_inverter3x3(self.array)[0]
            elif dimension == 4:
                self.inverse = ixp.symm_matrix_inverter4x4(self.array)[0]
        else: self.inverse = None

    @staticmethod
    def _sgn(sequence):
        cycle_length = 0
        for n, i in enumerate(sequence[:-1]):
            for j in sequence[(n + 1):]:
                if i == j: return 0
                cycle_length += i > j
        return (-1)**cycle_length

    @staticmethod
    def notation(function):
        fields = [str(arg) for arg in function.args]
        suffix = ''.join(['[' + n + ']' for n in fields[1:]])
        return fields[0] + suffix

    def __repr__(self):
        return '%s -> (dimension: %d, symmetry: %s)' % \
            (self.function, self.dimension, self.symmetry)

    __str__ = __repr__

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

def parse(sentence, expression=False, debug=False):
    """ Convert LaTeX Sentence to SymPy Expression

        :arg:    latex sentence (raw string)
        :arg:    expression mode [default: False]
        :arg:    debug mode [default: None]
        :return: symbolic expression or tensor namespace

        >>> parse(r'-(\\frac{2}{3} + 2\\sqrt[5]{x + 3})', expression=True)
        -2*(x + 3)**(1/5) - 2/3
        >>> parse(r's_n = (1 + 1/n)^n')
        ['s_n']
        >>> print(s_n)
        (1 + 1/n)**n

        >>> parse(r'% h^\\mu_\\mu [3]: nosym; h = h^\\mu_\\mu')
        ['hUD', 'h']
        >>> print(hUD)
        [[hUD00, hUD01, hUD02], [hUD10, hUD11, hUD12], [hUD20, hUD21, hUD22]]
        >>> print(h)
        hUD00 + hUD11 + hUD22

        >>> config = r'% g^{ij} [2]: metric, v_j [2], w_j [2];'
        >>> parse(config + r'u^i = g^{ij}(v_j + w_j)')
        ['gUU', 'gDD', 'vD', 'wD', 'uU']
        >>> print(gDD)
        [[gUU11/(gUU00*gUU11 - gUU01**2), -gUU01/(gUU00*gUU11 - gUU01**2)], [-gUU01/(gUU00*gUU11 - gUU01**2), gUU00/(gUU00*gUU11 - gUU01**2)]]
        >>> print(uU)
        [gUU00*vD0 + gUU00*wD0 + gUU01*vD1 + gUU01*wD1, gUU01*vD0 + gUU01*wD0 + gUU11*vD1 + gUU11*wD1]

        >>> config = r'% \\epsilon^{ijk} [3]: permutation, b_k [3];'
        >>> parse(config + r'a^i = \\epsilon^{ijk} b_{k,j}')
        ['epsilonUUU', 'bD', 'bD_dD', 'aU']
        >>> print(epsilonUUU)
        [[[0, 0, 0], [0, 0, 1], [0, -1, 0]], [[0, 0, -1], [0, 0, 0], [1, 0, 0]], [[0, 1, 0], [-1, 0, 0], [0, 0, 0]]]
        >>> print(aU)
        [-bD_dD12 + bD_dD21, bD_dD02 - bD_dD20, -bD_dD01 + bD_dD10]
    """
    namespace = Parser().parse(sentence, expression)
    if not isinstance(namespace, dict):
        return namespace
    if not debug:
        # update global scope with tensor namespace
        globals().update(namespace)
        for var in namespace:
            if isinstance(namespace[var], str):
                # evaluate each implied summation
                namespace[var] = eval(namespace[var])
                globals()[var] = namespace[var]
        Parser.namespace.update(namespace)
    # inject namespace into the previous stack frame
    frame = currentframe().f_back
    frame.f_globals.update(namespace)
    return list(namespace.keys())

if __name__ == "__main__":
    import doctest
    doctest.testmod()
