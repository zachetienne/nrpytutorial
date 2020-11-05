""" Convert LaTeX Sentence to SymPy Expression """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

from sympy import Function, Derivative, Symbol, Integer, Rational, Float, Pow, Add
from sympy import sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
from sympy import pi, exp, log, sqrt, expand, diff
from inspect import currentframe
from functional import uniquify
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
            [ ('WHITESPACE',    r'\s+'),
              ('DIMENSION',     r'[2-9][0-9]*D'),
              ('VARIABLE',      r'[a-zA-Z]+[UD]+'),
              ('STRING',        r'\"[^\"]+\"'),
              ('RATIONAL',      r'\-?[0-9]+\/\-?[1-9][0-9]*'),
              ('DECIMAL',       r'\-?[0-9]+\.[0-9]+'),
              ('INTEGER',       r'\-?[0-9]+'),
              ('PI',            r'\\pi'),
              ('COMMENT',       r'\%'),
              ('ARROW',         r'\-\>'),
              ('PLUS',          r'\+'),
              ('MINUS',         r'\-'),
              ('DIVIDE',        r'\/'),
              ('EQUAL',         r'\='),
              ('CARET',         r'\^'),
              ('UNDERSCORE',    r'\_'),
              ('COMMA',         r'\,'),
              ('COLON',         r'\:'),
              ('APOSTROPHE',    r'\''),
              ('LPAREN',        r'\('),
              ('RPAREN',        r'\)'),
              ('LBRACK',        r'\['),
              ('RBRACK',        r'\]'),
              ('LBRACE',        r'\{'),
              ('RBRACE',        r'\}'),
              ('LINE_BREAK',    r'\;|\\\\'),
              ('OPENING',       r'\\begin{align\*?}'),
              ('CLOSING',       r'\\end{align\*?}'),
              ('PAR_SYM',       r'\\partial'),
              ('COV_SYM',       r'\\nabla|D'),
              ('LIE_SYM',       r'\\mathcal\{L\}'),
              ('FUNC_CMD',      r'\\exp'),
              ('FRAC_CMD',      r'\\frac'),
              ('SQRT_CMD',      r'\\sqrt'),
              ('NLOG_CMD',      r'\\ln|\\log'),
              ('TRIG_CMD',      r'\\sinh|\\cosh|\\tanh|\\sin|\\cos|\\tan'),
              ('ASSIGN_MACRO',  r'assign'),
              ('DEFINE_MACRO',  r'define'),
              ('IGNORE_MACRO',  r'ignore'),
              ('PARSE_MACRO',   r'parse'),
              ('ALIAS_MACRO',   r'alias'),
              ('INDEX_KWRD',    r'index'),
              ('BASIS_KWRD',    r'basis'),
              ('DERIV_KWRD',    r'deriv'),
              ('DIACRITIC',     r'\\hat|\\tilde|\\bar'),
              ('VPHANTOM',      r'\\vphantom'),
              ('MATHIT',        r'\\mathit'),
              ('WEIGHT',        r'weight'),
              ('MODIFIER',      r'symbolic|variable'),
              ('SYMMETRY',      r'const|metric|' + symmetry),
              ('LETTER',        r'[a-zA-Z]|' + greek_pattern),
              ('COMMAND',       r'\\[a-zA-Z]+')]]))

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
            if token.lastgroup != 'WHITESPACE':
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
        <LATEX>         -> ( <ALIGN> | <CONFIG> | <ASSIGNMENT> ) { <LINE_BREAK> ( <ALIGN> | <CONFIG> | <ASSIGNMENT> ) }*
        <ALIGN>         -> <OPENING> ( <CONFIG> | <ASSIGNMENT> ) { <LINE_BREAK> ( <CONFIG> | <ASSIGNMENT> ) }* <CLOSING>
        <CONFIG>        -> '%' ( <PARSE> | <ALIAS> | <ASSIGN> | <DEFINE> )
        <PARSE>         -> <PARSE_MACRO> <ASSIGNMENT> { ',' <ASSIGNMENT> }*
        <ALIAS>         -> <ALIAS_MACRO> <STRING> <ARROW> <STRING> { ',' <STRING> <ARROW> <STRING> }*
        <ASSIGN>        -> <ASSIGN_MACRO> ( <SYMMETRY> | <WEIGHT> <NUMBER> ) ( <LETTER> | <VARIABLE> )
        <DEFINE>        -> <DEFINE_MACRO> ( <VARDEF> | <KEYDEF> ) { ',' ( <VARDEF> | <KEYDEF> ) }*
        <IGNORE>        -> <IGNORE_MACRO> <STRING> { ',' <STRING> }*
        <VARDEF>        -> [ <SYMMETRY> ] ( <LETTER> | <VARIABLE> ) [ '(' <DIMENSION> ')' ]
        <KEYDEF>        -> <BASIS_KWRD> <BASIS> | <DERIV_KWRD> <MODIFIER> | <INDEX_KWRD> <RANGE>
        <BASIS>         -> <BASIS_KWRD> '{' <LETTER> { ',' <LETTER> }* '}'
        <RANGE>         -> ( <LETTER> | '[' <LETTER> '-' <LETTER> ']' ) '=' <INTEGER> ':' <INTEGER>
        <ASSIGNMENT>    -> ( <TENSOR> | <OPERATOR> ) = <EXPRESSION>
        <EXPRESSION>    -> <TERM> { ( '+' | '-' ) <TERM> }*
        <TERM>          -> <FACTOR> { [ '/' ] <FACTOR> }*
        <FACTOR>        -> <BASE> { '^' <EXPONENT> }*
        <BASE>          -> [ '-' ] ( <ATOM> | <SUBEXPR> )
        <EXPONENT>      -> <BASE> | '{' <BASE> '}' | '{{' <BASE> '}}'
        <ATOM>          -> <COMMAND> | <OPERATOR> | <NUMBER> | <TENSOR>
        <SUBEXPR>       -> '(' <EXPRESSION> ')' | '[' <EXPRESSION> ']' | '{' <EXPRESSION> '}'
        <COMMAND>       -> <FUNC> | <FRAC> | <SQRT> | <NLOG> | <TRIG>
        <FUNC>          -> <FUNC_CMD> '(' <EXPRESSION> ')'
        <FRAC>          -> <FRAC_CMD> '{' <EXPRESSION> '}' '{' <EXPRESSION> '}'
        <SQRT>          -> <SQRT_CMD> [ '[' <INTEGER> ']' ] '{' <EXPRESSION> '}'
        <NLOG>          -> <NLOG_CMD> [ '_' ( <NUMBER> | '{' <NUMBER> '}' ) ] ( <NUMBER> | <TENSOR> | '(' <EXPRESSION> ')' )
        <TRIG>          -> <TRIG_CMD> [ '^' ( <NUMBER> | '{' <NUMBER> '}' ) ] ( <NUMBER> | <TENSOR> | '(' <EXPRESSION> ')' )
        <OPERATOR>      -> [ <VPHANTOM> '{' <MODIFIER> '}' ] ( <PARDRV> | <COVDRV> | <LIEDRV> )
        <PARDRV>        -> { <PAR_SYM> [ '^' <INTEGER> ] '_' <LETTER> }+ ( <TENSOR> | <SUBEXPR> )
        TODO: <COVDRV>  -> { ( <COV_SYM> | <DIACRITIC> '{' <COV_SYM> '}' ) ( '^' | '_' ) <LETTER> }+ ( <TENSOR> | <SUBEXPR> )
        TODO: <LIEDRV>  -> <LIE_SYM> '_' <SYMBOL> ( <TENSOR> | <SUBEXPR> )
        <NUMBER>        -> <RATIONAL> | <DECIMAL> | <INTEGER> | <PI>
        <TENSOR>        -> <SYMBOL> [ ( '_' <LOWER_INDEX> ) | ( '^' <UPPER_INDEX> [ '_' <LOWER_INDEX> ] ) ]
        <SYMBOL>        -> <LETTER> | <DIACRITIC> '{' <LETTER> '}' | <MATHIT> '{' <LETTER> { '_' | <LETTER> | <INTEGER> }* '}'
        <LOWER_INDEX>   -> <LETTER> | <INTEGER> | '{' { <LETTER> | <INTEGER> }* [ ( ',' | ';' ) { <LETTER> }+ ] '}'
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
        if 'ignore' not in self._namespace:
            self._namespace['ignore'] = ['\\left', '\\right', '\\,', '{}', '&']
        if 'metric' not in self._namespace:
            self._namespace['metric'] = {'': 'g', 'bar': 'g', 'hat': 'g', 'tilde': 'gamma'}
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
        for ignore in self._namespace['ignore']:
            sentence = sentence.replace(ignore, '')
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
        self._latex()
        return self._namespace

    # <LATEX> -> ( <ALIGN> | <CONFIG> | <ASSIGNMENT> ) { <LINE_BREAK> ( <ALIGN> | <CONFIG> | <ASSIGNMENT> ) }*
    def _latex(self):
        if self.peek('OPENING'):
            self._align()
        elif self.peek('COMMENT'):
            self._config()
        else: self._assignment()
        while self.accept('LINE_BREAK') or self.peek('COMMENT'):
            if self.peek('OPENING'):
                self._align()
            elif self.peek('COMMENT'):
                self._config()
            else: self._assignment()

    # <ALIGN> -> <OPENING> ( <CONFIG> | <ASSIGNMENT> ) { <LINE_BREAK> ( <CONFIG> | <ASSIGNMENT> ) }* <CLOSING>
    def _align(self):
        self.expect('OPENING')
        if self.peek('COMMENT'):
            self._config()
        else: self._assignment()
        while self.accept('LINE_BREAK') or self.peek('COMMENT'):
            if self.peek('COMMENT'):
                self._config()
            else: self._assignment()
        self.expect('CLOSING')

    # <CONFIG> -> '%' ( <PARSE> | <ASSIGN> | <DEFINE> )
    def _config(self):
        self.expect('COMMENT')
        if self.peek('PARSE_MACRO'):
            self._parse()
        elif self.peek('ALIAS_MACRO'):
            self._alias()
        elif self.peek('ASSIGN_MACRO'):
            self._assign()
        elif self.peek('DEFINE_MACRO'):
            self._define()
        elif self.peek('IGNORE_MACRO'):
            self._ignore()
        else:
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('unsupported macro at position %d' %
                position, sentence, position)

    # <PARSE> -> <PARSE_MACRO> <ASSIGNMENT> { ',' <ASSIGNMENT> }*
    def _parse(self):
        self.expect('PARSE_MACRO')
        self._assignment()
        while self.accept('COMMA'):
            self._assignment()

    # <ALIAS> -> <ALIAS_MACRO> <STRING> <ARROW> <STRING> { ',' <STRING> <ARROW> <STRING> }*
    def _alias(self):
        self.expect('ALIAS_MACRO')
        while True:
            old = self.lexer.lexeme[1:-1]
            self.expect('STRING')
            self.expect('ARROW')
            new = self.lexer.lexeme[1:-1]
            self.expect('STRING')
            sentence, position = self.lexer.sentence, self.lexer.mark()
            self.lexer.sentence = sentence[:position] + sentence[position:].replace(old, new)
            if not self.accept('COMMA'): break

    # <ASSIGN> -> <ASSIGN_MACRO> ( <SYMMETRY> | <WEIGHT> <NUMBER> ) ( <LETTER> | <VARIABLE> )
    def _assign(self):
        self.expect('ASSIGN_MACRO')
        symmetry, weight = None, None
        if self.peek('SYMMETRY'):
            symmetry = self.lexer.lexeme
            self.lexer.lex()
        elif self.accept('WEIGHT'):
            weight = self._number()
        else:
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('unexpected \'%s\' at position %d' %
                (sentence[position], position), sentence, position)
        if self.peek('LETTER'):
            symbol = self.lexer.lexeme
            self.expect('LETTER')
        elif self.peek('VARIABLE'):
            symbol = self.lexer.lexeme
            self.expect('VARIABLE')
        else:
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('unexpected \'%s\' at position %d' %
                (sentence[position], position), sentence, position)
        if symbol not in self._namespace:
            raise TensorError('cannot update undefined tensor \'%s\'' % symbol)
        tensor = self._namespace[symbol]
        if weight: tensor.weight = weight
        structure, dimension = tensor.structure, tensor.dimension
        if symmetry == 'metric':
            if symmetry == 'metric':
                diacritical = False
                for diacritic in ('bar', 'hat', 'tilde'):
                    if diacritic in symbol:
                        self._namespace['metric'][diacritic] = symbol.split(diacritic)[0]
                        diacritical = True
                if not diacritical:
                    self._namespace['metric'][''] = re.split(r'[UD]', symbol)[0]
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

    # <DEFINE> -> <DEFINE_MACRO> ( <VARDEF> | <KEYDEF> ) { ',' ( <VARDEF> | <KEYDEF> ) }*
    def _define(self):
        self.expect('DEFINE_MACRO')
        if any(self.peek(kwrd) for kwrd in ('BASIS_KWRD', 'DERIV_KWRD', 'INDEX_KWRD')):
            self._keydef()
        else: self._vardef()
        while self.accept('COMMA'):
            if any(self.peek(kwrd) for kwrd in ('BASIS_KWRD', 'DERIV_KWRD', 'INDEX_KWRD')):
                self._keydef()
            else: self._vardef()

    # <IGNORE> -> <IGNORE_MACRO> <STRING> { ',' <STRING> }*
    def _ignore(self):
        self.expect('IGNORE_MACRO')
        while True:
            string = self.lexer.lexeme[1:-1]
            self.expect('STRING')
            sentence, position = self.lexer.sentence, self.lexer.mark()
            self.lexer.sentence = sentence[:position] + sentence[position:].replace(string, '')
            if not self.accept('COMMA'): break

    # <VARDEF> -> [ <SYMMETRY> ] ( <LETTER> | <VARIABLE> ) [ '(' <DIMENSION> ')' ]
    def _vardef(self):
        symmetry = self.lexer.lexeme
        if self.peek('SYMMETRY'):
            self.lexer.lex()
        else: symmetry = None
        if self.peek('LETTER'):
            symbol = self.lexer.lexeme
            self.expect('LETTER')
        elif self.peek('VARIABLE'):
            symbol = self.lexer.lexeme
            self.expect('VARIABLE')
        else:
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('unexpected \'%s\' at position %d' %
                (sentence[position], position), sentence, position)
        if self.accept('LPAREN'):
            dimension = self.lexer.lexeme[:-1]
            self.expect('DIMENSION')
            dimension = int(dimension)
            self.expect('RPAREN')
        else: dimension = None
        if symmetry != 'const' and not dimension:
            raise TensorError('dimension only omittable for constant')
        if symmetry == 'const':
            self._namespace[symbol] = Function('Constant')(Symbol(symbol))
        else:
            if symmetry == 'metric':
                diacritical = False
                for diacritic in ('bar', 'hat', 'tilde'):
                    if diacritic in symbol:
                        self._namespace['metric'][diacritic] = symbol.split(diacritic)[0]
                        diacritical = True
                if not diacritical:
                    self._namespace['metric'][''] = re.split(r'[UD]', symbol)[0]
            tensor = Tensor(Function('Tensor')(symbol), dimension)
            self._define_tensor(tensor, symmetry, invertible=(symmetry == 'metric'),
                permutation=(symbol[:7] == 'epsilon'), kronecker=(symbol[:5] == 'delta'))

    # <KEYDEF> -> <BASIS_KWRD> <BASIS> | <DERIV_KWRD> <MODIFIER> | <INDEX_KWRD> <RANGE>
    def _keydef(self):
        if self.accept('BASIS_KWRD'):
            self._basis()
        elif self.accept('DERIV_KWRD'):
            modifier = self.lexer.lexeme
            self.expect('MODIFIER')
            self._namespace['deriv'] = modifier
        elif self.accept('INDEX_KWRD'):
            self._range()
        else:
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('unexpected keyword at position %d' %
                position, sentence, position)

    # <BASIS> -> <BASIS_KWRD> <LBRACK> <LETTER> [ ',' <LETTER> ]* <RBRACK>
    def _basis(self):
        self.expect('LBRACK')
        del self._namespace['basis'][:]
        while True:
            symbol = self._strip(self.lexer.lexeme)
            self.expect('LETTER')
            if symbol in self._namespace['basis']:
                sentence, position = self.lexer.sentence, self.lexer.mark()
                raise ParseError('duplicate basis symbol \'%s\' at position %d' %
                    (sentence[position], position), sentence, position)
            self._namespace['basis'].append(Symbol(symbol))
            if not self.accept('COMMA'): break
        self.expect('RBRACK')

    # <RANGE> -> ( <LETTER> | '[' <LETTER> '-' <LETTER> ']' ) '=' <INTEGER> ':' <INTEGER>
    def _range(self):
        if self.accept('LBRACK'):
            index_1 = self._strip(self.lexer.lexeme)
            self.expect('LETTER')
            self.expect('MINUS')
            index_2 = self._strip(self.lexer.lexeme)
            self.expect('LETTER')
            index = [chr(i) for i in range(ord(index_1), ord(index_2) + 1)]
            self.expect('RBRACK')
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

    # <ASSIGNMENT> -> ( <TENSOR> | <OPERATOR> ) = <EXPRESSION>
    def _assignment(self):
        pardrv = self.peek('PAR_SYM')
        covdrv = self.peek('COV_SYM')
        liedrv = self.peek('LIE_SYM')
        self.lexer.mark()
        if self.accept('DIACRITIC'):
            self.expect('LBRACE')
            covdrv = self.peek('COV_SYM')
            self.lexer.reset()
        if pardrv and self._namespace['deriv'] == 'symbolic':
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('cannot parse symbolic partial derivative on LHS' %
                position, sentence, position)
        LHS = self._pardrv() if pardrv \
         else self._covdrv('LHS') if covdrv \
         else self._liedrv('LHS') if liedrv \
         else self._tensor()
        indexed = LHS.func == Function('Tensor') and len(LHS.args) > 1
        self.expect('EQUAL')
        tree = ExprTree(self._expression())
        if not indexed:
            for subtree in tree.preorder():
                subexpr, rank = subtree.expr, len(subtree.expr.args)
                if subexpr.func == Function('Tensor') and rank > 1:
                    indexed = True
        function, RHS = LHS, expand(tree.root.expr)
        if indexed:
            # perform implied summation on indexed expression
            (LHS, RHS), dimension = self._summation(LHS, RHS)
        else:
            for subtree in tree.preorder():
                subexpr = subtree.expr
                if subexpr.func in (Function('Tensor'), Function('Constant')):
                    subtree.expr = subexpr.args[0]
                    del subtree.children[:]
            (LHS, RHS), dimension = (LHS.args[0], tree.reconstruct()), 0
        global_env = dict(sympy_env)
        global_env.update(self._namespace)
        for key in global_env:
            if isinstance(global_env[key], Tensor):
                global_env[key] = global_env[key].structure
            if isinstance(global_env[key], Function('Constant')):
                global_env[key] = global_env[key].args[0]
        # evaluate every implied summation and update namespace
        exec('%s = %s' % (LHS, RHS), global_env)
        symbol = LHS.split('[')[0] if indexed else str(LHS)
        tensor = Tensor(function, dimension, global_env[symbol])
        self._namespace.update({symbol: tensor})

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
        while any(self.peek(token) for token in ('DIVIDE', 'LPAREN', 'VPHANTOM', 'DIACRITIC',
                'RATIONAL', 'DECIMAL', 'INTEGER', 'PI', 'PAR_SYM', 'COV_SYM', 'LIE_SYM',
                'LETTER', 'COMMAND', 'FUNC_CMD', 'FRAC_CMD', 'SQRT_CMD', 'NLOG_CMD', 'TRIG_CMD')):
            if self.accept('DIVIDE'):
                expr /= self._factor()
            else: expr *= self._factor()
        return expr

    # <FACTOR> -> <BASE> { '^' <EXPONENT> }*
    def _factor(self):
        stack = [self._base()]
        while self.accept('CARET'):
            stack.append(self._exponent())
        if len(stack) == 1: stack.append(1)
        expr = stack.pop()
        for subexpr in reversed(stack):
            expr = exp(expr) if subexpr == Function('Tensor')('e') \
              else subexpr ** expr
        return expr

    # <BASE> -> [ '-' ] ( <ATOM> | <SUBEXPR> )
    def _base(self):
        sign = -1 if self.accept('MINUS') else 1
        if any(self.peek(i) for i in ('LPAREN', 'LBRACK', 'LBRACE')):
            expr = sign * self._subexpr()
            return expr
        return sign * self._atom()

    # <EXPONENT> -> <BASE> | '{' <BASE> '}' | '{{' <BASE> '}}'
    def _exponent(self):
        if self.accept('LBRACE'):
            if self.accept('LBRACE'):
                base = self._base()
                self.expect('RBRACE')
            else: base = self._base()
            self.expect('RBRACE')
            return base
        return self._base()

    # <ATOM> -> <COMMAND> | <OPERATOR> | <NUMBER> | <TENSOR>
    def _atom(self):
        self.lexer.mark()
        if self.accept('DIACRITIC'):
            self.expect('LBRACE')
            if self.peek('COV_SYM'):
                self.lexer.reset()
                return self._operator()
            self.lexer.reset()
        if any(self.peek(token) for token in
                ('COMMAND', 'FUNC_CMD', 'FRAC_CMD', 'SQRT_CMD', 'NLOG_CMD', 'TRIG_CMD')):
            return self._command()
        if any(self.peek(token) for token in
                ('VPHANTOM', 'PAR_SYM', 'COV_SYM', 'LIE_SYM')):
            return self._operator()
        if any(self.peek(token) for token in
                ('RATIONAL', 'DECIMAL', 'INTEGER', 'PI')):
            return self._number()
        if any(self.peek(token) for token in
                ('LETTER', 'DIACRITIC', 'MATHIT')):
            sentence, position = self.lexer.sentence, self.lexer.mark()
            tensor = Tensor(self._tensor(), None)
            symbol, indexing = tensor.symbol, tensor.indexing
            # reserved keyword for christoffel symbol
            if symbol[:5] == 'Gamma' and tensor.rank == 3:
                metric = self._namespace['metric'][symbol[5:-3]] + symbol[5:-3]
                if metric + 'DD' not in self._namespace:
                    raise ParseError('cannot generate christoffel symbol without defined metric \'%s\'' %
                        metric, sentence, position)
                sentence, position = self.lexer.sentence, self.lexer.mark()
                self.parse(self._generate_christoffel(tensor.function, self._namespace['metric']))
                self.lexer.initialize(sentence, position)
                self.lexer.lex()
            if symbol not in self._namespace:
                if indexing:
                    raise ParseError('cannot index undefined tensor \'%s\' at position %d' %
                        (symbol, position), sentence, position)
                self._define_tensor(Tensor(tensor.function, 0))
            return tensor.function
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    # <SUBEXPR> -> '(' <EXPRESSION> ')' | '[' <EXPRESSION> ']' | '{' <EXPRESSION> '}'
    def _subexpr(self):
        if self.accept('LPAREN'):
            expr = self._expression()
            self.expect('RPAREN')
        elif self.accept('LBRACK'):
            expr = self._expression()
            self.expect('RBRACK')
        elif self.accept('LBRACE'):
            expr = self._expression()
            self.expect('RBRACE')
        else:
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('unexpected \'%s\' at position %d' %
                (sentence[position], position), sentence, position)
        return expr

    # <COMMAND> -> <FUNC> | <FRAC> | <SQRT> | <NLOG> | <TRIG>
    def _command(self):
        command = self.lexer.lexeme
        if self.peek('FUNC_CMD'):
            return self._func()
        if self.peek('FRAC_CMD'):
            return self._frac()
        if self.peek('SQRT_CMD'):
            return self._sqrt()
        if self.peek('NLOG_CMD'):
            return self._nlog()
        if self.peek('TRIG_CMD'):
            return self._trig()
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unsupported command \'%s\' at position %d' %
            (command, position), sentence, position)

    # <FRAC> -> <FRAC_CMD> '{' <EXPRESSION> '}' '{' <EXPRESSION> '}'
    def _frac(self):
        self.expect('FRAC_CMD')
        self.expect('LBRACE')
        numerator = self._expression()
        self.expect('RBRACE')
        self.expect('LBRACE')
        denominator = self._expression()
        self.expect('RBRACE')
        return numerator / denominator

    # <SQRT> -> <SQRT_CMD> [ '[' <INTEGER> ']' ] '{' <EXPRESSION> '}'
    def _sqrt(self):
        self.expect('SQRT_CMD')
        if self.accept('LBRACK'):
            integer = self.lexer.lexeme
            self.expect('INTEGER')
            root = Rational(1, integer)
            self.expect('RBRACK')
        else: root = Rational(1, 2)
        self.expect('LBRACE')
        expr = self._expression()
        self.expect('RBRACE')
        if root == Rational(1, 2):
            return sqrt(expr)
        return Pow(expr, root)

    # <FUNC> -> <FUNC_CMD> '(' <EXPRESSION> ')'
    def _func(self):
        func = self._strip(self.lexer.lexeme)
        self.expect('FUNC_CMD')
        self.expect('LPAREN')
        expr = self._expression()
        self.expect('RPAREN')
        if func == 'exp':
            return exp(expr)
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unsupported function \'%s\' at position %d' %
            (func, position), sentence, position)

    # <NLOG> -> <NLOG_CMD> [ '_' ( <NUMBER> | '{' <NUMBER> '}' ) ] ( <NUMBER> | <TENSOR> | '(' <EXPRESSION> ')' )
    def _nlog(self):
        func = self._strip(self.lexer.lexeme)
        self.expect('NLOG_CMD')
        if func == 'log':
            if self.accept('UNDERSCORE'):
                if self.accept('LBRACE'):
                    base = self._number()
                    self.expect('RBRACE')
                else:
                    base = self._number()
                base = int(base)
            else: base = 10
        if any(self.peek(token) for token in
                ('RATIONAL', 'DECIMAL', 'INTEGER', 'PI')):
            expr = self._number()
        elif any(self.peek(token) for token in
                ('LETTER', 'DIACRITIC', 'MATHIT')):
            sentence, position = self.lexer.sentence, self.lexer.mark()
            tensor = Tensor(self._tensor(), None)
            symbol, indexing = tensor.symbol, tensor.indexing
            if symbol not in self._namespace:
                if indexing:
                    raise ParseError('cannot index undefined tensor \'%s\' at position %d' %
                        (symbol, position), sentence, position)
                self._define_tensor(Tensor(tensor.function, 0))
            expr = tensor.function
        elif self.accept('LPAREN'):
            expr = self._expression()
            self.expect('RPAREN')
        else:
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('unexpected \'%s\' at position %d' %
                (sentence[position], position), sentence, position)
        if func == 'ln': return log(expr)
        return log(expr, base)

    # <TRIG> -> <TRIG_CMD> [ '^' ( <NUMBER> | '{' <NUMBER> '}' ) ] ( <NUMBER> | <TENSOR> | '(' <EXPRESSION> ')' )
    def _trig(self):
        func = self._strip(self.lexer.lexeme)
        self.expect('TRIG_CMD')
        if self.accept('CARET'):
            if self.accept('LBRACE'):
                exponent = self._number()
                self.expect('RBRACE')
            else:
                exponent = self._number()
            exponent = int(exponent)
        else: exponent = 1
        if   func == 'cosh': trig = acosh if exponent == -1 else cosh
        elif func == 'sinh': trig = asinh if exponent == -1 else sinh
        elif func == 'tanh': trig = atanh if exponent == -1 else tanh
        elif func == 'cos':  trig = acos  if exponent == -1 else cos
        elif func == 'sin':  trig = asin  if exponent == -1 else sin
        elif func == 'tan':  trig = atan  if exponent == -1 else tan
        if any(self.peek(token) for token in
                ('RATIONAL', 'DECIMAL', 'INTEGER', 'PI')):
            expr = self._number()
        elif any(self.peek(token) for token in
                ('LETTER', 'DIACRITIC', 'MATHIT')):
            sentence, position = self.lexer.sentence, self.lexer.mark()
            tensor = Tensor(self._tensor(), None)
            symbol, indexing = tensor.symbol, tensor.indexing
            if symbol not in self._namespace:
                if indexing:
                    raise ParseError('cannot index undefined tensor \'%s\' at position %d' %
                        (symbol, position), sentence, position)
                self._define_tensor(Tensor(tensor.function, 0))
            expr = tensor.function
        elif self.accept('LPAREN'):
            expr = self._expression()
            self.expect('RPAREN')
        else:
            sentence, position = self.lexer.sentence, self.lexer.mark()
            raise ParseError('unexpected \'%s\' at position %d' %
                (sentence[position], position), sentence, position)
        if exponent == -1: return trig(expr)
        return trig(expr) ** exponent

    # <OPERATOR> -> [ <VPHANTOM> '{' <MODIFIER> '}' ] ( <PARDRV> | <COVDRV> | <LIEDRV> )
    def _operator(self):
        modifier = self._namespace['deriv']
        if self.accept('VPHANTOM'):
            self.expect('LBRACE')
            _modifier = self.lexer.lexeme
            self.expect('MODIFIER')
            self._namespace['deriv'] = _modifier
            self.expect('RBRACE')
        operator = self.lexer.lexeme
        if self.peek('PAR_SYM'):
            pardrv = self._pardrv()
            self._namespace['deriv'] = modifier
            return pardrv
        if self.peek('COV_SYM') or self.peek('DIACRITIC'):
            covdrv = self._covdrv('RHS')
            self._namespace['deriv'] = modifier
            return covdrv
        if self.peek('LIE_SYM'):
            liedrv = self._liedrv('RHS')
            self._namespace['deriv'] = modifier
            return liedrv
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unsupported operator \'%s\' at position %d' %
            (operator, position), sentence, position)

    # <PARDRV> -> { <PAR_SYM> [ '^' <INTEGER> ] '_' <LETTER> }+ ( <TENSOR> | <SUBEXPR> )
    def _pardrv(self):
        indexing, order = [], 1
        while self.accept('PAR_SYM'):
            if self.accept('CARET'):
                order = self.lexer.lexeme
                sentence, position = self.lexer.sentence, self.lexer.mark()
                self.expect('INTEGER')
                order = int(order)
            self.expect('UNDERSCORE')
            index = self._strip(self.lexer.lexeme)
            self.expect('LETTER')
            indexing.extend(order * [Symbol(index)])
        modifier = self._namespace['deriv']
        if all(index in self._namespace['basis'] for index in indexing):
            modifier = 'symbolic'
        if order > 1 and modifier != 'symbolic':
            raise ParseError('cannot specify derivative order unless symbolic',
                sentence, position)
        if any(self.peek(i) for i in ('LPAREN', 'LBRACK', 'LBRACE')):
            if modifier == 'symbolic':
                derivative = Derivative(self._subexpr(), *indexing)
                return derivative
            if modifier == 'variable':
                tree = ExprTree(self._subexpr())
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
        if modifier == 'symbolic':
            return Derivative(self._tensor(), *indexing)
        sentence, position = self.lexer.sentence, self.lexer.mark()
        function = self._tensor()
        symbol, order = str(function.args[0]), len(indexing)
        indexing = list(function.args[1:]) + indexing
        if symbol not in self._namespace:
            raise ParseError('cannot differentiate undefined tensor \'%s\' at position %d' %
                (symbol, position), sentence, position)
        tensor = self._namespace[symbol]
        symbol = symbol + ('' if '_d' in symbol else '_d') + order * 'D'
        function, rank = Function('Tensor')(symbol, *indexing), len(indexing) - order
        symmetry = 'sym%d%d' % (rank, rank + order - 1) if order == 2 else 'nosym'
        self._define_tensor(Tensor(function, tensor.dimension), symmetry)
        return function

    # <COVDRV> -> { ( <COV_SYM> | <DIACRITIC> '{' <COV_SYM> '}' ) ( '^' | '_' ) <LETTER> }+ ( <TENSOR> | <SUBEXPR> )
    def _covdrv(self, location):
        indexing, equation, diacritic = [], ['', ' = ', '', ''], ''
        sentence, position = self.lexer.sentence, self.lexer.mark()
        alphabet = (chr(97 + n) for n in range(26))
        while self.peek('COV_SYM') or self.peek('DIACRITIC'):
            lexeme = self._strip(self.lexer.lexeme)
            operator = '\\nabla'
            if self.accept('DIACRITIC'):
                diacritic = lexeme
                operator = '\\%s{\\nabla}' % diacritic
                self.expect('LBRACE')
                self.expect('COV_SYM')
                self.expect('RBRACE')
            else: self.expect('COV_SYM')
            metric = self._namespace['metric'][diacritic] + diacritic
            if metric + 'DD' not in self._namespace:
                raise ParseError('cannot generate covariant derivative without defined metric \'%s\'' %
                    metric, sentence, position)
            equation[0] += operator
            equation[3] += operator
            if self.accept('CARET'):
                index = self.lexer.lexeme
                equation[0] += '^' + index + ' '
                bound_index = next(x for x in alphabet if x != index)
                prefix = '\\' if len(self._namespace['metric'][diacritic]) > 1 else ''
                metric = '\\%s{%s}' % (diacritic, prefix + self._namespace['metric'][diacritic]) if diacritic \
                    else prefix + self._namespace['metric'][diacritic]
                equation[2] += '%s^{%s %s} ' % (metric, index, bound_index)
                equation[3] += '_' + bound_index + ' '
                index = self._strip(index)
                self.expect('LETTER')
                indexing.append((Symbol(index), 'U'))
            elif self.accept('UNDERSCORE'):
                index = _index = self.lexer.lexeme
                _indexing = [str(index[0]) for index in indexing]
                if _index in _indexing:
                    _index = next(x for x in alphabet if x not in _indexing)
                equation[0] += '_' + _index + ' '
                equation[3] += '_' + _index + ' '
                index = self._strip(index)
                self.expect('LETTER')
                indexing.append((Symbol(index), 'D'))
            else:
                sentence, position = self.lexer.sentence, self.lexer.mark()
                raise ParseError('unexpected \'%s\' at position %d' %
                    (sentence[position], position), sentence, position)
        marker_1 = self.lexer.mark()
        function = self._tensor()
        marker_2 = self.lexer.index - len(self.lexer.lexeme)
        equation[0] += self.lexer.sentence[marker_1:marker_2].strip()
        equation[3] += self.lexer.sentence[marker_1:marker_2].strip()
        if location == 'RHS':
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
                        symbol, sentence, position)
                tensor = Tensor(function, self._namespace[symbol].dimension)
                self.parse(self._generate_covdrv(tensor, indexing, diacritic))
                self.lexer.initialize(sentence, position)
                self.lexer.lex()
        symbol, suffix = str(function.args[0]), ''.join([index[1] for index in indexing])
        symbol = symbol + ('' if '_cd' in symbol else '_cd' + diacritic) + suffix
        indexing = list(function.args[1:]) + [index[0] for index in indexing]
        return Function('Tensor')(symbol, *indexing)

    # <LIEDRV> -> <LIE_SYM> '_' <SYMBOL> ( <TENSOR> | <SUBEXPR> )
    def _liedrv(self, location):
        self.expect('LIE_SYM')
        self.expect('UNDERSCORE')
        vector = self._strip(self._symbol())
        function = self._tensor()
        if location == 'RHS':
            sentence, position = self.lexer.sentence, self.lexer.mark()
            symbol = str(function.args[0])
            if symbol not in self._namespace:
                sentence, position = self.lexer.sentence, self.lexer.mark()
                raise ParseError('cannot differentiate undefined tensor \'%s\'' %
                    symbol, sentence, position)
            tensor = Tensor(function, self._namespace[symbol].dimension)
            tensor.weight = self._namespace[symbol].weight
            self.parse(self._generate_liedrv(tensor, vector))
            self.lexer.initialize(sentence, position)
            self.lexer.lex()
        symbol = str(function.args[0]) + '_ld' + vector
        return Function('Tensor')(symbol, *function.args[1:])

    # <NUMBER> -> <RATIONAL> | <DECIMAL> | <INTEGER> | <PI>
    def _number(self):
        number = self.lexer.lexeme
        if self.accept('RATIONAL'):
            rational = re.match(r'(\-?[0-9]+)\/(\-?[1-9][0-9]*)', number)
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
                        symbol, sentence, position)
                tensor = self._namespace[symbol]
                symbol = ''.join(symbol) + ('' if '_d' in ''.join(symbol) else '_d') + order * 'D'
                function, rank = Function('Tensor')(symbol, *indexing), len(indexing) - order
                symmetry = 'sym%d%d' % (rank, rank + order - 1) if order == 2 else 'nosym'
                self._define_tensor(Tensor(function, tensor.dimension), symmetry)
                return function
        self.lexer.mark()
        if self.accept('CARET'):
            if self.accept('LBRACE'):
                if self.accept('LBRACE'):
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

    # <SYMBOL> -> <LETTER> | <DIACRITIC> '{' <LETTER> '}' | <MATHIT> '{' <LETTER> { '_' | <LETTER> | <INTEGER> }* '}'
    def _symbol(self):
        lexeme = self.lexer.lexeme
        if self.accept('LETTER'):
            return lexeme
        if self.accept('DIACRITIC'):
            self.expect('LBRACE')
            symbol = self.lexer.lexeme + lexeme[1:]
            self.expect('LETTER')
            self.expect('RBRACE')
            return symbol
        if self.accept('MATHIT'):
            self.expect('LBRACE')
            symbol = [self.lexer.lexeme]
            if not self.accept('LETTER'):
                sentence, position = self.lexer.sentence, self.lexer.mark()
                raise ParseError('unexpected \'%s\' at position %d' %
                    (sentence[position], position), sentence, position)
            while any(self.peek(token) for token in
                    ('UNDERSCORE', 'LETTER', 'INTEGER')):
                symbol.extend([self.lexer.lexeme])
                self.lexer.lex()
            self.expect('RBRACE')
            return ''.join(symbol).replace('\\', '')
        sentence, position = self.lexer.sentence, self.lexer.mark()
        raise ParseError('unexpected \'%s\' at position %d' %
            (sentence[position], position), sentence, position)

    # <LOWER_INDEX> -> <LETTER> | <INTEGER> | '{' { <LETTER> | <INTEGER> }* [ ( ',' | ';' ) { <LETTER> }+ ] '}'
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
        if self.accept('LBRACE'):
            while self.peek('LETTER') or self.peek('INTEGER'):
                append_index()
            # TODO: SEMICOLON NOTATION FOR COVARIANT DERIVATIVE
            if self.accept('COMMA'):
                while self.peek('LETTER'):
                    order += 1
                    append_index()
            elif self.accept('SEMICOLON'):
                covariant = True
                while self.peek('LETTER'):
                    order += 1
                    append_index()
            self.expect('RBRACE')
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
        if self.accept('LBRACE'):
            while self.peek('LETTER') or self.peek('INTEGER'):
                append_index()
            self.expect('RBRACE')
            return indexing
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
            dimension = tensor.dimension
            if not symmetry and permutation:
                # instantiate permutation (Levi-Civita) symbol using parity
                index  = [chr(105 + n) for n in range(rank)]
                prefix = '[' * rank + 'sgn([' + ', '.join(index) + '])'
                suffix = ''.join(' for %s in range(%d)]' % (index[rank - i], dimension) for i in range(1, rank + 1))
                tensor.structure  = eval(prefix + suffix, {'sgn': sgn})
            elif not symmetry and kronecker:
                if rank != 2:
                    raise TensorError('cannot instantiate kronecker delta of rank ' + str(rank))
                tensor.structure = ixp.declare_indexedexp(rank=rank, dimension=dimension)
                for i in range(dimension): tensor.structure[i][i] = 1
            elif rank == 0:
                tensor.structure = Symbol(symbol)
            else:
                symmetry = 'sym01' if symmetry == 'metric' \
                      else None    if symmetry == 'nosym' \
                      else symmetry
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
    def _generate_christoffel(function, metric):
        symbol, indexing = '\\' + str(function.args[0])[:-3], function.args[1:]
        diacritic = 'bar'   if 'bar'   in symbol \
               else 'hat'   if 'hat'   in symbol \
               else 'tilde' if 'tilde' in symbol \
               else ''
        prefix = '\\' if len(metric[diacritic]) > 1 else ''
        metric = '\\%s{%s}' % (diacritic, prefix + metric[diacritic]) if diacritic \
            else prefix + metric[diacritic]
        if diacritic: symbol = '\\%s{%s}' % (diacritic, symbol[:-len(diacritic)])
        indexing = [('\\' if len(str(index)) > 1 else '') + str(index) for index in indexing]
        bound_index = next(x for x in (chr(97 + n) for n in range(26)) if x not in indexing)
        return (('{symbol}^{i1}_{{{i2}{i3}}} = \\frac{{1}}{{2}} {metric}^{{{i1} {bound_index}}}(\\partial_{i2} {metric}_{{{i3} {bound_index}}} + \\partial_{i3} {metric}_{{{bound_index} {i2}}} - \\partial_{bound_index} {metric}_{{{i2} {i3}}})')
                .format(i1 = indexing[0], i2 = indexing[1], i3 = indexing[2], symbol = symbol, metric = metric, bound_index = bound_index))

    @staticmethod
    def _generate_covdrv(tensor, indexing, diacritic=''):
        alphabet = (chr(97 + n) for n in range(26))
        indexing, order = tensor.indexing + indexing[::-1], len(indexing)
        indexing = [(str(idx), str(pos)) for idx, pos in indexing]
        for i, (index, position) in enumerate(indexing):
            _indexing = [index for index, _ in indexing]
            if index in _indexing[:i]:
                index = next(x for x in alphabet if x not in _indexing)
                indexing[i] = (index, position)
        LHS, pos = '', lambda x: -(order + 1) + x
        for deriv_index, _ in indexing[:-(order + 1):-1]:
            if len(deriv_index) > 1:
                deriv_index = '\\' + deriv_index
            LHS += ('\\%s{\\nabla}' % diacritic if diacritic else '\\nabla') + ('_%s ' % deriv_index)
        LHS += tensor.latex_format()
        def generate_RHS(symbol, order, indexing):
            if order == 0:
                _indexing = tensor.indexing
                tensor.indexing = indexing[:len(_indexing)]
                latex = tensor.latex_format()
                tensor.indexing = _indexing
                return latex
            deriv_index, _ = indexing[pos(order)]
            if len(deriv_index) > 1:
                deriv_index = '\\' + deriv_index
            latex = generate_RHS(symbol, order - 1, indexing)
            RHS = '\\partial_%s (%s)' % (deriv_index, latex)
            for index, position in indexing[:pos(order)]:
                alphabet = (chr(97 + n) for n in range(26))
                _indexing = [index for index, _ in indexing]
                bound_index = next(x for x in alphabet if x not in _indexing)
                _indexing = list(indexing)
                for i, _ in enumerate(indexing):
                    if index == indexing[i][0]:
                        _indexing[i] = (bound_index, indexing[i][1])
                latex = generate_RHS(symbol, order - 1, _indexing)
                if len(str(index)) > 1:
                    index = '\\' + str(index)
                RHS += ' + ' if position == 'U' else ' - '
                RHS += '\\%s{\\Gamma}' % diacritic if diacritic else '\\Gamma'
                if position == 'U':
                    RHS += '^%s_{%s %s} (%s)' % (index, bound_index, deriv_index, latex)
                else:
                    RHS += '^%s_{%s %s} (%s)' % (bound_index, index, deriv_index, latex)
            return RHS
        return LHS + ' = ' + generate_RHS(tensor.symbol, order, indexing)

    @staticmethod
    def _generate_liedrv(tensor, vector):
        if len(str(vector)) > 1:
            vector = '\\' + str(vector)
        indexing = [str(index[0]) for index in tensor.indexing]
        alphabet = (chr(97 + n) for n in range(26))
        for i, index in enumerate(indexing):
            if index in indexing[:i]:
                indexing[i] = next(x for x in alphabet if x not in indexing)
        LHS = '\\mathcal{L}_%s %s' % (vector, tensor.latex_format())
        bound_index = next(x for x in alphabet if x not in indexing)
        RHS = '%s^%s \\partial_%s %s' % (vector, bound_index, bound_index, tensor.latex_format())
        for index, position in tensor.indexing:
            _tensor = Tensor(tensor.function, tensor.dimension)
            _indexing = [bound_index if i == str(index) else i for i in indexing]
            _tensor.indexing = [(idx, pos) for idx, (_, pos) in zip(_indexing, tensor.indexing)]
            latex = _tensor.latex_format()
            if len(str(index)) > 1:
                index = '\\' + str(index)
            if position == 'U':
                RHS += ' - (\\partial_%s %s^%s) %s' % (bound_index, vector, index, latex)
            else:
                RHS += ' + (\\partial_%s %s^%s) %s' % (index, vector, bound_index, latex)
        if tensor.weight:
            latex = tensor.latex_format()
            RHS += ' + (%s)(\\partial_%s %s^%s) %s' % (tensor.weight, bound_index, vector, bound_index, latex)
        return LHS + ' = ' + RHS

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

    def __init__(self, function, dimension, structure=None, weight=None):
        self.function  = function
        self.dimension = dimension
        self.weight    = weight
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
        diacritic = 'bar'   if 'bar'   in self.symbol \
               else 'hat'   if 'hat'   in self.symbol \
               else 'tilde' if 'tilde' in self.symbol \
               else None
        symbol, indexing = '', list(self.indexing)
        if '_d' in self.symbol:
            _, pos_list = self.symbol.split('_d')
            for _ in pos_list:
                index = str(indexing.pop()[0])
                if len(index) > 1:
                    index = '\\' + index
                symbol += '\\partial_%s ' % index
        elif '_cd' in self.symbol:
            _, pos_list = self.symbol.split('_cd')
            for position in pos_list:
                index = str(indexing.pop()[0])
                if len(index) > 1:
                    index = '\\' + index
                if position == 'U':
                    symbol += '\\nabla^%s ' % index
                else:
                    symbol += '\\nabla_%s ' % index
        latex = [re.split('[UD]', self.symbol)[0], [], []]
        if len(latex[0]) > 1:
            latex[0] = '\\' + latex[0]
        if diacritic:
            latex[0] = '\\%s{%s}' % (diacritic, latex[0][:-len(diacritic)])
        latex[0] = symbol + latex[0]
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
        if self.rank == 0:
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
    for kwrd in ('basis', 'deriv', 'index', 'ignore', 'metric'):
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
