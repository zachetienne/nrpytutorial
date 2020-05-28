""" Convert LaTeX Sentence to SymPy Expression """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

from sympy import Symbol, Integer, Rational, Float, Pow
from sympy import var, sqrt
import re

class Lexer:
    """ LaTeX Lexer

        The following class will tokenize a sentence for parsing.
    """

    def __init__(self):
        greek = '|'.join([r'\\[aA]lpha', r'\\[bB]eta', r'\\[gG]amma', r'\\[dD]elta',
            r'\\[eE]psilon', r'\\[zZ]eta', r'\\[eE]ta', r'\\[tT]heta', r'\\[iI]ota',
            r'\\[kK]appa', r'\\[lL]ambda', r'\\[mM]u', r'\\[nN]u', r'\\[xX]i',
            r'\\[oO]mikron', r'\\[pP]i', r'\\[Rr]ho', r'\\[sS]igma', r'\\[tT]au',
            r'\\[uU]psilon', r'\\[pP]hi', r'\\[cC]hi', r'\\[pP]si', r'\\[oO]mega'])
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
              ('SPACE_DELIM',    r'\s+|(?:\\,)+'),
              ('SQRT_CMD',       r'\\sqrt'),
              ('FRAC_CMD',       r'\\frac'),
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
                raise ParseError('unexpected \'%s\' at position %d' % \
                    (self.sentence[self.index], self.index), self.sentence, self.index)
            elif token.lastgroup in ('BIGL_DELIM', 'BIGR_DELIM', \
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

class Parser:
    """ LaTeX Parser

        The following class will parse a tokenized sentence.
    
        LaTeX Grammar:
        <ROOT>          -> <VARIABLE> = <EXPR> | <EXPR>
        <EXPR>          -> [ - ] <TERM> { ( + | - ) <TERM> }
        <TERM>          -> <FACTOR> { [ ( / | ^ ) ] <FACTOR> }
        <FACTOR>        -> <OPERAND> | (<EXPR>) | [<EXPR>]
        <OPERAND>       -> <VARIABLE> | <NUMBER> | <COMMAND>
        <VARIABLE>      -> <SYMBOL> [ _( <SYMBOL> | <INTEGER> ]
        <NUMBER>        -> <RATIONAL> | <DECIMAL> | <INTEGER>
        <COMMAND>       -> <SQRT> | <FRAC>
        <SQRT>          -> \ sqrt [ [<INTEGER>] ] {<EXPR>}
        <FRAC>          -> \ frac {<EXPR>} {<EXPR>}
    """

    def __init__(self):
        self.lexer = Lexer()

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
            raise ParseError('unexpected \'%s\' at position %d' % \
                (sentence[position], position), sentence, position)
        return root

    def __root(self):
        if self.__peek('SYMBOL'):
            variable = eval(self.__variable())
            if self.__accept('EQUAL'):
                expr = eval(self.__expr())
                return {variable: expr}
            self.lexer.reset()
            self.lexer.lex()
        return eval(self.__expr())

    def __expr(self):
        sign = '-' if self.__accept('MINUS') else ''
        expr = sign + self.__term()
        while self.__peek('PLUS') or self.__peek('MINUS'):
            operator = self.lexer.lexeme
            self.lexer.lex()
            expr += operator + self.__term()
        return expr

    def __term(self):
        expr = self.__factor()
        while any(self.__peek(i) for i in ('LEFT_PAREN', 'LEFT_BRACKET', \
                'SYMBOL', 'RATIONAL', 'DECIMAL', 'INTEGER', 'DIVIDE', 'CARET',
                'SQRT_CMD', 'FRAC_CMD')):
            operator = self.lexer.lexeme if self.__peek('DIVIDE') \
                else '**' if self.__peek('CARET') else '*'
            if operator != '*': self.lexer.lex()
            expr += operator + self.__factor()
        return expr
    
    def __factor(self):
        if self.__accept('LEFT_PAREN'):
            expr = '(' + self.__expr() + ')'
            self.__expect('RIGHT_PAREN')
            return expr
        elif self.__accept('LEFT_BRACKET'):
            expr = '(' + self.__expr() + ')'
            self.__expect('RIGHT_BRACKET')
            return expr
        return self.__operand()
    
    def __operand(self):
        if self.__peek('SYMBOL'):
            return self.__variable()
        elif any(self.__peek(i) for i in \
            ('SQRT_CMD', 'FRAC_CMD', 'COMMAND')):
            return self.__command()
        return self.__number()
    
    def __variable(self):
        variable = self.lexer.lexeme
        if variable[0] == '\\':
            variable = variable[1:]
        self.__expect('SYMBOL')
        if self.__accept('UNDERSCORE'):
            if self.__peek('SYMBOL') or self.__peek('INTEGER'):
                subscript = self.lexer.lexeme
                if subscript[0] == '\\':
                    subscript = subscript[1:]
                variable += '_' + subscript
                self.lexer.lex(); var(variable)
                return 'Symbol(\'' + variable + '\')'
            sentence = self.lexer.sentence
            position = self.lexer.index - len(self.lexer.lexeme)
            raise ParseError('unexpected \'%s\' at position %d' % \
                (sentence[position], position), sentence, position)
        var(variable)
        return 'Symbol(\'' + variable + '\')'

    def __number(self):
        number = self.lexer.lexeme
        if self.__accept('RATIONAL'):
            rational = re.match(r'([1-9][0-9]*)\/([1-9][0-9]*)', number)
            if not rational:
                rational = re.match(r'\\frac{([0-9]+)}{([1-9]+)}', number)
            return 'Rational(%s, %s)' % (rational.group(1), rational.group(2))
        elif self.__accept('DECIMAL'):
            return 'Float(' + number + ')'
        elif self.__accept('INTEGER'):
            return 'Integer(' + number + ')'
        sentence = self.lexer.sentence
        position = self.lexer.index - len(self.lexer.lexeme)
        raise ParseError('unexpected \'%s\' at position %d' % \
            (sentence[position], position), sentence, position)
    
    def __command(self):
        command = self.lexer.lexeme
        if self.__accept('SQRT_CMD'):
            return self.__sqrt()
        elif self.__accept('FRAC_CMD'):
            return self.__frac()
        position = self.lexer.index - len(self.lexer.lexeme)
        raise ParseError('unsupported command \'%s\' at position %d' % \
            (command, position), self.lexer.sentence, position)
    
    def __sqrt(self):
        if self.__accept('LEFT_BRACKET'):
            root = 'Integer(' + self.lexer.lexeme + ')'
            self.__expect('INTEGER')
            self.__expect('RIGHT_BRACKET')
        else: root = 2
        self.__expect('LEFT_BRACE')
        expr = self.__expr()
        self.__expect('RIGHT_BRACE')
        return 'Pow(%s, Rational(1, %s))' % (expr, root)
    
    def __frac(self):
        self.__expect('LEFT_BRACE')
        numerator = self.__expr()
        self.__expect('RIGHT_BRACE')
        self.__expect('LEFT_BRACE')
        denominator = self.__expr()
        self.__expect('RIGHT_BRACE')
        return '(%s)/(%s)' % (numerator, denominator)
    
    def __peek(self, token_type):
        return self.lexer.token == token_type

    def __accept(self, token_type):
        if self.__peek(token_type):
            self.lexer.lex()
            return True
        return False
    
    def __expect(self, token_type):
        if not self.__accept(token_type):
            position = self.lexer.index - len(self.lexer.lexeme)
            raise ParseError('expected token %s at position %d' % \
                (token_type, position), self.lexer.sentence, position)

class ParseError(Exception):
    """ Invalid LaTeX Sentence """

    def __init__(self, message, sentence, position):
        super().__init__('%s\n%s^\n' % (sentence, (12 + position) * ' ') + message)

def parse(sentence):
    """ Convert LaTeX Sentence to SymPy Expression

        :arg:    LaTeX Sentence
        :return: SymPy Expression

        >>> from latex_parser import parse
        >>> parse(r'-\delta(x^(2n) - \\frac{2}{3})')
        -delta*(x**(2*n) - 2/3)
        >>> parse(r'x_1 = \\sqrt[5]{x + 3}')
        {x_1: (x + 3)**(1/5)}
    """
    return Parser().parse(sentence)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    