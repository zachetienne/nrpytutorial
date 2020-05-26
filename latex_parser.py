""" Convert LaTeX Sentence to SymPy Expression """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

from sympy.parsing.sympy_parser import parse_expr
import re

class Lexer:
    """ LaTeX Lexer

        The following class will tokenize a sentence for parsing.
    """

    def __init__(self):
        self.regex = re.compile('|'.join(['(?P<%s>%s)' % pattern for pattern in 
            [ ('RATIONAL',       r'(?:[0-9]+\/[1-9]+)|(?:\\frac{[0-9]+}{[1-9]+})'),
              ('DECIMAL',        r'[0-9]+\.[0-9]+'),
              ('INTEGER',        r'[1-9][0-9]*'),
              ('PLUS',           r'\+'),
              ('MINUS',          r'\-'),
              ('DIVIDE',         r'\/'),
              ('SUPERSCRIPT',    r'\^'),
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
              ('SPACE_DELIM',    r'(?:\s+)|(?:\\,)+'),
              ('BACKSLASH',      r'\\'),
              ('SQRT_CMD',       r'sqrt'),
              ('FRAC_CMD',       r'frac'),
              ('SYMBOL',         r'[a-zA-Z]') ]]))
    
    def initialize(self, sentence):
        """ Initialize Lexer

            :arg: sentence (raw string)
        """
        self.sentence = sentence
        self.token    = None
        self.word     = None
        self.index    = 0

    def tokenize(self):
        """ Tokenize Sentence

            :return: token iterator
        """
        while self.index < len(self.sentence):
            token = self.regex.match(self.sentence, self.index)
            if token is None:
                raise ParseError('%s\n%s^\n' % (self.sentence, (12 + self.index) * ' ') \
                    + 'unexpected \'%s\' at position %d' % (self.sentence[self.index], self.index))
            elif token.lastgroup in ('BIGL_DELIM', 'BIGR_DELIM', \
                    'LEFT_DELIM', 'RIGHT_DELIM', 'SPACE_DELIM'):
                self.index = token.end()
            else:
                self.index = token.end()
                self.word = token.group()
                yield token.lastgroup
    
    def lex(self):
        """ Retrieve Current Token

            :return: next token in iterator
        """
        try:
            self.token = next(self.tokenize())
        except StopIteration:
            self.token = None
        return self.token

class Parser:
    """ LaTeX Parser

        The following class will parse a tokenized sentence.
    
        LaTeX Grammar:
        <EXPRESSION>    -> [ - ] <TERM> { ( + | - ) <TERM> }
        <TERM>          -> <FACTOR> { [ ( / | ^ ) ] <FACTOR> }
        <FACTOR>        -> <OPERAND> | (<EXPRESSION>) | [<EXPRESSION>]
        <OPERAND>       -> <SYMBOL> | <NUMBER> | <COMMAND>
        <NUMBER>        -> <RATIONAL> | <DECIMAL> | <INTEGER>
        <COMMAND>       -> \ ( <SQRT> | <FRAC> | ... )
        <SQRT>          -> sqrt [ [<INTEGER>] ] {<EXPRESSION>}
        <FRAC>          -> frac {<EXPRESSION>} {<EXPRESSION>}
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
        return parse_expr(self.__expression())

    def __expression(self):
        sign = '-' if self.__accept('MINUS') else ''
        expr = sign + self.__term()
        while self.__peek('PLUS') or self.__peek('MINUS'):
            operator = self.lexer.word
            self.lexer.lex()
            expr += operator + self.__term()
        return expr

    def __term(self):
        expr = self.__factor()
        while any(self.__peek(i) for i in ('LEFT_PAREN', 'LEFT_BRACKET', 'DIVIDE', 'SUPERSCRIPT', \
                'SYMBOL', 'RATIONAL', 'DECIMAL', 'INTEGER', 'COMMAND')):
            operator = self.lexer.word if self.__peek('DIVIDE') \
                else '**' if self.__peek('SUPERSCRIPT') else '*'
            if operator != '*': self.lexer.lex()
            expr += operator + self.__factor()
        return expr
    
    def __factor(self):
        if self.__accept('LEFT_PAREN'):
            expr = '(' + self.__expression() + ')'
            self.__expect('RIGHT_PAREN')
            return expr
        elif self.__accept('LEFT_BRACKET'):
            expr = '(' + self.__expression() + ')'
            self.__expect('RIGHT_BRACKET')
            return expr
        return self.__operand()
    
    def __operand(self):
        operand = self.lexer.word
        if self.__accept('SYMBOL'):
            return operand
        elif self.__accept('BACKSLASH'):
            return self.__command()
        return self.__number()
    
    def __number(self):
        number = self.lexer.word
        if self.__accept('RATIONAL'):
            rational = re.match(r'([1-9][0-9]*)\/([1-9][0-9]*)', number)
            if not rational:
                rational = re.match(r'\\frac{([0-9]+)}{([1-9]+)}', number)
            return 'Rational(%s, %s)' % (rational.group(1), rational.group(2))
        elif self.__accept('DECIMAL'):
            return 'Float(' + str(number) + ')'
        elif self.__accept('INTEGER'):
            return number
        err_pos = self.lexer.index - len(self.lexer.word)
        raise ParseError('%s\n%s^\n' % (self.lexer.sentence, (12 + err_pos) * ' ') \
                + 'unexpected \'%s\' at position %d' % (self.lexer.sentence[err_pos], err_pos))
    
    def __command(self):
        if self.__accept('SQRT_CMD'):
            return self.__sqrt()
        elif self.__accept('FRAC_CMD'):
            return self.__frac()
        err_pos = self.lexer.index - len(self.lexer.word)
        raise ParseError('%s\n%s^\n' % (self.lexer.sentence, (12 + err_pos) * ' ') \
                + 'unsupported command at position %d' % err_pos)
    
    def __sqrt(self):
        if self.__accept('LEFT_BRACKET'):
            root = self.lexer.word
            self.__expect('INTEGER')
            self.__expect('RIGHT_BRACKET')
        else: root = 2
        self.__expect('LEFT_BRACE')
        expr = self.__expression()
        self.__expect('RIGHT_BRACE')
        return 'Pow(%s, Rational(1, %s))' % (expr, root)
    
    def __frac(self):
        self.__expect('LEFT_BRACE')
        numerator = self.__expression()
        self.__expect('RIGHT_BRACE')
        self.__expect('LEFT_BRACE')
        denominator = self.__expression()
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
            err_pos = self.lexer.index - len(self.lexer.word)
            raise ParseError('%s\n%s^\n' % (self.lexer.sentence, (12 + err_pos) * ' ') \
                    + 'expected token %s at position %d' % (token_type, err_pos))

class ParseError(Exception): pass

def parse(sentence):
    """ Convert LaTeX Sentence to SymPy Expression

        :arg:    LaTeX Sentence
        :return: SymPy Expression

        >>> from latex_parser import parse
        >>> parse(r'-a(b^(2a) - \\frac{2}{3}) + \\sqrt[5]{a + 3}')
        -a*(b**(2*a) - 2/3) + (a + 3)**(1/5)
    """
    return Parser().parse(sentence)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    