""" NRPy+ LaTeX to Sympy Parser """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

from sympy.parsing.sympy_parser import parse_expr
import re

class Lexer:
	""" LaTeX Lexer

		The following class will tokenize an expression for usage in parsing.
	"""

	def __init__(self):
		self.grammar = { r'(?:[0-9]+\/[1-9]+)|(?:\\frac{[0-9]+}{[1-9]+})' : 'RATIONAL',
						 r'[0-9]+\.[0-9]+' : 'DECIMAL',
						 r'[1-9][0-9]*'    : 'INTEGER',
						 r'\+'			   : 'PLUS',
						 r'\-'			   : 'MINUS',
						 r'\/'			   : 'DIVIDE',
						 r'\^'			   : 'SUPERSCRIPT',
						 r'\_'			   : 'SUBSCRIPT',
						 r'\('			   : 'LEFT_PAREN',
						 r'\)'			   : 'RIGHT_PAREN',
						 r'\{'			   : 'LEFT_BRACE',
						 r'\}'			   : 'RIGHT_BRACE',
						 r'\['			   : 'LEFT_BRACKET',
						 r'\]'			   : 'RIGHT_BRACKET',
						 r'\\'	   		   : 'COMMAND',
						 r'sqrt'           : 'CMD_SQRT',
						 r'[a-zA-Z]'       : 'SYMBOL' }
		self.regex = re.compile('|'.join(['(?P<%s>%s)' % \
			(self.grammar[pattern], pattern) for pattern in self.grammar]))
	
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
			if self.sentence[self.index].isspace():
				self.index += 1; continue
			if not token:
				raise RuntimeError('Unexpected \'' + self.sentence[self.index] + '\'')
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

		The following class will parse an expression according to a defined grammar.
	
		NRPy+ LaTeX Grammar
		<EXPRESSION> -> { - } <TERM> { ( + | - ) <TERM> }
		<TERM>		 -> <FACTOR> { { ( / | ^ | _ ) } <FACTOR> }
		<FACTOR>	 -> <SYMBOL> | <NUMBER> | <COMMAND> | \(<EXPRESSION>\)
		<SYMBOL>	 -> a | ... | z | A | ... | Z |
		<NUMBER>     -> <RATIONAL> | <DECIMAL> | <INTEGER>
		<COMMAND>    -> \ ( <SQRT> | ... )
		<SQRT>		 -> sqrt { [<INTEGER>] } \{<EXPRESSION>\}
	"""

	def __init__(self):
		self.lexer = Lexer()

	def parse(self, sentence):
		""" Parse Sentence

			:arg:    sentence (raw string)
			:return: parsed sentence
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
		while self.__peek('DIVIDE'):
			operator = self.lexer.word
			self.lexer.lex()
			expr += operator + self.__factor()
		while any(self.__peek(i) for i in ('COMMAND', 'LEFT_PAREN', \
				'SYMBOL', 'RATIONAL', 'DECIMAL', 'INTEGER')):
			expr += '*' + self.__factor()
		return expr
	
	def __factor(self):
		literal = self.lexer.word
		if self.__accept('SYMBOL'):
			return literal
		elif any(self.__peek(i) for i in ('RATIONAL', 'DECIMAL', 'INTEGER')):
			return self.__number()
		elif self.__accept('COMMAND'):
			return self.__command()
		elif self.__accept('LEFT_PAREN'):
			expr = '(' + self.__expression() + ')'
			self.__expect('RIGHT_PAREN')
			return expr
	
	def __number(self):
		number = self.lexer.word
		if self.__accept('RATIONAL'):
			rational = re.match(r'([1-9][0-9]*)\/([1-9][0-9]*)', number)
			if not rational:
				rational = re.match(r'\\frac{([0-9]+)}{([1-9]+)}', number)
			number = 'Rational(%s, %s)' % (rational.group(1), rational.group(2))
		elif self.__accept('DECIMAL'):
			number = 'Float(%s)' % number
		else: self.__expect('INTEGER')
		return number
	
	def __command(self):
		if self.__accept('CMD_SQRT'):
			return self.__sqrt()
	
	def __sqrt(self):
		if self.__accept('LEFT_BRACKET'):
			root = self.__number()
			self.__expect('RIGHT_BRACKET')
		else: root = 2
		self.__expect('LEFT_BRACE')
		expr = self.__expression()
		self.__expect('RIGHT_BRACE')
		return 'Pow(%s, Rational(1, %s))' % (expr, root)
	
	def __peek(self, token_type):
		return self.lexer.token == token_type

	def __accept(self, token_type):
		if self.__peek(token_type):
			self.lexer.lex()
			return True
		return False
	
	def __expect(self, token_type):
		if not self.__accept(token_type):
			raise RuntimeError('Expected \'' + token_type + '\'')

def parse(expression):
	""" Convert LaTeX Expression to SymPy Expression

		:arg:    LaTeX Expression
		:return: SymPy Expression

		>>> from latex_parser import parse
		>>> parse(r'-a(b - \\frac{2}{3}) + \\sqrt[5]{3}')
		-a*(b - 2/3) + 3**(1/5)
	"""
	return Parser().parse(expression)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
	