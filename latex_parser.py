from sympy.parsing.sympy_parser import parse_expr
import re

class Lexer:

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
						 r'int'			   : 'CMD_INT',
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
	""" NRPy+ LaTeX Grammar
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
		self.lexer.initialize(sentence)
		self.lexer.lex()
		return parse_expr(self.expression())

	def expression(self):
		sign = '-' if self.accept('MINUS') else ''
		expr = sign + self.term()
		while self.peek('PLUS') or self.peek('MINUS'):
			operator = self.lexer.word
			self.lexer.lex()
			expr += operator + self.term()
		return expr

	def term(self):
		expr = self.factor()
		while self.peek('DIVIDE'):
			operator = self.lexer.word
			self.lexer.lex()
			expr += operator + self.factor()
		while any(self.peek(i) for i in ('COMMAND', 'LEFT_PAREN', \
				'SYMBOL', 'RATIONAL', 'DECIMAL', 'INTEGER')):
			expr += '*' + self.factor()
		return expr
	
	def factor(self):
		literal = self.lexer.word
		if self.accept('SYMBOL'):
			return literal
		elif any(self.peek(i) for i in ('RATIONAL', 'DECIMAL', 'INTEGER')):
			return self.number()
		elif self.accept('COMMAND'):
			return self.command()
		elif self.accept('LEFT_PAREN'):
			expr = '(' + self.expression() + ')'
			self.expect('RIGHT_PAREN')
			return expr
	
	def number(self):
		number = self.lexer.word
		if self.accept('RATIONAL'):
			rational = re.match(r'([1-9][0-9]*)\/([1-9][0-9]*)', number)
			if not rational:
				rational = re.match(r'\\frac{([0-9]+)}{([1-9]+)}', number)
			number = 'Rational(%s, %s)' % (rational.group(1), rational.group(2))
		elif self.accept('DECIMAL'):
			number = 'Float(%s)' % number
		else: self.expect('INTEGER')
		return number
	
	def command(self): pass
	
	def peek(self, token_type):
		return self.lexer.token == token_type

	def accept(self, token_type):
		if self.peek(token_type):
			self.lexer.lex()
			return True
		return False
	
	def expect(self, token_type):
		if not self.accept(token_type):
			raise RuntimeError('Expected \'' + token_type + '\'')
	