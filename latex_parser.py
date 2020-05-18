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

	def __init__(self):
		self.lexer = Lexer()

	def parse(self, sentence):
		self.lexer.initialize(sentence)
		self.lexer.lex()
		return parse_expr(self.expression())
	
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
	