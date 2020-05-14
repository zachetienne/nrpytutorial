import re

class Lexer:

	def __init__(self):
		self.grammar = { r'[a-zA-Z]'       : 'SYMBOL',
						 r'(?:[0-9]+\/[1-9]+)|(?:\\frac{[0-9]+}{[1-9]+})' : 'RATIONAL',
						 r'[0-9]+\.[0-9]+' : 'DECIMAL',
						 r'[1-9][0-9]*'    : 'INTEGER',
						 r'\+'			   : 'ADDITION',
						 r'\-'			   : 'SUBTRACTION',
						 r'\/'			   : 'DIVISION',
						 r'\^'			   : 'SUPERSCRIPT',
						 r'\_'			   : 'SUBSCRIPT',
						 r'\('			   : 'LEFT_PAREN',
						 r'\)'			   : 'RIGHT_PAREN',
						 r'{'			   : 'LEFT_BRACKET',
						 r'}'			   : 'RIGHT_BRACKET',
						 r'\\[a-zA-Z]+'	   : 'COMMAND' }
		self.regex = re.compile('|'.join(['(?P<%s>%s)' % \
			(self.grammar[pattern], pattern) for pattern in self.grammar]))
	
	def initialize(self, sentence):
		""" Initialize Lexer

			:arg: sentence (raw string)
		"""
		self.sentence = sentence
		self.index = 0

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
			yield token.lastgroup
		yield None
	
	def token(self):
		""" Retrieve Current Token

			:return: next token in iterator
		"""
		return next(self.tokenize())

class Parser:

	def parse(self, sentence):
		self.lexer = Lexer()
		self.lexer.initialize(sentence)
		self.token = self.lexer.tokenize()

	def peek(self, token_type):
		return self.token == token_type

	def accept(self, token_type):
		if self.peek(token_type):
			self.token = self.tokenize()
			return True
		return False
	
	def expect(self, token_type):
		if not self.accept(token_type):
			raise RuntimeError('Expected \'' + token_type + '\'')
	