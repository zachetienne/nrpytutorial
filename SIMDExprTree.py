__author__ = 'Ken Sible'

class ExprTree:
    """ SymPy Expression Tree (Nodal Structure) 
    
        >>> tree = ExprTree(cos(a + b)**2)
        >>> print(tree)
        ExprTree(cos(a + b)**2)
        >>> repr(tree) # Preorder (Default)
        [cos(a + b)**2, cos(a + b), a + b, a, b, 2]
        >>> tree.replace({a + b: x})
        cos(x)**2
        >>> [str(subtree.expr.func).split('.')[-1][:-2] \
                 for subtree in tree.preorder()]
        ['Pow', 'cos', 'Symbol', 'Integer']
    """

    def __init__(self, expr):
        self.root = self.Node(expr)
        self.build(self.root)
    
    def build(self, node, clear=False):
        """ Build expression (sub)tree.

            :arg:   root node of (sub)tree
            :arg:   clear children (default: False)

            >>> root = tree.Node(cos(a + b)**2)
            >>> tree.build(root)
            >>> repr(tree)
            [cos(a + b)**2, cos(a + b), a + b, a, b, 2]
            >>> tree.root.expr = sin(ab)**2
            >>> tree.build(root, clear=True)
            >>> repr(tree)
            [sin(ab)**2, sin(ab), ab, a, b, 2]
        """
        if clear: node.children.clear()
        for arg in node.expr.args:
            subtree = self.Node(arg)
            node.append(subtree)
            self.build(subtree)

    def preorder(self, node=None):
        """ Generate iterator for preorder traversal.

            :arg:    root node of (sub)tree
            :return: iterator

            >>> tree = ExprTree(cos(ab)**2)
            >>> for i, subtree in enumerate(tree.preorder()):
                    if subtree.expr.func == Mul:
                        print((i, subtree.expr))
            (2, ab)
        """
        if node == None:
            node = self.root
        yield node
        for child in node.children:
            for subtree in self.preorder(child):
                yield subtree
    
    def postorder(self, node=None):
        """ Generate iterator for postorder traversal.

            :arg:    root node of (sub)tree
            :return: iterator

            >>> tree = ExprTree(cos(ab)**2)
            >>> for i, subtree in enumerate(tree.postorder()):
                    if subtree.expr.func == Mul:
                        print((i, subtree.expr))
            (2, ab)
        """
        if node == None:
            node = self.root
        for child in node.children:
            for subtree in self.postorder(child):
                yield subtree
        yield node
    
    def replace(self, rule):
        """
        Replace subexpression(s) from dictionary mapping.

        :arg:    replacement dictionary {original: substitution}
        :return: root expression

        >>> tree = ExprTree(cos(a + b)**2)
        >>> tree.replace({a + b: x})
        cos(x)**2
        >>> tree = ExprTree(a*b + c)
        >>> tree.replace({a*b: x, c: y})
        x + y
        """
        from sympy.parsing.sympy_parser import parse_expr
        for expr in rule:
            self.root.expr = parse_expr(str(self.root.expr).\
                replace(str(expr), str(rule[expr])))
        self.build(self.root, clear=True)
        return self.root.expr
    
    def reconstruct(self, evaluate=True):
        """
        Reconstruct root expression from expression tree.

        :arg:    evaluate root expression
        :return: root expression

        >>> tree = ExprTree(cos(a + b)**2)
        >>> tree.root.children[0] = sin(ab)
        >>> tree.reconstruct()
        >>> print(tree)
        ExprTree(sin(a + b)**2)
        """
        for subtree in self.postorder(self.root):
            if subtree.children:
                expr_list = [node.expr for node in subtree.children]
                subtree.expr = subtree.expr.func(\
                    *expr_list, evaluate=evaluate)
        return self.root.expr

    class Node:
        def __init__(self, expr):
            self.expr = expr
            self.children = []

        def append(self, node):
            self.children.append(node)
    
    def __repr__(self):
        return str([node.expr for node in self.preorder(self.root)])

    def __str__(self):
        return f'ExprTree({self.root.expr})'