""" SymPy (N-Ary) Expression Tree

The following script will extend the expression tree from SymPy, 
allowing direct node manipulation for subexpression replacement.
The expression tree structure within SymPy expressions stores
subexpressions inside immutable tuples, preventing the client from
modifying the expression tree. Therefore, the client must depend on
build-in functions, such as xreplace, for subexpression replacement,
which might be suboptimal for their specific purpose. The ExprTree class
is implemented as an n-ary tree data structure for SymPy expressions,
equipped with a build method for constructing the expression tree, 
a reconstruct method for reconstructing the root expression, a replace
method for subexpression replacement, and preorder/postorder traversal
iterators (or generators). The __repr__ representation of the expression
tree will return a string of the expressions using the preorder traversal,
while the __str__ representation will return a string of the class name 
and root expression. The Node subclass has a field for an expression and
a field for subexpression children (implemented as a mutable list).
"""
# Author: Ken Sible
# Email:  ksible@outlook.com

__author__ = 'Ken Sible'

class ExprTree:
    """ SymPy (N-Ary) Expression Tree
    
        >>> from sympy.abc import a, b, x
        >>> from sympy import cos
        >>> tree = ExprTree(cos(a + b)**2)
        >>> print(tree)
        ExprTree(cos(a + b)**2)
        >>> repr(tree)
        '[cos(a + b)**2, cos(a + b), a + b, a, b, 2]'
        >>> tree.replace({a + b: x})
        cos(x)**2
    """

    def __init__(self, expr):
        self.root = self.Node(expr)
        self.build(self.root)
    
    def build(self, node, clear=False):
        """ Build expression (sub)tree.

            :arg:   root node of (sub)tree
            :arg:   clear children (default: False)

            >>> from sympy.abc import a, b
            >>> from sympy import cos, sin
            >>> tree = ExprTree(cos(a + b)**2)
            >>> tree.root.expr = sin(a*b)**2
            >>> tree.build(tree.root, clear=True)
            >>> repr(tree)
            '[sin(a*b)**2, sin(a*b), a*b, a, b, 2]'
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

            >>> from sympy.abc import a, b
            >>> from sympy import cos, Mul
            >>> tree = ExprTree(cos(a*b)**2)
            >>> for i, subtree in enumerate(tree.preorder()):
            ...     if subtree.expr.func == Mul:
            ...         print((i, subtree.expr))
            (2, a*b)
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

            >>> from sympy.abc import a, b
            >>> from sympy import cos, Mul
            >>> tree = ExprTree(cos(a*b)**2)
            >>> for i, subtree in enumerate(tree.postorder()):
            ...     if subtree.expr.func == Mul:
            ...         print((i, subtree.expr))
            (2, a*b)
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

        >>> from sympy.abc import a, b, x
        >>> from sympy import cos
        >>> tree = ExprTree(a*b + cos(a + b)**2)
        >>> tree.replace({a*b: x, a + b: x})
        x + cos(x)**2
        """
        from sympy.parsing.sympy_parser import parse_expr
        for expr in rule:
            self.root.expr = parse_expr(str(self.root.expr).\
                replace(str(expr), str(rule[expr])))
        self.build(self.root, clear=True)
        return self.root.expr
    
    def reconstruct(self, evaluate=False):
        """
        Reconstruct root expression from expression tree.

        :arg:    evaluate root expression (default: False)
        :return: root expression

        >>> from sympy.abc import a, b
        >>> from sympy import cos, sin
        >>> tree = ExprTree(cos(a + b)**2)
        >>> tree.root.children[0].expr = sin(a + b)
        >>> tree.reconstruct()
        sin(a + b)**2
        """
        for subtree in self.postorder():
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
        return str([node.expr for node in self.preorder()])

    def __str__(self):
        return 'ExprTree(%s)' % str(self.root.expr)

if __name__ == "__main__":
    import doctest
    doctest.testmod()