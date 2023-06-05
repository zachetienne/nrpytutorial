""" Functional Programming Toolkit """
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

import sys


def pipe(x, *f):
    """ Pipe Operator

        >>> pipe(range(5, 0, -1), reversed, list)
        [1, 2, 3, 4, 5]

        >>> pipe([3, 2, 2, 4, 5, 1], sorted, set, list)
        [1, 2, 3, 4, 5]
    """
    if not f: return x
    return pipe(f[0](x), *f[1:])


def repeat(f, x, n):
    """ Repeat Function

        >>> list(repeat(flatten, [1, 2, [3, [4]], 5], 2))
        [1, 2, 3, 4, 5]
    """
    if n == 0: return x
    return repeat(f, f(x), n - 1)


def chain(*iterable):
    """ Chain Iterable(s)

        >>> list(chain([1], [2, 3], [4, 5]))
        [1, 2, 3, 4, 5]
    """
    for iter_ in iterable:
        try:
            iter(iter_)
        except TypeError:
            iter_ = [iter_]
        for element in iter_:
            yield element


def flatten(iterable):
    """ Flatten Iterable

        >>> list(flatten([1, [2, 3], [4, 5]]))
        [1, 2, 3, 4, 5]
    """
    return chain(*iterable)


def reduce(f, iterable, initializer=None):
    """ Reduction Operation

        >>> reduce(lambda x, y: x + y, [1, 2, 3, 4, 5])
        15

        >>> reduce(lambda x, y: x + y, ['w', 'o', 'r', 'd'])
        'word'

        >>> x = [1, 2, [3, 4], 'aabb']
        >>> reduce(lambda i, _: i + 1, [1] + x[1:])
        4
    """
    iterable = iter(iterable)
    result = next(iterable) if initializer is None \
        else initializer
    for element in iterable:
        result = f(result, element)
    return result


def uniquify(iterable):
    """ Uniquify Iterable

        >>> uniquify(([1, 1, 2, 3, 3, 3, 4, 5, 5]))
        [1, 2, 3, 4, 5]
    """
    return reduce(lambda l, x: l if x in l else l + [x], iterable, [])


def product(*iterable, **kwargs):
    """ Cartesian Product

        >>> list(product(['a', 'b'], [1, 2, 3]))
        [('a', 1), ('a', 2), ('a', 3), ('b', 1), ('b', 2), ('b', 3)]

        >>> list(product([1, 2, 3], repeat=2))
        [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)]

        >>> for i, j in product(['a', 'b'], range(1, 3)):
        ...     print('%s: %d' % (i, j))
        a: 1
        a: 2
        b: 1
        b: 2
    """
    if 'repeat' in kwargs:
        if kwargs['repeat'] > 1 and len(iterable) == 1:
            iterable = kwargs['repeat'] * iterable
    f = lambda A, B: [list(flatten([x] + [y])) for x in A for y in B]
    for prod in reduce(f, iterable):
        yield tuple(prod)


if __name__ == "__main__":
    import doctest

    sys.exit(doctest.testmod()[0])
