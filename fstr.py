from inspect import currentframe
import re
import sys
from here import here

def f(s):
    """
    Mimic the functionality of formatted strings in Python3. Convert curly brackets in s
    to expressions.
    >>> f('3+2={3+2}.')
    '3+2=5.'
    >>> f('3+2={3+2}')
    '3+2=5'
    >>> f('{"="*3} test {"="*3}')
    '=== test ==='
    >>> f('{{hello}}')
    '{hello}'

    >>> metric="gxx"
    >>> f('\\texttt{{{metric}}}')
    '\\texttt{gxx}'
    """
    globs = currentframe().f_back.f_globals
    locs= currentframe().f_back.f_locals
    count = 0
    ns = u''
    w = ''
    i = 0
    while i < len(s):
        c = s[i]
        if i + 1 < len(s):
            nc = s[i+1]
        else:
            nc = ""
        i += 1

        if c == '{' and nc == '{':
            ns += '{'
            i += 1
        elif c == '}' and nc == '}':
            ns += '}'
            i += 1
        elif c == '{':
            count = 1
            j = i
            while i < len(s):
                if s[i] == '{':
                    count += 1
                elif s[i] == '}':
                    count -= 1
                if count == 0:
                    break
                i += 1
            ns += str(eval(s[j:i],globs,locs))
            i += 1
        else:
           ns += c
    return ns

if __name__ == "__main__":
    import doctest
    doctest.testmod()
