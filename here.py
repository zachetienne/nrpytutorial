from __future__ import print_function
from colored import colored
import os
import re
import sys

_here = os.path.realpath(os.getcwd())


def here(*args):
    herell(False, *args)


def herecc(*args):
    herell(True, *args)


def herell(usecc, *args):
    import inspect
    stack = inspect.stack()
    frame = stack[2]
    if usecc:
        herestr = re.sub(r"^herecc\((.*)\)$", r"HERE: \1:", frame.code_context[0].strip())
    else:
        herestr = "HERE:"
    if isinstance(frame, tuple):
        frame = frame[0]
        fname = _here
        line = frame.f_lineno
    else:
        fname = os.path.realpath(frame.filename)
        line = frame.lineno
    if fname.startswith(_here):
        fname = fname[len(_here) + 1:]
    assert isinstance(fname, str)
    assert isinstance(line, int)
    nargs = [colored(herestr, "cyan"), fname + ":" + colored(line, "yellow")] + list(args)  # , flush=True)
    print(*nargs)
    sys.stdout.flush()
    sys.stderr.flush()
    frame = None
    stack = None


if __name__ == "__main__":
    here(_here)
    herecc(_here)
