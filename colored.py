import sys

def not_colored(a,_):
    return repr(a)

colors = {
  "red":"\033[31m",
  "green":"\033[32m",
  "yellow":"\033[33m",
  "blue":"\033[34m",
  "magenta":"\033[35m",
  "cyan":"\033[36m",
}
reset = "\033[0m"

def colored(arg,c):
    assert isinstance(c, str)
    assert c in colors
    s = str(arg)
    return colors[c] + s + reset

if hasattr(sys.stdout,"isatty"):
    is_tty = sys.stdout.isatty()
else:
    is_tty = False

is_jupyter = type(sys.stdout).__name__ == 'OutStream' and  type(sys.stdout).__module__ == 'ipykernel.iostream'
if (not is_tty) and (not is_jupyter):
    colored = not_colored

if __name__ == "__main__":
    if installed:
        print(colored("Colored was installed", "green"))
    else:
        print("Colored was NOT installed")
