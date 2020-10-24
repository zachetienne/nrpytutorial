# As documented in the NRPy+ tutorial notebook
# Tutorial-cmdline_helper.ipynb, this Python script
# provides a multi-platform means to run executables,
# remove files, and compile code.
# *** Note that this is a stripped-down version of
#     cmdline_helper; the full version exists in NRPy+

# Basic functions:
# Execute_input_string(): Executes an input string and redirects
#            output from stdout & stderr to desired destinations.
# delete_existing_files(file_or_wildcard):
#          Runs del file_or_wildcard in Windows, or
#                rm file_or_wildcard in Linux/MacOS

# Authors: Brandon Clark
#          Zach Etienne
#          zachetie **at** gmail **dot* com
#          Kevin Lituchy

import io, os, shlex, subprocess, sys, time


# Execute_input_string(): Executes an input string and redirects
#            output from stdout & stderr to desired destinations.
def Execute_input_string(input_string, file_to_redirect_stdout=os.devnull, verbose=True):

    if verbose:
        print("(EXEC): Executing `"+input_string+"`...")
    start = time.time()
    # https://docs.python.org/3/library/subprocess.html
    if os.name != 'nt':
        args = shlex.split(input_string)
    else:
        args = input_string

    # https://stackoverflow.com/questions/18421757/live-output-from-subprocess-command
    filename = "tmp.txt"
    with io.open(filename, 'w') as writer, io.open(filename, 'rb', buffering=-1) as reader, io.open(file_to_redirect_stdout, 'wb') as rdirect:
        process = subprocess.Popen(args, stdout=rdirect, stderr=writer)
        while process.poll() is None:
            # https://stackoverflow.com/questions/21689365/python-3-typeerror-must-be-str-not-bytes-with-sys-stdout-write/21689447
            sys.stdout.write(reader.read().decode('utf-8'))
            time.sleep(0.2)
        # Read the remaining
        sys.stdout.write(reader.read().decode('utf-8'))
    delete_existing_files(filename)
    end = time.time()
    if verbose:
        print("(BENCH): Finished executing in "+str(end-start)+" seconds.")

# delete_existing_files(file_or_wildcard):
#          Runs del file_or_wildcard in Windows, or
#                rm file_or_wildcard in Linux/MacOS
def delete_existing_files(file_or_wildcard):
    delete_string = ""
    if os.name == "nt":
        delete_string += "del " + file_or_wildcard
    else:
        delete_string += "rm -f " + file_or_wildcard
    os.system(delete_string)

def output_Jupyter_notebook_to_LaTeXed_PDF(notebookname,location_of_template_file=os.path.join("."),verbose=True):
    Execute_input_string(r"jupyter nbconvert --to latex --template "
                         +os.path.join(location_of_template_file,"latex_nrpy_style.tplx")
                         +r" --log-level='WARN' "+notebookname+".ipynb",verbose=False)
    for _i in range(3):  # _i is an unused variable.
        Execute_input_string(r"pdflatex -interaction=batchmode "+notebookname+".tex",verbose=False)
    delete_existing_files(notebookname+".out "+notebookname+".aux "+notebookname+".log")
    if verbose:
        import textwrap
        wrapper = textwrap.TextWrapper(initial_indent="",subsequent_indent="    ",width=75)
        print(wrapper.fill("Created "+notebookname+".tex, and compiled LaTeX file to PDF file "+notebookname+".pdf"))
