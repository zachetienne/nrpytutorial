# As documented in the NRPy+ tutorial notebook
# Tutorial-cmdline_helper.ipynb, this Python script
# provides a multi-platform means to run executables,
# remove files, and compile code.

# Basic functions:
# check_executable_exists(): Check to see whether an executable exists.
#                            Error out or return False if not exists;
#                            return True if executable exists in PATH.
# C_compile(): Compile C code using gcc.
# Execute(): Execute generated executable file, using taskset
#            if available. Calls Execute_input_string() to
#            redirect output from stdout & stderr to desired
#            destinations.
# Execute_input_string(): Executes an input string and redirects
#            output from stdout & stderr to desired destinations.
# delete_existing_files(file_or_wildcard):
#          Runs del file_or_wildcard in Windows, or
#                rm file_or_wildcard in Linux/MacOS

# Authors: Brandon Clark
#          Zach Etienne
#          zachetie **at** gmail **dot* com
#          Kevin Lituchy

import io, os, shlex, subprocess, sys, time, multiprocessing, getpass, platform

# check_executable_exists(): Check to see whether an executable exists.
#                            Error out or return False if not exists;
#                            return True if executable exists in PATH.
def check_executable_exists(exec_name,error_if_not_found=True):
    cmd = "where" if os.name == "nt" else "which"
    try:
        subprocess.check_output([cmd, exec_name])
    except subprocess.CalledProcessError:
        if error_if_not_found:
            print("Sorry, cannot execute the command: " + exec_name)
            sys.exit(1)
        else:
            return False
    return True

# C_compile(): Write a function to compile the Main C code into an executable file
def C_compile(main_C_output_path, main_C_output_file, compile_mode="optimized", custom_compile_string="", additional_libraries=""):
    print("Compiling executable...")
    # Step 1: Check for gcc compiler
    check_executable_exists("gcc")

    if additional_libraries != "":
        additional_libraries = " " + additional_libraries

    # Step 2: Delete existing version of executable
    if os.name == "nt":
        main_C_output_file += ".exe"
    delete_existing_files(main_C_output_file)

    # Step 3: Compile the executable
    if compile_mode=="safe":
        compile_string = "gcc -std=gnu99 -O2 -g -fopenmp "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"+additional_libraries
        Execute_input_string(compile_string, os.devnull)
        # Check if executable exists (i.e., compile was successful), if not, try with more conservative compile flags.
        if not os.path.isfile(main_C_output_file):
            # Step 3.A: Revert to more compatible gcc compile option
            print("Most safe failed. Removing -fopenmp:")
            compile_string = "gcc -std=gnu99 -O2 "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"+additional_libraries
            Execute_input_string(compile_string, os.devnull)
        if not os.path.isfile(main_C_output_file):
            print("Sorry, compilation failed")
            sys.exit(1)
    elif compile_mode=="icc":
        check_executable_exists("icc")
        compile_string = "icc -std=gnu99 -O2 -xHost -qopenmp -unroll "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"+additional_libraries
        Execute_input_string(compile_string, os.devnull)
        # Check if executable exists (i.e., compile was successful), if not, try with more conservative compile flags.
        if not os.path.isfile(main_C_output_file):
            print("Sorry, compilation failed")
            sys.exit(1)
    elif compile_mode=="custom":
        Execute_input_string(custom_compile_string, os.devnull)
        # Check if executable exists (i.e., compile was successful), if not, try with more conservative compile flags.
        if not os.path.isfile(main_C_output_file):
            print("Sorry, compilation failed")
            sys.exit(1)
    elif compile_mode=="optimized":
        compile_string = "gcc -std=gnu99 -Ofast -fopenmp -march=native -funroll-loops "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"+additional_libraries
        Execute_input_string(compile_string, os.devnull)
        # Check if executable exists (i.e., compile was successful), if not, try with more conservative compile flags.
        if not os.path.isfile(main_C_output_file):
            # Step 3.A: Revert to more compatible gcc compile option
            print("Most optimized compilation failed. Removing -march=native:")
            compile_string = "gcc -std=gnu99 -Ofast -fopenmp -funroll-loops "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"+additional_libraries
            Execute_input_string(compile_string, os.devnull)
        if not os.path.isfile(main_C_output_file):
            # Step 3.B: Revert to maximally compatible gcc compile option
            print("Next-to-most optimized compilation failed. Moving to maximally-compatible gcc compile option:")
            compile_string = "gcc -std=gnu99 -O2 "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"+additional_libraries
            Execute_input_string(compile_string, os.devnull)
        # Step 3.C: If there are still missing components within the compiler, say compilation failed
        if not os.path.isfile(main_C_output_file):
            print("Sorry, compilation failed")
            sys.exit(1)
    else:
        print("Sorry, compile_mode = \""+compile_mode+"\" unsupported.")
        sys.exit(1)

    print("Finished compilation.")

# Execute(): Execute generated executable file, using taskset
#            if available. Calls Execute_input_string() to
#            redirect output from stdout & stderr to desired
#            destinations.
def Execute(executable, executable_output_arguments="", file_to_redirect_stdout=os.devnull,verbose=True):
    # Step 1: Delete old version of executable file
    if file_to_redirect_stdout != os.devnull:
        delete_existing_files(file_to_redirect_stdout)

    # Step 2: Build the script for executing the desired executable
    execute_string = ""
    # When in Windows...
    # https://stackoverflow.com/questions/1325581/how-do-i-check-if-im-running-on-windows-in-python
    if os.name == "nt":
        # ... do as the Windows do
        # https://stackoverflow.com/questions/49018413/filenotfounderror-subprocess-popendir-windows-7
        execute_prefix = "cmd /c " # Run with cmd /c executable [options] on Windows
    else:
        execute_prefix = "./"      # Run with ./executable [options] on Linux & Mac
    taskset_exists = check_executable_exists("taskset", error_if_not_found=False)
    if taskset_exists:
        execute_string += "taskset -c 0"
        if getpass.getuser() != "jovyan": # on mybinder, username is jovyan, and taskset -c 0 is the fastest option.
            # If not on mybinder and taskset exists:
            has_HT_cores = False # Does CPU have hyperthreading cores?
            if platform.processor() != '': # If processor string returns null, then assume CPU does not support hyperthreading.
                                           # This will yield correct behavior on ARM (e.g., cell phone) CPUs.
                has_HT_cores=True
            if has_HT_cores == True:
                # NOTE: You will observe a speed-up by using only *PHYSICAL* (as opposed to logical/hyperthreading) cores:
                N_cores_to_use = int(multiprocessing.cpu_count()/2) # To account for hyperthreading cores
            else:
                N_cores_to_use = int(multiprocessing.cpu_count()) # Use all cores if none are hyperthreading cores.
                                                                  # This will happen on ARM (e.g., cellphone) CPUs
            for i in range(N_cores_to_use-1):
                execute_string += ","+str(i+1)
        execute_string += " "
    execute_string += execute_prefix+executable+" "+executable_output_arguments

    # Step 3: Execute the desired executable
    Execute_input_string(execute_string, file_to_redirect_stdout,verbose)

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

# https://stackoverflow.com/questions/1274405/how-to-create-new-folder
def mkdir(newpath):
    if not os.path.exists(os.path.join(newpath)):
        os.makedirs(os.path.join(newpath))

def output_Jupyter_notebook_to_LaTeXed_PDF(notebookname,location_of_template_file=os.path.join("."),verbose=True):
    Execute_input_string(r"jupyter nbconvert --to latex --template-file "
                         +os.path.join(location_of_template_file,"latex_nrpy_style.tplx")
                         +r" --log-level='WARN' "+notebookname+".ipynb",verbose=False)
    for _i in range(3):  # _i is an unused variable.
        Execute_input_string(r"pdflatex -interaction=batchmode "+notebookname+".tex",verbose=False)
    delete_existing_files(notebookname+".out "+notebookname+".aux "+notebookname+".log")
    if verbose:
        import textwrap
        wrapper = textwrap.TextWrapper(initial_indent="",subsequent_indent="    ",width=75)
        print(wrapper.fill("Created "+notebookname+".tex, and compiled LaTeX file to PDF file "+notebookname+".pdf"))
