# As documented in the NRPy+ tutorial module
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

import io, os, shlex, subprocess, sys, time, multiprocessing

# check_executable_exists(): Check to see whether an executable exists. 
#                            Error out or return False if not exists;
#                            return True if executable exists in PATH.
def check_executable_exists(exec_name,error_if_not_found=True):
    cmd = "where" if os.name == "nt" else "which"
    try: 
        subprocess.check_output([cmd, exec_name])
    except:
        if error_out==True:
            print("Sorry, cannot execute the command: "+exec_name)
            sys.exit(1)
        else:
            return False
    return True

# C_compile(): Write a function to compile the Main C code into an executable file
def C_compile(main_C_output_path, main_C_output_file):
    print("Compiling executable...")
    # Step 1: Check for gcc compiler
    check_executable_exists("gcc")

    # Step 2: Delete existing version of executable
    delete_string = ""
    if os.name == "nt":
        main_C_output_file += ".exe"
    delete_existing_files(main_C_output_file)
    
    # Step 3: Construct the script to run the compilation
    compile_string = "gcc -Ofast -fopenmp -march=native "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"
    Execute_input_string(compile_string, os.devnull)
    # Check if -fopenmp exists, if not, run a slower mode within the compiler
    if not os.path.isfile(main_C_output_file):
        # Step 3.A: Revert to maximally compatible gcc compile option
        print("Optimized compilation failed. Moving to GCC Compatibility (slower) mode:")
        compile_string = "gcc -O2 "+str(main_C_output_path)+" -o "+str(main_C_output_file)+" -lm"
        Execute_input_string(compile_string, os.devnull)
    # Step 3.B: If there are still missing components within the compiler, say compilation failed
    if not os.path.isfile(main_C_output_file):
        print("Sorry, compilation failed")
        sys.exit(1)
    print("Finished compilation.")

# Execute(): Execute generated executable file, using taskset 
#            if available. Calls Execute_input_string() to
#            redirect output from stdout & stderr to desired
#            destinations.
def Execute(executable, executable_output_arguments = "", file_to_redirect_stdout = os.devnull):
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
        execute_string += "cmd /c "+executable.replace("_Playground", "")
    
    taskset_exists = check_executable_exists("taskset",error_if_not_found=False)
    if taskset_exists == True:
        execute_string += "taskset -c 0"
        N_physical_cores = int(multiprocessing.cpu_count()/2) # To account for hyperthreading
        for i in range(N_physical_cores-1):
            execute_string += ","+str(i+1)
        execute_string += " "
    execute_string += os.path.join(".", executable)+" "+executable_output_arguments

    # Step 3: Execute the desired executable
    Execute_input_string(execute_string, file_to_redirect_stdout)

# Execute_input_string(): Executes an input string and redirects 
#            output from stdout & stderr to desired destinations.
def Execute_input_string(input_string, file_to_redirect_stdout):
    print("Executing `"+input_string+"`...")
    start = time.time()
    # https://docs.python.org/3/library/subprocess.html
    args = shlex.split(input_string)
    # https://stackoverflow.com/questions/18421757/live-output-from-subprocess-command
    filename = "tmp.txt"
    with io.open(filename, 'wb') as writer, io.open(filename, 'rb', 1) as reader, io.open(file_to_redirect_stdout, 'w') as rdirect:
        process = subprocess.Popen(args, stdout=rdirect, stderr=writer)
        while process.poll() is None:
            sys.stdout.write(reader.read())
            time.sleep(0.2)
        # Read the remaining
        sys.stdout.write(reader.read())
    delete_existing_files(filename)
# Old Python 2 version (broken for stdout redirect): if (sys.version_info < (3, 0)):
#     else:
#         with open(os.devnull, 'w') as f:  # replace 'w' with 'wb' for Python 3<-- Doesn't work!
#             process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#             for c in iter(lambda: process.stderr.read(1), ''):  # replace '' with b'' for Python 3<-- Doesn't work!
#                 sys.stdout.write(c)
#                 f.write(c)
#             for c in iter(lambda: process.stdout.read(1), ''):  # replace '' with b'' for Python 3<-- Doesn't work!
#                 sys.stdout.write(f.read())
    end = time.time()
    print("Finished executing in "+str(end-start)+" seconds.")

# delete_existing_files(file_or_wildcard): 
#          Runs del file_or_wildcard in Windows, or
#                rm file_or_wildcard in Linux/MacOS
def delete_existing_files(file_or_wildcard):
    delete_string = ""
    if os.name == "nt":
        delete_string += "del "+file_or_wildcard
    else:
        delete_string += "rm -f "+file_or_wildcard
    os.system(delete_string)