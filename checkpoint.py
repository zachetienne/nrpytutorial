from outputC import Cfunction
from loop import loop, simple_loop

def indent(loopstring, ntabs=1):
	ind = "\n" + "    " * ntabs
	returnstring = ""
	for line in iter(loopstring.splitlines()):
		returnstring += ind + line
	returnstring += "\n"

	return returnstring

# def initial_data_from_file():

desc = "Set up the initial data using data from a file, for checkpointing purposes"
name = "initial_data_from_file"
params = "const paramstruct *restrict params, REAL * restrict in_gfs, REAL t_initial"

outer = loop("n_gf", "0", "NUM_EVOL_GFS", "1", "#pragma omp parallel for")
preloop = indent("""char in_filename[50];
sprintf(in_filename, "checkpoint_data/gf_%d_tfinal_%.2f.txt", n_gf, t_initial);
FILE *in;
in = fopen(in_filename, "r+");""")

postloop = indent("fclose(in);")

loopopts = "AllPoints,DisableOpenMP"
loopbody = 'fscanf(in, "%lf", &in_gfs[IDX4S(n_gf, i0, i1, i2)]);'
inner = indent(simple_loop(loopopts, loopbody), 0)

body = indent(outer[0] + preloop + inner + postloop + outer[1])

prototype, func = Cfunction(desc=desc, name=name, params=params, body=body)
print(func)	

