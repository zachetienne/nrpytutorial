# As documented in the NRPy+ tutorial module
#   Tutorial-SpinWeighted_Spherical_Harmonics.ipynb,
#   this module will produce the required C code for
#   computing spin weighted spherical hamronics for s = -2,
#   for l in [0,maximum_l] and m in [-l, l], using
#   the Goldberg formula

# Authors: Brandon Clark
#          Zachariah B. Etienne
#          zachetie **at** gmail **dot* com

# The Goldberg formula can be found at the following citation:
#    https://aip.scitation.org/doi/10.1063/1.1705135
# Wikipedia also has an article on Spin-Weighted Spherical Hamronics:
#   https://en.wikipedia.org/w/index.php?title=Spin-weighted_spherical_harmonics&oldid=853425244

# Step 1: Initialize needed Python/NRPy+ modules
from outputC import outputC       # NRPy+: Core C code output module
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import os                         # Python built-in: Multiplatform operating system functions

def SpinWeight_minus2_SphHarmonics(maximum_l=8,filename=os.path.join("SpinWeight_minus2_SphHarmonics","SpinWeight_minus2_SphHarmonics.h")):
    # Step 2: Defining the Goldberg function

    # Step 2.a: Declare SymPy symbols:
    th, ph = sp.symbols('th ph',real=True)

    # Step 2.b: Define the Goldberg formula for spin-weighted spherical harmonics
    #           (https://aip.scitation.org/doi/10.1063/1.1705135);
    #           referenced & described in Wikipedia Spin-weighted spherical harmonics article:
    #           https://en.wikipedia.org/w/index.php?title=Spin-weighted_spherical_harmonics&oldid=853425244
    def Y(s, l, m, th, ph, GenerateMathematicaCode=False):
        Sum = 0
        for r in range(l-s + 1):
            if GenerateMathematicaCode == True:
                # Mathematica needs expression to be in terms of cotangent, so that code validation below
                #    yields identity with existing Mathematica notebook on spin-weighted spherical harmonics.
                Sum +=  sp.binomial(l-s, r)*sp.binomial(l+s, r+s-m)*(-1)**(l-r-s)*sp.exp(sp.I*m*ph)*sp.cot(th/2)**(2*r+s-m)
            else:
                # SymPy C code generation cannot handle the cotangent function, so define cot(th/2) as 1/tan(th/2):
                Sum +=  sp.binomial(l-s, r)*sp.binomial(l+s, r+s-m)*(-1)**(l-r-s)*sp.exp(sp.I*m*ph)/sp.tan(th/2)**(2*r+s-m)

        return (-1)**m*sp.simplify(sp.sqrt(sp.factorial(l+m)*sp.factorial(l-m)*(2*l+1)/(4*sp.pi*sp.factorial(l+s)*sp.factorial(l-s)))*sp.sin(th/2)**(2*l)*Sum)

    # Step 3: (DISABLED FOR NOW; PASSES TEST).
    #         Code Validation against Mathematica notebook:
    #         https://demonstrations.wolfram.com/versions/source.jsp?id=SpinWeightedSphericalHarmonics&version=0012

    # # For the l=0 case m=0, otherwise there is a divide-by-zero in the Y() function above.
    # print("FullSimplify[Y[-2, 0, 0, th, ph]-"+str(sp.mathematica_code(sp.simplify(Y(-2, 0, 0, th, ph,GenerateMathematicaCode=True))))+"] \n") # Agrees with Mathematica notebook for l = 0

    # # Check the other cases
    # for l in range(1,9): # Agrees with Mathematica notebook for  l = 1, 2, 4, 5, 6, 7, 8;
    #     print("FullSimplify[Y[-2, "+str(l)+", m, th, ph]-("+
    #           str(sp.mathematica_code(sp.simplify(Y(-2, l, m, th, ph, GenerateMathematicaCode=True)))).replace("binomial","Binomial").replace("factorial","Factorial")+")] \n")


    # Step 4: Generating C Code function for computing
    #         s=-2 spin-weighted spherical harmonics,
    #         using NRPy+'s outputC() function.

    outCparams = "preindent=3,outCfileaccess=a,outCverbose=False,includebraces=True"

    with open(filename, "w") as file:
        file.write("""
void SpinWeight_minus2_SphHarmonics(const int l, const int m, const REAL th, const REAL ph,
                                   REAL *reYlmswm2_l_m, REAL *imYlmswm2_l_m) {
if(l<0 || l>"""+str(maximum_l)+""" || m<-l || m>+l) {
    printf("ERROR: SpinWeight_minus2_SphHarmonics handles only l=[0,"""+str(maximum_l)+"""] and only m=[-l,+l] is defined.\\n");
    printf("       You chose l=%d and m=%d, which is out of these bounds.\\n",l,m);
    exit(1);
}\n""")

        file.write("switch(l) {\n")
        for l in range(maximum_l+1): # Output values up to and including l=8.
            file.write("    case "+str(l)+":\n")
            file.write("        switch(m) {\n")
            for m in range(-l,l+1):
                file.write("            case "+str(m)+":\n")
                Y_m2_lm = Y(-2, l, m, th, ph)
                Cstring = outputC([sp.re(Y_m2_lm),sp.im(Y_m2_lm)],["*reYlmswm2_l_m","*imYlmswm2_l_m"],
                                  "returnstring",outCparams)
                file.write(Cstring)
                file.write("                  return;\n")
            file.write("        }  /* End switch(m) */\n")
        file.write("    } /* End switch(l) */\n")
        file.write("} /* End function SpinWeight_minus2_SphHarmonics() */\n")
