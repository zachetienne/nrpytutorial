# This code calculates the invariant scalars I, J, J1, J2, J3, and J4. It does so by following the papers
# arXiv:gr-qc/0407013 and arXiv:0704.1756, and the example set by the Kranc-generated ETK thorn which
# can be found at https://bitbucket.org/einsteintoolkit/einsteinanalysis/src. While WeylScalarInvariants_Cartesian()
# does not depend on WeylScalars_Cartesian(), when implementing C code generated from the outputs of these two
# functions, the psis from WeylScalars_Cartesian() must be calculated before the invariants generated here.

# As documented in Step 6 of the NRPy+ tutorial module
#             Tutorial-WeylScalarsInvariants-Cartesian

# Step P1: import all needed modules from NRPy+:
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends
import grid as gri                # NRPy+: Functions having to do with numerical grids

# Step P2: Initialize WeylScalar Invariants parameters
# None at this time.

def WeylScalarInvariants_Cartesian():
    # Step 1: Declare Weyl scalar psi's as gridfunctions, as they will be read from memory.
    psi4r,psi4i,psi3r,psi3i,psi2r,psi2i,psi1r,psi1i,psi0r,psi0i = gri.register_gridfunctions("AUX",["psi4r","psi4i",
                                                                                                    "psi3r","psi3i",
                                                                                                    "psi2r","psi2i",
                                                                                                    "psi1r","psi1i",
                                                                                                    "psi0r","psi0i"])

    # We will first set the complex psis as the appropriate combinations of their real and imaginary parts.
    # This will allow us to take advantage of SymPy's functions for handling complex numbers.
    psi4 = psi4r + sp.I * psi4i
    psi3 = psi3r + sp.I * psi3i
    psi2 = psi2r + sp.I * psi2i
    psi1 = psi1r + sp.I * psi1i
    psi0 = psi0r + sp.I * psi0i

    # We declare the variants as global to access them from other codes. We will also declare them as gridfunctions.
    global curvIr,curvIi,curvJr,curvJi,J1curv,J2curv,J3curv,J4curv
    curvIr, curvIi, curvJr, curvJi, J1curv, J2curv, J3curv, J4curv = gri.register_gridfunctions("AUX",["curvIr","curvIi",
                                                                                                       "curvJr","curvJi",
                                                                                                       "J1curv","J2curv",
                                                                                                       "J3curv","J4curv"])

    # The equations for the real and complex parts of I and J, from arXiv:gr-qc/0407013, equations (2.2a) and (2.2b):
    # I &= 3 \psi_2^2 - 4 \psi_1 \psi_3 + \psi_4 \psi_0 \\
    # J &= det |\psi_4 & \psi_3 & \psi_2|
    #          |\psi_3 & \psi_2 & \psi_1|
    #          |\psi_2 & \psi_1 & \psi_0|

    curvIr = sp.re(3*psi2*psi2 - 4*psi1*psi3 + psi4*psi0)
    curvIi = sp.im(3*psi2*psi2 - 4*psi1*psi3 + psi4*psi0)
    curvJr = sp.re(psi4 * (psi2*psi0 - psi1*psi1) -
                   psi3 * (psi3*psi0 - psi1*psi2) +
                   psi2 * (psi3*psi1 - psi2*psi2) )
    curvJi = sp.im(psi4 * (psi2*psi0 - psi1*psi1) -
                   psi3 * (psi3*psi0 - psi1*psi2) +
                   psi2 * (psi3*psi1 - psi2*psi2) )

    # Next, the real scalars J_1, J_2, J_3, and J_4 from equations B5-B8 of arXiv:0704.1756.
    # These equations are based directly on those used in the Mathematica notebook that generates WeylScal4
    # (available at https://bitbucket.org/einsteintoolkit/einsteinanalysis/src), modified so that Python can
    # interpret them. Those equations were generated in turn using xTensor from equations B5-B8.
    J1curv =-16*(3*psi2i**2-3*psi2r**2-4*psi1i*psi3i+4*psi1r*psi3r+psi0i*psi4i-psi0r*psi4r)

    J2curv = 96*(-3*psi2i**2*psi2r+psi2r**3+2*psi1r*psi2i*psi3i+2*psi1i*psi2r*psi3i-psi0r*psi3i**2+
                2*psi1i*psi2i*psi3r-2*psi1r*psi2r*psi3r-2*psi0i*psi3i*psi3r+psi0r*psi3r**2-
                2*psi1i*psi1r*psi4i+psi0r*psi2i*psi4i+psi0i*psi2r*psi4i-psi1i**2*psi4r+psi1r**2*psi4r+
                psi0i*psi2i*psi4r-psi0r*psi2r*psi4r)

    J3curv = 64*(9*psi2i**4-54*psi2i**2*psi2r**2+9*psi2r**4-24*psi1i*psi2i**2*psi3i+48*psi1r*psi2i*psi2r*psi3i+
                24*psi1i*psi2r**2*psi3i+16*psi1i**2*psi3i**2-16*psi1r**2*psi3i**2+
                24*psi1r*psi2i**2*psi3r+48*psi1i*psi2i*psi2r*psi3r-24*psi1r*psi2r**2*psi3r-64*psi1i*psi1r*psi3i*psi3r-
                16*psi1i**2*psi3r**2+16*psi1r**2*psi3r**2+6*psi0i*psi2i**2*psi4i-12*psi0r*psi2i*psi2r*psi4i-
                6*psi0i*psi2r**2*psi4i-8*psi0i*psi1i*psi3i*psi4i+8*psi0r*psi1r*psi3i*psi4i+8*psi0r*psi1i*psi3r*psi4i+
                8*psi0i*psi1r*psi3r*psi4i+psi0i**2*psi4i**2-psi0r**2*psi4i**2-6*psi0r*psi2i**2*psi4r-
                12*psi0i*psi2i*psi2r*psi4r+6*psi0r*psi2r**2*psi4r+8*psi0r*psi1i*psi3i*psi4r+8*psi0i*psi1r*psi3i*psi4r+
                8*psi0i*psi1i*psi3r*psi4r-8*psi0r*psi1r*psi3r*psi4r-4*psi0i*psi0r*psi4i*psi4r-psi0i**2*psi4r**2+
                psi0r**2*psi4r**2)

    J4curv = -640*(-15*psi2i**4*psi2r+30*psi2i**2*psi2r**3-3*psi2r**5+10*psi1r*psi2i**3*psi3i+
                  30*psi1i*psi2i**2*psi2r*psi3i-30*psi1r*psi2i*psi2r**2*psi3i-10*psi1i*psi2r**3*psi3i-
                  16*psi1i*psi1r*psi2i*psi3i**2-3*psi0r*psi2i**2*psi3i**2-8*psi1i**2*psi2r*psi3i**2+
                  8*psi1r**2*psi2r*psi3i**2-6*psi0i*psi2i*psi2r*psi3i**2+3*psi0r*psi2r**2*psi3i**2+
                  4*psi0r*psi1i*psi3i**3+4*psi0i*psi1r*psi3i**3+10*psi1i*psi2i**3*psi3r-
                  30*psi1r*psi2i**2*psi2r*psi3r-30*psi1i*psi2i*psi2r**2*psi3r+10*psi1r*psi2r**3*psi3r-
                  16*psi1i**2*psi2i*psi3i*psi3r+16*psi1r**2*psi2i*psi3i*psi3r-6*psi0i*psi2i**2*psi3i*psi3r+
                  32*psi1i*psi1r*psi2r*psi3i*psi3r+12*psi0r*psi2i*psi2r*psi3i*psi3r+6*psi0i*psi2r**2*psi3i*psi3r+
                  12*psi0i*psi1i*psi3i**2*psi3r-12*psi0r*psi1r*psi3i**2*psi3r+16*psi1i*psi1r*psi2i*psi3r**2+
                  3*psi0r*psi2i**2*psi3r**2+8*psi1i**2*psi2r*psi3r**2-8*psi1r**2*psi2r*psi3r**2+
                  6*psi0i*psi2i*psi2r*psi3r**2-3*psi0r*psi2r**2*psi3r**2-12*psi0r*psi1i*psi3i*psi3r**2-
                  12*psi0i*psi1r*psi3i*psi3r**2-4*psi0i*psi1i*psi3r**3+4*psi0r*psi1r*psi3r**3-
                  6*psi1i*psi1r*psi2i**2*psi4i+2*psi0r*psi2i**3*psi4i-6*psi1i**2*psi2i*psi2r*psi4i+
                  6*psi1r**2*psi2i*psi2r*psi4i+6*psi0i*psi2i**2*psi2r*psi4i+6*psi1i*psi1r*psi2r**2*psi4i-
                  6*psi0r*psi2i*psi2r**2*psi4i-2*psi0i*psi2r**3*psi4i+12*psi1i**2*psi1r*psi3i*psi4i-
                  4*psi1r**3*psi3i*psi4i-2*psi0r*psi1i*psi2i*psi3i*psi4i-2*psi0i*psi1r*psi2i*psi3i*psi4i-
                  2*psi0i*psi1i*psi2r*psi3i*psi4i+2*psi0r*psi1r*psi2r*psi3i*psi4i-2*psi0i*psi0r*psi3i**2*psi4i+
                  4*psi1i**3*psi3r*psi4i-12*psi1i*psi1r**2*psi3r*psi4i-2*psi0i*psi1i*psi2i*psi3r*psi4i+
                  2*psi0r*psi1r*psi2i*psi3r*psi4i+2*psi0r*psi1i*psi2r*psi3r*psi4i+2*psi0i*psi1r*psi2r*psi3r*psi4i-
                  2*psi0i**2*psi3i*psi3r*psi4i+2*psi0r**2*psi3i*psi3r*psi4i+2*psi0i*psi0r*psi3r**2*psi4i-
                  psi0r*psi1i**2*psi4i**2-2*psi0i*psi1i*psi1r*psi4i**2+psi0r*psi1r**2*psi4i**2+
                  2*psi0i*psi0r*psi2i*psi4i**2+psi0i**2*psi2r*psi4i**2-psi0r**2*psi2r*psi4i**2-3*psi1i**2*psi2i**2*psi4r+
                  3*psi1r**2*psi2i**2*psi4r+2*psi0i*psi2i**3*psi4r+12*psi1i*psi1r*psi2i*psi2r*psi4r-
                  6*psi0r*psi2i**2*psi2r*psi4r+3*psi1i**2*psi2r**2*psi4r-3*psi1r**2*psi2r**2*psi4r-
                  6*psi0i*psi2i*psi2r**2*psi4r+2*psi0r*psi2r**3*psi4r+4*psi1i**3*psi3i*psi4r-12*psi1i*psi1r**2*psi3i*psi4r-
                  2*psi0i*psi1i*psi2i*psi3i*psi4r+2*psi0r*psi1r*psi2i*psi3i*psi4r+2*psi0r*psi1i*psi2r*psi3i*psi4r+
                  2*psi0i*psi1r*psi2r*psi3i*psi4r-psi0i**2*psi3i**2*psi4r+psi0r**2*psi3i**2*psi4r-
                  12*psi1i**2*psi1r*psi3r*psi4r+4*psi1r**3*psi3r*psi4r+2*psi0r*psi1i*psi2i*psi3r*psi4r+
                  2*psi0i*psi1r*psi2i*psi3r*psi4r+2*psi0i*psi1i*psi2r*psi3r*psi4r-2*psi0r*psi1r*psi2r*psi3r*psi4r+
                  4*psi0i*psi0r*psi3i*psi3r*psi4r+psi0i**2*psi3r**2*psi4r-psi0r**2*psi3r**2*psi4r-
                  2*psi0i*psi1i**2*psi4i*psi4r+4*psi0r*psi1i*psi1r*psi4i*psi4r+2*psi0i*psi1r**2*psi4i*psi4r+
                  2*psi0i**2*psi2i*psi4i*psi4r-2*psi0r**2*psi2i*psi4i*psi4r-4*psi0i*psi0r*psi2r*psi4i*psi4r+
                  psi0r*psi1i**2*psi4r**2+2*psi0i*psi1i*psi1r*psi4r**2-psi0r*psi1r**2*psi4r**2-
                  2*psi0i*psi0r*psi2i*psi4r**2-psi0i**2*psi2r*psi4r**2+psi0r**2*psi2r*psi4r**2)
