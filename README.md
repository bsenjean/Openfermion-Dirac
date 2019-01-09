# Dirac interfaced with Openfermion

# Requirements

- Dirac 
- Openfermion
- all requirements for the aforementioned programs.

# Installation

INSTALL DIRAC: http://diracprogram.org/doku.php

INSTALL OPENFERMION : https://github.com/quantumlib/OpenFermion

INSTALL Openfermion-Dirac interface:

Clone the complete repository:
```
$ git clone https://github.com/bsenjean/Openfermion-Dirac.git
```

Build the code:
```
$ cd /path/to/Openfermion-Dirac/
$ pip install -e .
$ cd utils/
$ gfortran dirac_openfermion_mointegral_export.F90 -o dirac_openfermion_mointegral_export.x
```

In /path/to/Openfermion-Dirac/openfermion_dirac/_run_dirac.py change the following :
- /path/to/dirac/build/pam to your own path to pam (pam is the run_script of the Dirac program)
- /path/to/Openfermion-Dirac/utils/dirac_openfermion_mointegral_export.x
  that constructs the FCIDUMP integral file from the one_body_integrals (MRCONEE) and the two_body_integrals (MDCINT)

# Examples

Different examples are furnished in the examples/ repository, as well as a tutorial.

# Known limits and problems

From the side of Dirac: 
- Charged systems : sometimes (like HeH+) one cannot run the CCSD calculation with Dirac directly.
                  This is most certainely due to the fact that the program does not know how to
                  distribute the remaining electrons to make the reference HF state...
                  It could be fixed by explicitely saying where to distribute the electrons in
                  the input. This is yet to be done.
                  For Be+, the energy given by Dirac is the one of Be2+. This is again because
                  CCSD does not know what to do with this open shell system, and just considers
                  of 2 electrons instead of 3.

- For Beryllium : CCSD calculation does not work.

From the side of Openfermion:
- Fock space :  The second-quantized Hamiltonian, and so the qubit Hamiltonian,
              operates in the Fock space. Hence, the ground-state is not
              always the expected one. This is especially true for ionic species, 
              such as HeH+ for instance. One should add an option to constrain charge 
              (and spin) conservation to get only the states we are inerested of. 
              This is yet to be done.
