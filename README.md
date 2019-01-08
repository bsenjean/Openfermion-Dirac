Dirac interfaced with Openfermion

##############################
1) Install                   #
##############################

- Dirac 
- Openfermion
- Psi4
- Openfermion-Psi4
- all requirements for the aforementioned programs.

In Dirac:
$> cd path/to/dirac/src/openfermion_interface 
$> pip install -e .

##############################
2) Use                       #
##############################

Add the path/to/dirac/build/ directory to your bash_profile to use
- pam, the run_script of dirac
- dirac_openfermion_mointegrals_export, created to construct the FCIDUMP 
  from the one_body_integrals (MRCONEE) and the two_body_integrals (MDCINT)

##############################
3) Examples                  #
##############################

Different examples are furnished in src/openfermion_interface/work_dir/

##############################
4) Known limits and problems #
##############################

- From the side of Dirac: 
Charged systems : sometimes (like HeH+) one cannot run the CCSD calculation with Dirac directly.
                  This is most certainely due to the fact that the program does not know how to
                  distribute the remaining electrons to make the reference HF state...
                  It could be fixed by explicitely saying where to distribute the electrons in
                  the input. This is yet to be done.
                  For Be+, the energy given by Dirac is the one of Be2+. This is again because
                  CCSD does not know what to do with this open shell system, and just get rid
                  of 2 electrons instead of one.

CCSD energy : for small system, this energy can be wrong because of a bug. This should be
              fixed in the next release of Dirac, i.e. very soon.

- From the side of Openfermion:
Fock space : Solving the Qubit Hamiltonian gives you all possible states with all possible
                  number of particles and repartitions of them. It contains all the solutions
                  in the Fock space (extended Hilbert space). This is because the second-quantized
                  Hamiltonian operates in the Fock space. Hence, the ground-state is not
                  always the expected one, especially for Be2+ or HeH+ for instance. One should
                  add an option to constrain charge (and spin) conservation to get only the states
                  we are inerested of. This is yet to be done.
