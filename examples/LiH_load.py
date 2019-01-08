# Test set up for generating Hamiltonian for H2.

from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.hamiltonians import MolecularData
from openfermion.transforms import jordan_wigner, project_onto_sector, bravyi_kitaev
from openfermion.utils import count_qubits,eigenspectrum
from openfermion.ops import InteractionOperator
import os
import subprocess
import sys

# Set molecule parameters.
basis = 'sto-3g'
bond_length = 2.0
multiplicity = 1
charge = 0
data_directory=os.getcwd()
geometry = [('Li', (0., 0., 0.)), ('H', (0., 0., bond_length))]

print()
print('#'*40)
print('NONREL Dirac calculation')
print('#'*40)
print()
run_scf = 1
if run_scf==1:
 description = 'R' + str(bond_length) + '_scf_dirac'
run_ccsd = 1
if run_ccsd==1:
 description = 'R' + str(bond_length) + '_ccsd_dirac'

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               multiplicity=multiplicity,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)


if os.path.exists("{}/{}.hdf5".format(data_directory,molecule.name)) is False:
      print("No file found. You should first run a calculation with save=True (see LiH_save.py")
      sys.exit(0)

print("HDF5 file found, loading from file. Name : {}", format(molecule.get_from_file('name')))
print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_from_file('hf_energy')))
print('MP2 energy of {} Hartree.'.format(molecule.get_from_file('mp2_energy')))
print('CCSD energy of {} Hartree.'.format(molecule.get_from_file('ccsd_energy')))
E_core=molecule.get_from_file('nuclear_repulsion')
one_body_coeff = molecule.get_from_file('one_body_coefficients')
two_body_coeff = molecule.get_from_file('two_body_coefficients')
fermionic_ham_print = molecule.get_from_file('print_molecular_hamiltonian')
print('Core energy : {} Hartree.'.format(E_core))
#print('Molecular Hamiltonian : {}'.format(fermionic_ham_print))

""" If one actually wants to use it to construct the qubit Hamiltonian, do :
fermionic_ham = InteractionOperator(float(E_core),one_body_coeff,two_body_coeff)
qubit_ham = jordan_wigner(fermionic_ham)
evs = min(eigenspectrum(qubit_ham))
"""
