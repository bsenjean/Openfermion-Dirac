from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.transforms import jordan_wigner
from openfermion.utils import eigenspectrum
from openfermion.ops import InteractionOperator
import os
import sys

# Set molecule parameters.
basis = 'STO-3G'
bond_length = 2.0
multiplicity = 1
charge = 0
data_directory=os.getcwd()
geometry = [('Li', (0., 0., 0.)), ('H', (0., 0., bond_length))]

print()
print('#'*40)
print('NONREL CCSD loading')
print('#'*40)
print()
description = 'R' + str(bond_length) + '_ccsd'

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               multiplicity=multiplicity,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)


if os.path.exists("{}/{}.hdf5".format(data_directory,molecule.name)) is False:
      print("No file found. You should first run a calculation with save=True (see LiH_save.py)")
      sys.exit(0)

print("HDF5 file found, loading from file. Name : {}".format(molecule.get_from_file('name')))
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
