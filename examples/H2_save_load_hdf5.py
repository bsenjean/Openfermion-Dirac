# Test set up for generating Hamiltonian for H2.

from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.transforms import jordan_wigner, project_onto_sector, bravyi_kitaev
from openfermion.utils import eigenspectrum
from openfermion.ops import InteractionOperator
import os
import subprocess
import ast

# Set molecule parameters.
basis = 'sto-3g'
bond_length = 1.0
multiplicity = 1
charge = 0
data_directory=os.getcwd()
geometry = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length))]
delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
delete_FCIDUMP = False
save = True

print('#'*40)
print('NONREL Dirac calculation')
print('#'*40)
print()
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
      print("No file found, running Dirac calculation")
      molecule = run_dirac(molecule,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    delete_FCIDUMP=delete_FCIDUMP,
                    run_ccsd=run_ccsd,
                    save=save)

if save is False:
      fermion_hamiltonian = molecule.get_molecular_hamiltonian()[0]
      qubit_hamiltonian_dirac = jordan_wigner(fermion_hamiltonian)
      evs_dirac = eigenspectrum(qubit_hamiltonian_dirac)
      print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
      print('MP2 energy of {} Hartree.'.format(molecule.get_energies()[1]))
      print('CCSD energy of {} Hartree.'.format(molecule.get_energies()[2]))
      print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(min(evs_dirac)))

else:
      print("HDF5 file found, loading from file. Name : {}".format(molecule.get_from_file('name')))
      print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_from_file('hf_energy')))
      print('MP2 energy of {} Hartree.'.format(molecule.get_from_file('mp2_energy')))
      print('CCSD energy of {} Hartree.'.format(molecule.get_from_file('ccsd_energy')))
      E_core=molecule.get_from_file('nuclear_repulsion')
      one_body_coeff = molecule.get_from_file('one_body_coefficients')
      two_body_coeff = molecule.get_from_file('two_body_coefficients')
      fermionic_ham_print = molecule.get_from_file('print_molecular_hamiltonian')
      fermionic_ham = InteractionOperator(float(E_core),one_body_coeff,two_body_coeff)
      qubit_ham = jordan_wigner(fermionic_ham)
      evs = min(eigenspectrum(qubit_ham))
      print("Ground-state of the qubit Hamiltonian is {}.".format(evs))
