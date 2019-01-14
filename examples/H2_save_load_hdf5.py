from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.transforms import jordan_wigner
from openfermion.utils import eigenspectrum
from openfermion.ops import InteractionOperator
import os

# Set molecule parameters.
basis = 'STO-3G'
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
run_ccsd = True
point_nucleus = True
if run_ccsd:
 description = 'R' + str(bond_length) + '_ccsd'
else:
 description = 'R' + str(bond_length) + '_scf'

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               multiplicity=multiplicity,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)

if os.path.exists("{}/{}.hdf5".format(data_directory,molecule.name)) is False:
      print("No file found, running Dirac calculation")
      molecule = run_dirac(molecule,
                    point_nucleus=point_nucleus,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    delete_FCIDUMP=delete_FCIDUMP,
                    run_ccsd=run_ccsd,
                    save=save)

if save is False:
      molecular_hamiltonian = molecule.get_molecular_hamiltonian()[0]
      qubit_hamiltonian = jordan_wigner(molecular_hamiltonian)
      evs = eigenspectrum(qubit_hamiltonian)
      print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
      print('MP2 energy of {} Hartree.'.format(molecule.get_energies()[1]))
      print('CCSD energy of {} Hartree.'.format(molecule.get_energies()[2]))
      print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(min(evs)))

else:
      print("HDF5 file found, loading from file. Name : {}".format(molecule.get_from_file('name')))
      print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_from_file('hf_energy')))
      print('MP2 energy of {} Hartree.'.format(molecule.get_from_file('mp2_energy')))
      print('CCSD energy of {} Hartree.'.format(molecule.get_from_file('ccsd_energy')))
      E_core=molecule.get_from_file('nuclear_repulsion')
      one_body_coeff = molecule.get_from_file('one_body_coefficients')
      two_body_coeff = molecule.get_from_file('two_body_coefficients')
      molecular_ham_print = molecule.get_from_file('print_molecular_hamiltonian')
      molecular_ham = InteractionOperator(float(E_core),one_body_coeff,two_body_coeff)
      qubit_ham = jordan_wigner(molecular_ham)
      evs = min(eigenspectrum(qubit_ham))
      print("Ground-state of the qubit Hamiltonian is {}.".format(evs))
