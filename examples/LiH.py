# Test set up for generating Hamiltonian for H2.

from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.hamiltonians import MolecularData
from openfermion.transforms import jordan_wigner, project_onto_sector, bravyi_kitaev
from openfermion.utils import count_qubits,eigenspectrum
import os
import subprocess

# Set molecule parameters.
basis = 'sto-3g'
bond_length = 2.0
multiplicity = 1
charge = 0
data_directory=os.getcwd()

delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
delete_FCIDUMP = False
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

molecule = run_dirac(molecule,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    delete_FCIDUMP=delete_FCIDUMP,
                    run_scf=run_scf,
                    run_ccsd=run_ccsd)

fermion_hamiltonian = molecule.get_molecular_hamiltonian()[0]
qubit_hamiltonian_dirac = jordan_wigner(fermion_hamiltonian)
evs_dirac = eigenspectrum(qubit_hamiltonian_dirac)
print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs_dirac))
