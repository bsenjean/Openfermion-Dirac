from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.transforms import jordan_wigner, bravyi_kitaev, symmetry_conserving_bravyi_kitaev, get_fermion_operator
from openfermion.utils import eigenspectrum
import os

# Set molecule parameters.
basis = 'STO-3G'
charge = 1
data_directory=os.getcwd()

delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
fcidump = True
geometry = [('Be', (0., 0., 0.))]

print()
print('#'*40)
print('NONREL Dirac calculation')
print('#'*40)
print()
run_ccsd = True
point_nucleus = True
description = 'ccsd'

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)

molecule = run_dirac(molecule,
                    fcidump=fcidump,
                    point_nucleus=point_nucleus,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    run_ccsd=run_ccsd)

molecular_hamiltonian = molecule.get_molecular_hamiltonian()[0]
number_orbs = len(molecule.get_integrals_FCIDUMP()[1])
print('size spinorbs : {}'.format(number_orbs))
qubit_hamiltonian = jordan_wigner(molecular_hamiltonian)
evs = eigenspectrum(qubit_hamiltonian)
print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
print('MP2 energy of {} Hartree. (WRONG)'.format(molecule.get_energies()[1]))
print('CCSD energy of {} Hartree. (WRONG)'.format(molecule.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs))

fermion_operator = get_fermion_operator(molecular_hamiltonian)
qubit_hamiltonian = bravyi_kitaev(fermion_operator)
evs = eigenspectrum(qubit_hamiltonian)
print('Solving the Qubit Hamiltonian (Bravyi-Kitaev): \n {}'.format(evs))
#symmetry_conserving_bravyi_kitaev(fermionicoperator,number_of_active_orbs,number_of_active_elec)
qubit_hamiltonian = symmetry_conserving_bravyi_kitaev(fermion_operator,number_orbs,3)
evs = eigenspectrum(qubit_hamiltonian)
print('Solving the Qubit Hamiltonian (Bravyi-Kitaev) Symmetry conserving: \n {}'.format(evs))

molecule = run_dirac(molecule,
                    point_nucleus=point_nucleus,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    run_fci=True)

print('FCI energy of {} Hartree. (CORRECT)'.format(molecule.get_energies()[3]))
