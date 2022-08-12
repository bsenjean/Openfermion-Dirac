from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.transforms import jordan_wigner, get_fermion_operator, symmetry_conserving_bravyi_kitaev
from openfermion.linalg import eigenspectrum
import os

# Set molecule parameters.
basis = 'STO-3G'
bond_length = 1
charge = 1
data_directory=os.getcwd()

delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
geometry = [('He', (0., 0., 0.)), ('H', (0., 0., bond_length))]

print()
print('#'*40)
print('REL CCSD Dirac calculation')
print('#'*40)
print()
run_ccsd = True
relativistic = True
description = 'R' + str(bond_length) + '_ccsd'

# This manual option is not necessary anymore and has been fixed since then.
#manual_option="""**RELCCSD
#*CCENER
#.NOSDT"""

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               charge=charge,
                               description=description,
                               relativistic=relativistic,
                               data_directory=data_directory)

molecule = run_dirac(molecule,
                    fcidump=True,
                    relativistic=relativistic,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    run_ccsd=run_ccsd)

molecular_hamiltonian = molecule.get_molecular_hamiltonian()[0]
qubit_hamiltonian = jordan_wigner(molecular_hamiltonian)
evs = eigenspectrum(qubit_hamiltonian)
print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs))

#symmetry_conserving_bravyi_kitaev(fermionicoperator,number_of_active_orbs,number_of_active_elec)
number_orbs = len(molecule.get_integrals_FCIDUMP()[1])
fermion_operator = get_fermion_operator(molecular_hamiltonian)
qubit_hamiltonian = symmetry_conserving_bravyi_kitaev(fermion_operator,number_orbs,2)
evs = eigenspectrum(qubit_hamiltonian)
print('Solving the Qubit Hamiltonian (Bravyi-Kitaev) Symmetry conserving: \n {}'.format(evs))

print()
print('#'*40)
print('REL FCI Dirac calculation')
print('#'*40)
print()
run_fci = True
relativistic = True
description = 'R' + str(bond_length) + '_fci'

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               charge=charge,
                               description=description,
                               relativistic=relativistic,
                               data_directory=data_directory)

molecule = run_dirac(molecule,
                    relativistic=relativistic,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    run_fci=run_fci)

print('FCI energy of {} Hartree.'.format(molecule.get_energies()[3]))
