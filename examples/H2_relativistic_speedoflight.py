from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.transforms import jordan_wigner
from openfermion.linalg import eigenspectrum
import os

# Set molecule parameters.
basis = 'STO-3G'
bond_length = 1.0
charge = 0
data_directory=os.getcwd()

delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
fcidump = True
geometry = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length))]

speed_of_light = 3.0

print()
print('#'*60)
print('REL Dirac calculation with speed of light = ', speed_of_light, ' a.u.')
print('#'*60)
print()

run_ccsd = True
description = 'R' + str(bond_length) + '_ccsd'
relativistic=True

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               charge=charge,
                               speed_of_light=speed_of_light,
                               description=description,
                               relativistic=relativistic,
                               data_directory=data_directory)

molecule = run_dirac(molecule,
                    fcidump=fcidump,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    run_ccsd=run_ccsd,
                    relativistic=relativistic,
                    speed_of_light=speed_of_light)

molecular_hamiltonian = molecule.get_molecular_hamiltonian()[0]
qubit_hamiltonian = jordan_wigner(molecular_hamiltonian)
evs = eigenspectrum(qubit_hamiltonian)
print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs))
