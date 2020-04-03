from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.transforms import jordan_wigner
from openfermion.utils import eigenspectrum
import os

# Set molecule parameters.
basis = 'STO-3G'
bond_length = 0.5
charge = 1
data_directory=os.getcwd()

delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
fcidump = True
geometry = [('He', (0., 0., 0.)), ('H', (0., 0., bond_length))]

print()
print('#'*40)
print('NONREL Dirac calculation, cc-pVDZ basis for H')
print('#'*40)
print()

run_ccsd = True
if run_ccsd:
 description = 'R' + str(bond_length) + '_ccsd'
else:
 description = 'R' + str(bond_length) + '_scf'

manual_option = "**RELCCSD\n*CCENER\n.NOSDT"

basis="special" # necessary keyword to specify that we are using a special basis
special_basis=["STO-3G","H BASIS cc-pVDZ"] # list of two string : [default basis, special basis for a given atomic species]
point_nucleus = True
molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               special_basis=special_basis,
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
                    run_ccsd=run_ccsd,
                    manual_option=manual_option)

molecular_hamiltonian = molecule.get_molecular_hamiltonian()[0]
qubit_hamiltonian = jordan_wigner(molecular_hamiltonian)
evs = eigenspectrum(qubit_hamiltonian)
print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule.get_energies()[2]))
print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs))
