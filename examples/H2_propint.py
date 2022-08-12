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
delete_MDPROP = True
# a ghost nuclei "foo" is added to lower the symmetry.
geometry = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length))]
propint = "ZDIPLEN"

print()
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

properties = ['MOLGRD','DIPOLE','QUADRUPOLE','EFG','POLARIZABILITY']
#properties = False


molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)

molecule = run_dirac(molecule,
                    fcidump=True,
                    propint=propint,
                    properties=properties,
                    point_nucleus=point_nucleus,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    delete_MDPROP=delete_MDPROP,
                    run_ccsd=run_ccsd)

print("spinorbs = ",molecule.get_integrals_FCIDUMP()[1])
print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule.get_energies()[2]))
print('Dipole moment: {}.\n'.format(molecule.get_elecdipole()))
print('Quadrupole moment: {}.\n'.format(molecule.get_elecquadrupole()))
print('Polarizability: {}.\n'.format(molecule.get_elecpolarizability()))
print('property integrals in the MO basis: {}'.format(molecule.get_propint()))

property_hamiltonian = molecule.get_property_hamiltonian()
qubit_hamiltonian = jordan_wigner(property_hamiltonian)
print(property_hamiltonian)
print(qubit_hamiltonian)
