from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.transforms import jordan_wigner
from openfermion.utils import eigenspectrum
import os

# Set molecule parameters.
basis = 'STO-3G'
bond_length = 0.7414
multiplicity = 1
charge = 0
data_directory=os.getcwd()

delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
fcidump=True
geometry = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length))]

print()
print('#'*40)
print('NONREL Dirac calculation')
print('#'*40)
print()
point_nucleus = True
run_ccsd=False
if run_ccsd:
 description = 'R' + str(bond_length) + '_ccsd'
else:
 description = 'R' + str(bond_length) + '_scf'
save = True

properties = ['MOLGRD','DIPOLE','QUADRUPOLE','EFG','POLARIZABILITY']

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               multiplicity=multiplicity,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)

molecule = run_dirac(molecule,
                    fcidump=fcidump,
                    point_nucleus=point_nucleus,
                    properties=properties,
                    save=save,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT)

print('Hartree-Fock energy of {} Hartree. From the hdf5 file: {}'.format(molecule.get_energies()[0],molecule.get_from_file('hf_energy')))
print('Dipole moment: {}. From the hdf5 file: {}.'.format(molecule.get_elecdipole(),molecule.get_from_file('elec_dipole')))
print('Quadrupole moment: {}. From the hdf5 file: {}'.format(molecule.get_elecquadrupole(),molecule.get_from_file('elec_quadrupole')))
print('Polarizability: {}. From the hdf5 file: {}'.format(molecule.get_elecpolarizability(),molecule.get_from_file('elec_polarizability')))
print("The gradient can be read in the output. I've not define an option to extract it from the ouput yet.")
