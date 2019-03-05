from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.transforms import jordan_wigner
from openfermion.utils import eigenspectrum
import os

# Set molecule parameters.
basis = 'special'
special_basis = ["STO-3G","foo NOBASIS"]
bond_length = 1.0
multiplicity = 1
charge = 0
data_directory=os.getcwd()

delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
delete_FCIDUMP = False
# a ghost nuclei "foo" is added to lower the symmetry.
geometry = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length)), ('foo', (0.,0., 10.))]
operator = [["ZDIPLEN","COMFACTOR",0.001]]#["XDIPLEN","COMFACTOR",0.0],["YDIPLEN","COMFACTOR",0.0]]
for i in range(len(operator)):
 for j in range(len(operator[i])):
  print(operator[i][j])

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

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               special_basis=special_basis,
                               multiplicity=multiplicity,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)

molecule = run_dirac(molecule,
                    operator=operator,
                    point_nucleus=point_nucleus,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    delete_FCIDUMP=delete_FCIDUMP,
                    run_ccsd=run_ccsd)

print("spinorbs = ",molecule.get_integrals_FCIDUMP()[1])
print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule.get_energies()[2]))
