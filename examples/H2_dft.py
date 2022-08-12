from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.transforms import jordan_wigner
from openfermion.linalg import eigenspectrum
import os

# Set molecule parameters.
basis = 'STO-3G'
bond_length = 3.0
charge = 0
data_directory=os.getcwd()

delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
geometry = [('H', (0., 0., 0.)), ('H', (0., 0., bond_length))]

print()
print('#'*40)
print('NONREL Dirac calculation')
print('#'*40)
print()
run_dft="LDA"
point_nucleus = True
description = 'R' + str(bond_length) + '_' + run_dft

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)

molecule = run_dirac(molecule,
                    point_nucleus=point_nucleus,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    NONREL=True,
                    get="DFCOEF",
                    run_dft=run_dft)

print('DFT energy of {} Hartree.'.format(molecule.get_energies()[0]))
