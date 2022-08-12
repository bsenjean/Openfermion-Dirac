from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.transforms import jordan_wigner
from openfermion.linalg import eigenspectrum
import numpy as np
import os

# Set molecule parameters.
basis = 'cc-pVDZ'
R = 1.5
angle = 60
angle_rad = angle*np.pi/180.0
charge = 0
data_directory=os.getcwd()

delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
geometry = [('H', (0., R*np.sin(angle_rad), R*np.cos(angle_rad))), ('H', (0., R*np.sin(angle_rad), -R*np.cos(angle_rad))), 
            ('H', (0., -R*np.sin(angle_rad), -R*np.cos(angle_rad))), ('H', (0., -R*np.sin(angle_rad), R*np.cos(angle_rad)))]

print()
print('#'*40)
print('NONREL Dirac calculation')
print('#'*40)
print()
run_ccsd = True
point_nucleus = True
description = 'R' + str(R) + '_T' + str(angle) +  '_ccsd'

molecule = MolecularData_Dirac(geometry=geometry,
                               basis=basis,
                               charge=charge,
                               description=description,
                               data_directory=data_directory)

molecule = run_dirac(molecule,
                    fcidump=True,
                    point_nucleus=point_nucleus,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    run_ccsd=run_ccsd)

print('Spinorbs = ', molecule.get_integrals_FCIDUMP()[1])
print('size spinorbs : {}'.format(len(molecule.get_integrals_FCIDUMP()[1])))
#molecular_hamiltonian = molecule.get_molecular_hamiltonian()[0]
#qubit_hamiltonian = jordan_wigner(molecular_hamiltonian)
#evs = eigenspectrum(qubit_hamiltonian)
print('Hartree-Fock energy of {} Hartree.'.format(molecule.get_energies()[0]))
print('MP2 energy of {} Hartree.'.format(molecule.get_energies()[1]))
print('CCSD energy of {} Hartree.'.format(molecule.get_energies()[2]))
#print('Solving the Qubit Hamiltonian (Jordan-Wigner): \n {}'.format(evs))
