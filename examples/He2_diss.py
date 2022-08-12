from openfermion_dirac import MolecularData_Dirac, run_dirac
from openfermion.transforms import jordan_wigner
from openfermion.linalg import eigenspectrum
import os

# Set molecule parameters.
basis = 'cc-pVDZ'
bond_length = 1.0
charge = 0
data_directory=os.getcwd()

delete_input = True
delete_xyz = True
delete_output = False
delete_MRCONEE = True
delete_MDCINT = True
run_ccsd = True
point_nucleus = True

bond_length_interval = 0.1
n_points = 50

# Generate molecule at different bond lengths.
hf_energies = []
mp2_energies= []
ccsd_energies = []
fci_energies = []
bond_lengths = []
for point in range(n_points):
    bond_length = 2.8 + bond_length_interval * point
    bond_lengths += [bond_length]
    description = str(round(bond_length,2))
    print(description)
    geometry = [('He', (0., 0., 0.)), ('He', (0., 0., bond_length))]
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
                    run_ccsd=run_ccsd)
    hf_energies += [molecule.get_energies()[0]]
    mp2_energies += [molecule.get_energies()[1]]
    ccsd_energies += [molecule.get_energies()[2]]

    molecule = run_dirac(molecule,
                    point_nucleus=point_nucleus,
                    delete_input=delete_input,
                    delete_xyz=delete_xyz,
                    delete_output=delete_output,
                    delete_MRCONEE=delete_MRCONEE,
                    delete_MDCINT=delete_MDCINT,
                    run_fci=True)
    fci_energies += [molecule.get_energies()[3]]
with open("He2_diss.dat","w") as f:
  f.write('%15s %15s %15s %15s %15s\n' % ("bond_length","HF","MP2","CCSD","FCI"))
  for i in range(len(bond_lengths)):
      f.write('%15.6f %15.8f %15.8f %15.8f %15.8f\n' % (float(bond_lengths[i]),float(hf_energies[i]),float(mp2_energies[i]),float(ccsd_energies[i]),float(fci_energies[i])))
  f.close()
