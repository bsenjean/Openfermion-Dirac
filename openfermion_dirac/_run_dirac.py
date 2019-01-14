# OpenFermion plugin to interface with Dirac.
# Copyright 2017 The OpenFermion Developers.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""Functions to prepare Dirac input and run calculations. This program is inspired from _run_psi4.py of the OpenFermion-Psi4 interface."""
from __future__ import absolute_import

import os
import re
import subprocess
import warnings

class ActiveOrbitalsError(Exception):
    pass
class SpecialBasisError(Exception):
    pass
class SpeedOfLightError(Exception):
    pass

def create_geometry_string(geometry):
    """This function converts MolecularData geometry to Dirac geometry.
       This function is taken from Openfermion-Psi4 interface.

    Args:
        geometry: A list of tuples giving the coordinates of each atom.
            example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))]. Distances in
            angstrom. Use atomic symbols to specify atoms.

    Returns:
        geo_string: A string giving the geometry for each atom on a line, e.g.:
            H 0. 0. 0.
            H 0. 0. 0.7414
    """
    geo_string = ''
    for item in geometry:
        atom = item[0]
        coordinates = item[1]
        line = '{} {} {} {}'.format(atom,
                                    coordinates[0],
                                    coordinates[1],
                                    coordinates[2])
        if len(geo_string) > 0:
            geo_string += '\n'
        geo_string += line
    return geo_string


def generate_dirac_input(molecule,
                        symmetry,
                        run_scf,
                        run_ccsd,
                        relativistic,
                        speed_of_light,
                        active):
    """This function creates and saves a Dirac input file.

    Args:
        molecule: An instance of the MolecularData class.
        symmetry: Boolean to specify the use of symmetry
        run_scf: Boolean to run SCF calculation.
        run_ccsd: Boolean to run CCSD calculation. (note that SCF and MP2
                  energies are done as well with CCSD)
        relativistic: Boolean to specify relativistic calculation or not
        speed_of_light: Real value for the speed of light (137 a.u.) in atomic unit,
                        to be changed if wanted in order to increase or decrease relativistic
                        effects.
        active: A list of 3 real numbers select active orbitals.
                first number : lowest energy
                second number : highest energy
                third number : minimum gap required between the neighbor energy
                               of the lowest and highest energy set previously.
                               If the gap is lower than this number, the lowest
                               (or highest) energy is shifted down (or up) until
                               the gap is larger to the third number.

    Returns:
        input_file: A string giving the name of the saved input file, and the xyz file.
    """
    # Create Dirac geometry string.
    geo_string = create_geometry_string(molecule.geometry)
    xyz_file = molecule.filename + '.xyz'
    with open(xyz_file, 'w') as f:
     f.write(str(molecule.n_atoms)+'\n')
     f.write(molecule.filename + ' # anything can be in this line\n')
     f.write(geo_string)
    f.close()

    # Check if the keywords are well defined
    if active is False:
       pass
    elif not isinstance(active, list):
       raise ActiveOrbitalsError('keyword active should be a list')
    elif len(active) != 3:
       raise ActiveOrbitalsError('keyword active needs 3 numbers in a list')
    elif (active[0] >= active[1]):
       raise ActiveOrbitalsError('the second argument of active should be higher than the first argument')
 
    if molecule.basis == "special" and molecule.special_basis is None:
       raise SpecialBasisError('special_basis should be specified')
    elif molecule.basis == "special" and len(molecule.special_basis) != 2:
       raise SpecialBasisError('special_basis should be a list of two strings corresponding to the default and special basis, respectively')
    elif (molecule.basis != "special") and (molecule.special_basis is not None):
       raise SpecialBasisError('special_basis was specified without setting basis = "special"')

    if speed_of_light is not False and relativistic is False:
       raise SpeedOfLightError('A given speed of light has been specified without setting relativistic to True')

    # Write input file and return handle.
    input_file = molecule.filename + '.inp'
    with open(input_file, 'w') as f:
      f.write("**DIRAC\n")
      f.write(".4INDEX\n")
      f.write(".WAVE FUNCTION\n")
      f.write("**WAVE FUNCTION\n")
      f.write(".SCF\n")
      if run_ccsd:
       f.write(".RELCCSD\n")
      if relativistic is False:
       f.write("**HAMILTONIAN\n")
       f.write(".NONREL\n")
       f.write("**INTEGRALS\n")
       f.write(".NUCMOD\n")
       f.write(" 1\n")
      if speed_of_light is not False:
       f.write("**GENERAL\n")
       f.write(".CVALUE\n")
       f.write(" " + str(speed_of_light) + "\n")
      f.write("**MOLTRA\n")
      f.write(".MDCINT\n")
      if active is not False:
       f.write(".ACTIVE\n")
       f.write("energy " + str(active[0]) + " " + str(active[1]) + " " + str(active[2]) + "\n")
      f.write("**MOLECULE\n")
      f.write("*CHARGE\n")
      f.write(".CHARGE\n")
      f.write(" " + str(molecule.charge) + "\n")
      if symmetry is False:
       f.write("*SYMMETRY\n")
       f.write(".NOSYM\n")
      f.write("*BASIS\n")
      f.write(".DEFAULT\n")
      if molecule.basis == "special":
       f.write(molecule.special_basis[0] + "\n")
       f.write(".SPECIAL\n")
       f.write(molecule.special_basis[1] + "\n")
      else:
       f.write(molecule.basis + "\n")
      f.write("*END OF INPUT\n")
    f.close()

    return input_file, xyz_file

def rename(molecule):
    output_file_dirac = molecule.filename + "_" + molecule.name + '.out'
    output_file = molecule.filename + '.out'
    os.rename("FCIDUMP", "FCIDUMP_" + molecule.name)
    os.rename(output_file_dirac,output_file)

def clean_up(molecule, delete_input=True, delete_xyz=True, delete_output=False, delete_MRCONEE=True,
             delete_MDCINT=True, delete_FCIDUMP=False):
    os.remove("FCITABLE")
    input_file = molecule.filename + '.inp'
    xyz_file = molecule.filename + '.xyz'
    output_file_dirac = molecule.filename + "_" + molecule.name + '.out'
    output_file = molecule.filename + '.out'
    run_directory = os.getcwd()
    for local_file in os.listdir(run_directory):
        if local_file.endswith('.clean'):
            os.remove(run_directory + '/' + local_file)
    try:
        os.remove('timer.dat')
    except:
        pass
    if delete_input:
        os.remove(input_file)
    if delete_output:
        os.remove(output_file)
    if delete_xyz:
        os.remove(xyz_file)
    if delete_MRCONEE:
        os.remove("MRCONEE")
    if delete_MDCINT:
        os.remove("MDCINT")
    if delete_FCIDUMP:
        os.remove("FCIDUMP_" + molecule.name)


def run_dirac(molecule,
             symmetry=True,
             run_scf=True,
             run_ccsd=False,
             relativistic=False,
             speed_of_light=False,
             active=False,
             delete_input=False,
             delete_output=False,
             delete_xyz=False,
             delete_MRCONEE=False,
             delete_MDCINT=False,
             delete_FCIDUMP=False,
             save=False):
    """This function runs a Dirac calculation.

    Args:
        molecule: An instance of the MolecularData class.
        symmetry: Optional boolean to remove symmetry in the calculation
        run_scf: Optional boolean to run SCF calculation.
        run_ccsd: Optional boolean to run CCSD calculation.
        relativistic: Optional boolean to run relativistic calculation
        speed_of_light: Optional real to give another value to the speed of light
        active: Optional list of 3 real numbers select active orbitals.
                first number : lowest energy
                second number : highest energy
                third number : minimum gap required between the lowest (highest)
        delete_input: Optional boolean to delete Dirac input file.
        delete_output: Optional boolean to delete Dirac output file.
        delete_xyz: Optional boolean to delete Dirac xyz file.
        delete_MRCONEE: Optional boolean to delete Dirac MRCONEE file.
        delete_MDCINT: Optional boolean to delete Dirac MDCINT file.
        delete_FCIDUMP: Optional boolean to delete Dirac FCIDUMP file.

    Returns:
        molecule: The updated MolecularData object.
    """
    # Prepare input.
    input_file, xyz_file = generate_dirac_input(molecule,
                        symmetry,
                        run_scf,
                        run_ccsd,
                        relativistic,
                        speed_of_light,
                        active)

    # Run Dirac
    print('Starting Dirac calculation\n')
    subprocess.check_call('pam --mol=' + xyz_file + ' --inp=' + input_file + ' --get="MRCONEE MDCINT" --silent --noarch > output_script', shell=True)

    # run dirac_openfermion_mointegral_export.x
    print('Creation of the FCIDUMP file\n')
    subprocess.check_call("dirac_openfermion_mointegral_export.x >> output_script",shell=True)

    rename(molecule)

    if save:
     try:
        print("Saving the results")
        molecule.save()
     except:
        warnings.warn('Error in saving results.',
                      Warning)

    # Clean-up
    clean_up(molecule, delete_input, delete_xyz, delete_output, delete_MRCONEE, delete_MDCINT, delete_FCIDUMP)
    return molecule
