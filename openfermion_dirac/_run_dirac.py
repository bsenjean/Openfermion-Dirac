#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

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
class PropertyError(Exception):
    pass
class InputError(Exception):
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
                        run_dft,
                        run_ccsd,
                        relativistic,
                        point_nucleus,
                        speed_of_light,
                        active,
                        properties,
                        operator,
                        propint,
                        manual_option):
    """This function creates and saves a Dirac input file.

    Args:
        molecule: An instance of the MolecularData class.
        symmetry: Boolean to specify the use of symmetry
        run_ccsd: Boolean to run CCSD calculation. (note that SCF and MP2
                  energies are done as well with CCSD)
        relativistic: Boolean to specify relativistic calculation or not
        point_nucleus : Boolean to specify the use of the nuclear model of point nucleus,
                        instead of Gaussian charge distribution (default).
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
      if properties is not False:
       f.write(".PROPERTIES\n")
      f.write("**WAVE FUNCTION\n")
      f.write(".SCF\n")
      if run_ccsd:
       f.write(".RELCCSD\n")
      f.write("**HAMILTONIAN\n")
      if not relativistic:
       f.write(".LEVY-LEBLOND\n")
      if run_dft is not False:
       f.write(".DFT\n")
       f.write(run_dft + "\n")
      if operator is not False:
       for i in range(len(operator)):
        f.write(".OPERATOR\n")
        for j in range(len(operator[i])):
         f.write(" " + str(operator[i][j]) + "\n")
      if properties is not False:
       f.write("**PROPERTY\n")
       for prop in properties:
        f.write("." + prop + "\n")
      if point_nucleus:
       f.write("**INTEGRALS\n")
       f.write(".NUCMOD\n")
       f.write(" 1\n")
      if speed_of_light is not False:
       f.write("**GENERAL\n")
       f.write(".CVALUE\n")
       f.write(" " + str(speed_of_light) + "\n")
      f.write("**MOLTRA\n")
      f.write(".MDCINT\n")
      if propint is not False:
       f.write(".PRPTRA\n")
      f.write(".ACTIVE\n")
      if active is not False:
       f.write("energy " + str(active[0]) + " " + str(active[1]) + " " + str(active[2]) + "\n")
      else:
       f.write(" all\n")
      if propint is not False:
#       f.write("# need to add a comment here for H2, don't know why\n")
#       f.write(".PRPTRA\n")
       f.write("*PRPTRA\n")
       f.write(".OPERATOR\n")
       f.write(" " + propint + "\n")
      f.write("**MOLECULE\n")
      f.write("*CHARGE\n")
      f.write(".CHARGE\n")
      f.write(" " + str(molecule.charge) + "\n")
      if not symmetry:
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
      if manual_option is not False:
       f.write(manual_option + "\n")
      f.write("*END OF INPUT\n")

    return input_file, xyz_file

def rename(molecule,propint):
    output_file_dirac = molecule.name + "_" + molecule.name + '.out'
    output_file = molecule.filename + '.out'
    os.rename("FCIDUMP", molecule.data_directory + "/" + "FCIDUMP_" + molecule.name)
    if propint is not False:
     os.rename("PROPINT", molecule.data_directory + "/" + "PROPINT_" + molecule.name)
    os.rename(output_file_dirac,output_file)

def clean_up(molecule, delete_input=True, delete_xyz=True, delete_output=False, delete_MRCONEE=True,
             delete_MDCINT=True, delete_MDPROP=False, delete_FCIDUMP=False):
    os.remove("FCITABLE")
    input_file = molecule.filename + '.inp'
    xyz_file = molecule.filename + '.xyz'
    output_file_dirac = molecule.name + "_" + molecule.name + '.out'
    output_file = molecule.filename + '.out'
    run_directory = molecule.data_directory + "/"
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
    if delete_MDPROP:
        os.remove("MDPROP")
    if delete_FCIDUMP:
        os.remove(molecule.data_directory + "/" + "FCIDUMP_" + molecule.name)


def run_dirac(molecule,
             symmetry=True,
             run_dft=False,
             run_ccsd=False,
             relativistic=False,
             point_nucleus=False,
             speed_of_light=False,
             active=False,
             save=False,
             properties=False,
             operator=False,
             propint=False,
             manual_option=False,
             delete_input=False,
             delete_output=False,
             delete_xyz=False,
             delete_MRCONEE=False,
             delete_MDCINT=False,
             delete_MDPROP=False,
             delete_FCIDUMP=False):
    """This function runs a Dirac calculation.

    Args:
        molecule: An instance of the MolecularData class.
        symmetry: Optional boolean to remove symmetry in the calculation
        run_ccsd: Optional boolean to run DFT calculation.
        run_ccsd: Optional boolean to run CCSD calculation.
        relativistic: Optional boolean to run relativistic calculation
        point_nucleus : Boolean to specify the use of the nuclear model of point nucleus,
                        instead of Gaussian charge distribution (default).
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
        delete_MDPROP: Optional boolean to delete Dirac MDPROP file.
        delete_FCIDUMP: Optional boolean to delete Dirac FCIDUMP file.

    Returns:
        molecule: The updated MolecularData object.
    """
    # Prepare input.
    input_file, xyz_file = generate_dirac_input(molecule,
                        symmetry,
                        run_dft,
                        run_ccsd,
                        relativistic,
                        point_nucleus,
                        speed_of_light,
                        active,
                        properties,
                        operator,
                        propint,
                        manual_option)

    # Run Dirac
    print('Starting Dirac calculation\n')
    if propint is not False:
      subprocess.check_call("pam --mol=" + xyz_file + " --inp=" + input_file + " --get='MRCONEE MDCINT MDPROP' --silent --noarch", shell=True)
    else:
      subprocess.check_call("pam --mol=" + xyz_file + " --inp=" + input_file + " --get='MRCONEE MDCINT' --silent --noarch", shell=True)

    # run dirac_openfermion_mointegral_export.x
    print('\nCreation of the FCIDUMP file\n')
    subprocess.check_call("dirac_openfermion_mointegral_export.x fcidump",shell=True)
    
    if propint is not False:
      print('\nCreation of the PROPINT file\n')
      subprocess.check_call("dirac_openfermion_mointegral_export.x propint " + str(propint),shell=True)

    rename(molecule,propint)

    if save:
     try:
        print("\nSaving the results\n")
        molecule.save()
     except:
        warnings.warn('Error in saving results.',
                      Warning)

    # Clean-up
    clean_up(molecule, delete_input, delete_xyz, delete_output, delete_MRCONEE, delete_MDCINT, delete_MDPROP, delete_FCIDUMP)
    return molecule
