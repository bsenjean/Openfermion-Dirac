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
                        run_fci,
                        openshell,
                        relativistic,
                        NONREL,
                        point_nucleus,
                        uncontract,
                        speed_of_light,
                        active,
                        properties,
                        operator,
                        propint,
                        fcidump,
                        manual_option):
    """This function creates and saves a Dirac input file.

    Args:
        molecule: An instance of the MolecularData class.
        symmetry: Boolean to specify the use of symmetry
        run_ccsd: Boolean to run CCSD calculation. (note that SCF and MP2
                  energies are done as well with CCSD)
        run_fci: Boolean to run FCI calculation.
        openshell: Perform a Average-of-configuration open-shell Hartree-Fock.
                   List of 3 values, 
                   - the first one for the number of closed-shell electrons
                   - the second one for the number of open-shell
                   - the third is for the fractional occupation (#active electrons/#active spinor)
        relativistic: Boolean to specify relativistic calculation or not
        NONREL: Use NONREL Hamiltonian instead of Levy-Leblond.
        point_nucleus: Boolean to specify the use of the nuclear model of point nucleus,
                        instead of Gaussian charge distribution (default).
        uncontract : Boolean to specify the use of uncontracted basis set
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
          or active can also be a string
          (the first for the starting active orbital and the second for the last active orbital).
          The program will check if the list has size two or three, depending on the input provided.
        properties: list of strings to ask for calculation of specific molecular properties
        operator: list of list of operators and their prefactor to be added in the Hamiltonian (like external electric field...)
        propint: Property integrals transformation 
        fcidump: Boolean if the FCIDUMP file has to be written or not. If not, the MRCONEE and MDCINT files are not needed and one can skip 2- and 4-index transformation.

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
    elif not (isinstance(active, list) and (len(active) == 3)) and not isinstance(active, str):
       raise ActiveOrbitalsError('keyword active should be a list of three real numbers or a string (see tutorial)')
 
    if molecule.basis == "special" and molecule.special_basis is None:
       raise SpecialBasisError('special_basis should be specified')
    elif molecule.basis == "special" and len(molecule.special_basis) != 2:
       raise SpecialBasisError('special_basis should be a list of two strings corresponding to the default and special basis, respectively')
    elif (molecule.basis != "special") and (molecule.special_basis is not None):
       raise SpecialBasisError('special_basis was specified without setting basis = "special"')

    if speed_of_light is not False and relativistic is False:
       raise SpeedOfLightError('A given speed of light has been specified without setting relativistic to True')

    if run_fci:
       warnings.warn('run_fci is simply implemented to get the FCI ground-state energy for now. Using it to get anything else is highly not recommended.\
 This option has not been tested meticulously and it might lead to errors depending on the system that you want to treat.',
                      Warning)

    nelec=int(molecule.n_electrons)
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
      if (openshell is not False):
         f.write("*SCF\n")
         f.write(".CLOSED SHELL\n"+openshell[0]+"\n")
         f.write(".OPEN SHELL\n"+openshell[1]+"\n"+openshell[2]+"\n")
      if run_fci:
       f.write(".DIRRCI\n")
      f.write("**HAMILTONIAN\n")
      if not relativistic:
       if NONREL:
          f.write(".NONREL\n")
       else:
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
      f.write("**INTEGRALS\n")
      if point_nucleus:
       f.write(".NUCMOD\n")
       f.write(" 1\n")
      if uncontract:
       f.write("*READIN\n")
       f.write(".UNCONTRACT\n")
      if speed_of_light is not False:
       f.write("**GENERAL\n")
       f.write(".CVALUE\n")
       f.write(" " + str(speed_of_light) + "\n")
      f.write("**MOLTRA\n")
      if run_ccsd or fcidump or run_fci:
        f.write(".MDCINT\n")
      else:
        # specifically ask not to skip the 4-index transformation of the effective Fock matrix. Only if fcidump is not True and CCSD and FCI are not performed.
        f.write(".NO4IND\n")
      if propint is not False:
       f.write(".PRPTRA\n")
      if run_fci:
       f.write(".SCHEME\n")
       f.write(" 4\n")
      f.write(".ACTIVE\n")
      if active is not False:
       if isinstance(active,str):
         f.write(active + "\n")
       elif isinstance(active,list):
         f.write("energy " + str(active[0]) + " " + str(active[1]) + " " + str(active[2]) + "\n")
      else:
       f.write(" all\n")
      if propint is not False:
       f.write("*PRPTRA\n")
       f.write(".OPERATOR\n")
       f.write(" " + propint + "\n")
      f.write("**MOLECULE\n")
      f.write("*CHARGE\n")
      f.write(".CHARGE\n")
      f.write(" " + str(molecule.charge) + "\n")
      if not symmetry or run_fci:
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
      if run_fci:
       if nelec%2==0:
        f.write("\n &RASORB  NELEC=" + str(nelec) + ", NRAS1=" + str(int(nelec/2)) + "," + str(int(nelec/2)) + ", MAXH1=" + str(nelec) + ", MAXE3=" + str(nelec) + "  &END\n")
       else:
        f.write("\n &RASORB  NELEC=" + str(nelec) + ", NRAS1=" + str(int(nelec//2+1)) + "," + str(int(nelec//2)) + ", MAXH1=" + str(nelec) + ", MAXE3=" + str(nelec) + "  &END\n")
       f.write(" &CIROOT  NROOTS=1 &END\n")
       f.write(" &DIRECT  MAXITER=20 &END\n")

    return input_file, xyz_file

def rename(molecule,fcidump,propint):
    output_file_dirac = molecule.filename + "_" + molecule.name + '.out'
    output_file = molecule.filename + '.out'
    if fcidump:
     os.rename(molecule.data_directory + "/" + "FCIDUMP", molecule.data_directory + "/" + "FCIDUMP_" + molecule.name)
    if propint is not False:
     os.rename(molecule.data_directory + "/" + "PROPINT", molecule.data_directory + "/" + "PROPINT_" + molecule.name)
    os.rename(output_file_dirac,output_file)

def clean_up(molecule, fcidump, delete_input, delete_xyz, delete_output, delete_MRCONEE,
             delete_MDCINT, delete_MDPROP):
    input_file = molecule.filename + '.inp'
    xyz_file = molecule.filename + '.xyz'
    output_file = molecule.filename + '.out'
    if delete_input:
      try:
        os.remove(input_file)
      except:
        pass
    if delete_output:
      try:
        os.remove(output_file)
      except:
        pass
    if delete_xyz:
      try:
        os.remove(xyz_file)
      except:
        pass
    if delete_MRCONEE:
      try:
        os.remove("MRCONEE")
      except:
        pass
    if delete_MDCINT:
      try:
        os.remove("MDCINT")
      except:
        pass
    if delete_MDPROP:
      try:
        os.remove("MDPROP")
      except:
        pass
    if fcidump:
      try:
        os.remove("FCITABLE")
      except:
        pass

def run_dirac(molecule,
             symmetry=True,
             run_dft=False,
             run_ccsd=False,
             run_fci=False,
             openshell=False,
             fcidump=False,
             relativistic=False,
             NONREL=False,
             point_nucleus=False,
             uncontract=False,
             speed_of_light=False,
             active=False,
             save=False,
             properties=False,
             operator=False,
             propint=False,
             manual_option=False,
             get=False,
             restart=False,
             delete_input=False,
             delete_output=False,
             delete_xyz=False,
             delete_MRCONEE=False,
             delete_MDCINT=False,
             delete_MDPROP=False):
    """This function runs a Dirac calculation.

    Args:
        molecule: An instance of the MolecularData class.
        symmetry: Optional boolean to remove symmetry in the calculation
        run_dft: Optional boolean to run DFT calculation.
        run_ccsd: Optional boolean to run CCSD calculation.
        run_fci: Optional boolean to run FCI calculation.
        fcidump: Optional boolean to create the FCIDUMP file containing
                 all the one- and two-electron integrals.
                 Necessary to construct the molecular Hamiltonian
                 in OpenFermion.
        relativistic: Optional boolean to run relativistic calculation
        NONREL: Use NONREL Hamiltonian instead of Levy-Leblond.
        point_nucleus : Boolean to specify the use of the nuclear model of point nucleus,
                        instead of Gaussian charge distribution (default).
        speed_of_light: Optional real to give another value to the speed of light
        active: A list of 3 real numbers select active orbitals.
                first number : lowest energy
                second number : highest energy
                third number : minimum gap required between the neighbor energy
                               of the lowest and highest energy set previously.
                               If the gap is lower than this number, the lowest
                               (or highest) energy is shifted down (or up) until
                               the gap is larger to the third number.
          or active can also be a string
          (the first for the starting active orbital and the second for the last active orbital).
          The program will check if the list has size two or three, depending on the input provided.
        save: Boolean to save calculation in a hdf5 file.
        properties: list of strings to ask for calculation of specific molecular properties
        operator: list of list of operators and their prefactor to be added in the Hamiltonian (like external electric field...)
        propint: Property integrals transformation 
        get: string to ask specific files from DIRAC, like DFCOEF.
        restart: restart a calculation with a DFCOEF file (binary file containing all MO coeffs and basis functions)
        delete_input: Optional boolean to delete Dirac input file.
        delete_output: Optional boolean to delete Dirac output file.
        delete_xyz: Optional boolean to delete Dirac xyz file.
        delete_MRCONEE: Optional boolean to delete Dirac MRCONEE file.
        delete_MDCINT: Optional boolean to delete Dirac MDCINT file.
        delete_MDPROP: Optional boolean to delete Dirac MDPROP file.

    Returns:
        molecule: The updated MolecularData object.
    """
    if run_fci and fcidump:
       raise InputError('Both run_fci and fcidump cannot be True. For now, run_fci should just be used to generate the Dirac-output.')

    # Prepare input.
    input_file, xyz_file = generate_dirac_input(molecule,
                        symmetry,
                        run_dft,
                        run_ccsd,
                        run_fci,
                        openshell,
                        relativistic,
                        NONREL,
                        point_nucleus,
                        uncontract,
                        speed_of_light,
                        active,
                        properties,
                        operator,
                        propint,
                        fcidump,
                        manual_option)

    # Run Dirac
    get_all = ""
    if fcidump:
       get_all += "MRCONEE MDCINT"
    if propint is not False:
       get_all += " MDPROP"
    if get is not False:
       get_all += " "+get
    if not restart:
        subprocess.check_call("pam --mol=" + xyz_file + " --inp=" + input_file + " --get='" + get_all + "' --silent --noarch", shell=True, cwd=molecule.data_directory)
    else:
        subprocess.check_call("pam --mol=" + xyz_file + " --inp=" + input_file + " --get='" + get_all + "' --put='DFCOEF' --silent --noarch", shell=True, cwd=molecule.data_directory)

    # run dirac_mointegral_export.x
    if fcidump:
      subprocess.check_call("dirac_mointegral_export.x fcidump",shell=True,cwd=molecule.data_directory)
    
    if propint is not False:
      subprocess.check_call("dirac_mointegral_export.x propint " + str(propint),shell=True,cwd=molecule.data_directory)

    rename(molecule,fcidump,propint)

    if save:
     try:
        molecule.save()
     except:
        warnings.warn('Error in saving results. Check that fcidump is set to True.',
                       Warning)

    # Clean-up
    clean_up(molecule, fcidump, delete_input, delete_xyz, delete_output, delete_MRCONEE, delete_MDCINT, delete_MDPROP)
    return molecule
