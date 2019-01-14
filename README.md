# Dirac interfaced with Openfermion

# Requirements

- Dirac 
- Openfermion
- all requirements for the aforementioned programs.

# Installation

INSTALL DIRAC: http://diracprogram.org/doku.php

INSTALL OPENFERMION : https://github.com/quantumlib/OpenFermion

INSTALL Openfermion-Dirac interface:

Clone the complete repository:
```
$ git clone https://github.com/bsenjean/Openfermion-Dirac.git
```

Build the code:
```
$ cd /path/to/Openfermion-Dirac/
$ pip install -e .
$ cd utils/
$ gfortran dirac_openfermion_mointegral_export.F90 -o dirac_openfermion_mointegral_export.x
```

In /path/to/Openfermion-Dirac/openfermion_dirac/_run_dirac.py change the following :
- Set your own path to pam (pam is the run_script of the Dirac program), which is called in the subprocess,
or set directly the /path/to/dirac/build/pam to your bash_profile.
- Same for dirac_openfermion_mointegral_export.x which is in /path/to/Openfermion-Dirac/utils/, either set
your own path to it in the subprocess in _run_dirac.py, or set the path to your bash_profile.
(Note that dirac_openfermion_mointegral_export.x will be in the release of Dirac2019).

# Use

Different examples are furnished in the examples/ repository in python, as well as a tutorial in tutorial/. If one wants to play more with the tutorial, use jupyter notebook to do so.
