from setuptools import find_packages, setup

setup(name='openfermion_interface',
      version='0.1',
      description='OpenFermion/Dirac interface',
      url='https://github.com/bsenjean/Openfermion-Dirac',
      author='Bruno',
      author_email='senjean@lorentz.leidenuniv.nl',
      license='Apache 2',
      packages=['openfermion_dirac'],
      install_requires=["openfermion"]
      )
