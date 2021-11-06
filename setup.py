from distutils.core import setup
from pkg_resources import parse_requirements
import os

install_reqs = parse_requirements(open('requirements.txt', 'r'))
reqs = [str(ir) for ir in install_reqs]

setup(name='NCBItools', version='0.1', py_modules=['NCBItools'], install_requires = reqs)
