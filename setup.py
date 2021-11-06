from distutils.core import setup
from pip.req import parse_requirements
import os

install_reqs = parse_requirements(os.path.dirname(__file__))
reqs = [str(ir.req) for ir in install_reqs]

setup(name='NCBItools', version='0.1', py_modules=['NCBItools'], install_requires = reqs)
