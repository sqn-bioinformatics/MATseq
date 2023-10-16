from setuptools import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize("removes_PCR_duplicates.pyx"))
