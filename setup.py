from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "cEM_zip",
    ext_modules = cythonize("Jazzlib/*.pyx"),
)

#python setup.py build_ext --inplace
