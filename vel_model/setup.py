from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

setup(
      name='vel_model_cython',
      ext_modules=cythonize("vel_model_cython.pyx",
      compiler_directives={'cpow' : True,
                            'language_level' : "3",  # added for python3
                           }
      ),
      include_dirs=[numpy.get_include()]
      )