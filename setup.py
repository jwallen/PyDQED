#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   PyDQED - A Python wrapper for the DQED constrained nonlinear optimization code
#
#   Copyright (c) 2011 by Joshua W. Allen (joshua.w.allen@gmail.com)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import numpy

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# Turn on HTML annotation file generation
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

# The Cython extension modules to compile
ext_modules = [
    Extension(
        'pydqed', 
        ['pydqed.pyx'], 
        include_dirs=['.', numpy.get_include()], 
        libraries=['gfortran'], 
        extra_objects=['build/dqed.o'],
    ),
]

# Run the setup command
setup(name='PyDQED',
    version='0.1.0',
    description='A Python wrapper to the DQED bounded constrained nonlinear optimization code',
    author='Joshua W. Allen',
    author_email='joshua.w.allen@gmail.com',
    url='http://github.com/jwallen/PyDQED',
    py_modules=['pydqed'],
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
