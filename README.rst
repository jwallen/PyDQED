******************************************************************************
PyDQED - A Python wrapper for the DQED constrained nonlinear optimization code
******************************************************************************

Introduction
============

PyDQED provides a means for Python code to utilize DQED, a bounded, constrained 
nonlinear optimization code. DQED is 
`publicly <http://people.sc.fsu.edu/~jburkardt/f_src/dqed/dqed.html>`_
`available <http://www.netlib.org/opt/dqed.f>`_, and is distributed with PyDQED.
PyDQED provides a Python extension type for DQED, which in turn provides a
Pythonic means of setting the solver options, providing residual and jacobian
functions, and running the solver.

License
=======

Copyright (c) 2011 by Joshua W. Allen (jwallen@mit.edu).

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the 'Software'),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

Dependencies
============

PyDQED has been tested on Python versions 2.5 and 2.6. It may or may not work
for other Python versions.

There are several Python dependencies that you must install before installing 
PyDQED:

* `Python <http://www.python.org/>`_ (versions 2.5.x and 2.6.x are known to work)

* `NumPy <http://numpy.scipy.org/>`_ (version 1.3.0 or later is recommended)

* `Cython <http://www.cython.org/>`_ (version 0.12.1 or later is recommended)

In addition, you will also need a Fortran compiler and a C compiler that
produce object files that can interoperate. The ``gfortran`` and ``gcc`` 
compiles from the GNU Compiler Collection are known to work. On Windows the
`MinGW <http://www.mingw.org/>`_ compiler collection provides these compilers.

The code for DQED has been provided with the PyDQED package. The license for 
DQED is different than that of the PyDQED wrapper code. **You are responsible 
for knowing and abiding by all licenses associated with DQED as well as with 
PyDQED.**

Installation
============

Windows
-------

The provided batch script will compile DQED and the PyDQED wrapper code. This 
script presumes that you have the 32-bit version of the MinGW C and Fortran 
compilers installed. Once you have run the batch script, you can install PyDQED
into your Python packages if you desire by running the following command from 
the base package directory:

    > python setup.py install

Linux
-----

A Makefile has been provided that can be used to compile DQED and the PyDQED 
wrapper code. To use, invoke the following command from the base package 
directory::

    $ make

This command will build PyDQED in-place, rather than installing it to your
Python package directory. If you wish to formally install PyDQED, run the
following command from the base package directory after the ``make`` command
(you may need root privileges for this)::

    # python setup.py install

Mac
---

Very much like the linux, but if you have installed gfortran via homebrew,
(`brew install gfortran`) then you may need to tell the makefile where to 
find the libraries::

    $ export LIBRARY_PATH=/usr/local/lib/gcc
    $ make
    $ make install

If you have a new Lion (OS X 10.7) installation without gcc-4.2, your numpy
distutils will probably try to link with clang, which ignores the LIBRARY_PATH.
If you link with gcc, it may get confused by the i386 and x86_64 bundle.
You can fix both of these problems by first typing the line::

    $ export LDSHARED='gcc -bundle -undefined dynamic_lookup -arch x86_64' 


