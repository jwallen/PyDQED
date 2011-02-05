************
Installation
************

Prerequisites
=============

PyDQED is currently available for the `Python <http://www.python.org/>`_ 2.x 
series. In particular, Python 2.5 and later are known to work. In addition to
the packages in the Python standard library, PyDQED depends on the following 
packages:

* `NumPy <http://numpy.scipy.org/>`_ version 1.3.0 or later

* `Cython <http://www.cython.org/>`_ version 0.12.1 or later

In addition, you will also need a Fortran compiler and a C compiler that
produce object files that can interoperate [#f1]_. The ``gfortran`` and ``gcc`` 
compilers from the `GNU Compiler Collection <http://gcc.gnu.org/>`_ are known 
to work, and are probably your best bet no matter what operating system you 
are using. On Windows the `MinGW <http://www.mingw.org/>`_ compiler collection 
provides these compilers.

The DQED code is provided with the PyDQED package; you do not need to download 
it separately.

.. [#f1] The Fortran interfaces are exposed to Python via C, so the installer
    needs to be able to link object files from Fortran and C for this to work.

Installing PyDQED
=================

If you are running an operating system other than Windows, refer to the 
section directly below. Windows users get their own special installation
procedure, described subsequently.

Unix-like Systems
-----------------

A Makefile has been provided that can be used to compile DQED and the PyDQED 
wrapper code. To use, invoke the following command from the root package 
directory::

    $ make

This command will build PyDQED in-place, rather than installing it to your
Python package directory. At this point you can optionally run the unit test 
script ``test.py`` and/or any of the provided examples to confirm that PyDQED
was compiled successfully.

If you wish to formally install PyDQED, run the following command from the root 
package directory after the ``make`` command (you may need root privileges for 
this)::

    $ python setup.py install

You may wish to write a file `make.inc` that sets certain variables used by
the Makefiles (e.g. the Fortran compiler). An example of such a file, 
`make.inc.example`, has been provided.

Windows
-------

.. warning:: 

    A regression in Cython 0.14 has resulted in it being unusable for PyDQED
    in Windows. Cython 0.13 is known to work, and more recent released of
    Cython should also correct for this regression. See
    `this thread <http://www.mail-archive.com/cython-dev@codespeak.net/msg10367.html>`_
    from the ``cython-dev`` mailing list for more information.

A batch script ``make.bat`` has been provided in the root package directory.
Double-clicking this script will cause DQED to be compiled into a static 
library and the PyDQED wrapper code to be compiled. 

.. note:: 
    
    The batch script presumes that you have the 32-bit version of the MinGW
    C and Fortran compilers installed. You may need to add the location of
    the MinGW ``bin`` directory to the ``PATH`` environment variable if this
    has not already been done.

At this point you can optionally run the unit test script ``test.py`` and/or 
any of the provided examples to confirm that PyDQED was compiled successfully.

.. warning::

    When using MinGW with Cython on Windows, you may encounter the error
    "Unable to find vcvarsall.bat". A workaround for this issue is available
    from the Cython FAQ at
    `this page <http://wiki.cython.org/FAQ#HowdoIworkaroundthe.22unabletofindvcvarsall.bat.22errorwhenusingMinGWasthecompiler.28onWindows.29.3F>`_.
    In particular the ``pydistutils.cfg`` file approach should work.

If you wish to formally install PyDQED, run the following command from the root 
package directory after the batch script completes successfully (you may need
administrator privileges for this)::

    > python setup.py install

