################################################################################
#
#   PyDQED - A Python wrapper for the DQED constrained nonlinear optimization code
#
#   Copyright (c) 2011 by Joshua W. Allen (jwallen@mit.edu)
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

"""
This `Cython <http://www.cython.org/>`_ module exposes the DQED bounded 
constrained nonlinear optimization code to Python and provides a Python 
extension type, the :class:`DQED` base class, for working with DQED.

To use DQED, write a Python class or Cython extension type that derives from
the :class:`DQED` class and implement the :meth:`evaluate` method.
Run by calling the :meth:`initialize` method to set the solver parameters, then 
by using the :meth:`solve` method to execute the optimization or solve.

You can implement your derived class in pure Python, but for a significant
speed boost consider using Cython to compile the module. You can see the
proper Cython syntax for the residual and jacobian methods by looking at the
corresponding methods in the :class:`DASSL` base class.
"""

import math
import numpy
cimport numpy
cimport cython

################################################################################

# Expose the (double-precision) DQED function prototype to this module
cdef extern from "dqed.h":
    void dqed_(
        void* dqedev,           # The function that evaluates the nonlinear equations and constraints
        int* mequa,             # The number of equations to be solved
        int* nvars,             # The number of unknown variables
        int* mcon,              # The number of general constraints (not including simple bounds)
        int* ind,               # The type of simple bounds to use for each variable
        double* bl,             # The lower bound to use for each variable (if applicable)
        double* bu,             # The upper bound to use for each variable (if applicable)
        double* x,              # The initial guess (on input) or final solution (on output)
        double* fjac,           # The final values of the Jacobian matrix for the constrants and equations
        int* ldfjac,            # The leading dimension of the Jacobian matrix
        double* fnorm,          # The value of the Euclidean norm at the solution
        int* igo,               # A flag containing the current or final status of the optimization/solve
        int* iopt,              # Integer parameters modifying how DQED executes
        double* ropt,           # Double-precision parameters modifying how DQED executes
        int* iwork,             # Work space for integer values
        double* rwork           # Work space for double-precision values
    )

################################################################################

class DQEDError(Exception):
    """
    An exception class for exceptions relating to use of DQED. Pass a string
    message describing the circumstances of the exceptional behavior.
    """
    pass

################################################################################

cdef class DQED:
    """
    A base class for using the DQED bounded constrained nonlinear optimization 
    code by R. Hanson and F. Krogh.
    """
    
    cdef public int Nvars, Ncons, Neq
    cdef numpy.ndarray rwork, iwork, ropt, iopt, ind, bl, bu, x, fnorm, fjac

    def __init__(self):
        self.Nvars = 0
        self.Ncons = 0
        self.Neq = 0
    
    cpdef initialize(self, int Neq, int Nvars, int Ncons, list bounds=None, 
        double tolf=1e-5, double told=1e-5, double tolx=1e-5, int maxIter=100, 
        bint verbose=False):
        """
        Initialize the DQED solver. The required parameters are:

        * `Neq` - The number of algebraic equations.
        
        * `Nvars` - The number of unknown variables.
        
        * `Ncons` - The number of constraint equations.
        
        The optional parameters are:
        
        * `bounds` - A list of 2-tuples giving the lower and upper bound for
          each unknown variable. Use ``None`` if there is no bound in one or
          either direction. If provided, you must give bounds for every unknown
          variable.
        
        * `tolf` - The tolerance used for stopping when the norm of the 
          residual has absolute length less than `tolf`, i.e.
          :math:`\| \vector{f} \| \le \epsilon_f.
        
        * `told` - The tolerance used for stopping when changes to the unknown
          variables has absolute length less than `told`, i.e.
          :math:`\| \Delta \vector{x} \| \le \epsilon_d.
        
        * `tolx` - The tolerance used for stopping when changes to the unknown
          variables has relative length less than `tolx`, i.e.
          :math:`\| \Delta \vector{x} \| \le \epsilon_x \cdot \| \vector{x} \|`.
        
        * `maxIter` - The maximum number of iterations to use
        
        * `verbose` - ``True`` to have DQED print extra information about the
          solve, ``False`` to only see printed output when the solver has an
          error.
        
        """
        
        self.Nvars = Nvars
        self.Ncons = Ncons
        self.Neq = Neq

        # Determine required sizes of work arrays
        Nt = 5                             # if not using option 15
        Na = Ncons + 2 * Nvars + Nt + 1    # if not using option 14
        liwork = 3 * Ncons + 9 * Nvars + 4 * Nt + Na + 10
        lrwork = Na * (Na + 6) + Nvars * (Nt + 33) + (Neq + max(Neq,Nvars) + 14) * Nt + 9 * Ncons + 26 + 12

        # Allocate work arrays
        self.rwork = numpy.zeros(lrwork, numpy.float64)
        self.iwork = numpy.zeros(liwork, numpy.int32)
        
        # Tell how much storage we gave the solver
        self.iwork[0] = lrwork
        self.iwork[1] = liwork

        # Allocate options arrays
        self.ropt = numpy.zeros(3, numpy.float64)
        self.iopt = numpy.zeros(13, numpy.int32)
        
        self.iopt[0] = 1         # Change the verbosity of the printed output
        self.iopt[1] = 1 if verbose else 0   # The desired output verbosity
        
        self.iopt[2] = 2         # Change the number of iterations
        self.iopt[3] = maxIter   # Maximum number of iterations
        
        self.iopt[4] = 4         # Set the option to change the value of TOLF
        self.iopt[5] = 1         # Where in ROPT to look for the new value
        self.ropt[0] = tolf      # New value for TOLF
        
        self.iopt[6] = 5         # Set the option to change the value of TOLX
        self.iopt[7] = 2         # Where in ROPT to look for the new value
        self.ropt[1] = tolx      # New value for TOLX
        
        self.iopt[8] = 6         # Set the option to change the value of TOLD
        self.iopt[9] = 3         # Where in ROPT to look for the new value
        self.ropt[2] = told      # New value for TOLD
        
        self.iopt[10] = 17       # Do not allow the flag IGO to return the value IGO=3
        self.iopt[11] = 1        # Forces a full model step
        self.iopt[12] = 99       # No further options are changed

        # Set up bounds arrays
        bounds = bounds or [(None,None) for i in range(Nvars)]
        if len(bounds) != Nvars + Ncons:
            raise DQEDError('If bounds are specified, they must be specified for every unknown variable.')
        self.ind = numpy.zeros((Nvars+Ncons), numpy.int32)
        self.bl = numpy.zeros((Nvars+Ncons), numpy.float64)
        self.bu = numpy.zeros((Nvars+Ncons), numpy.float64)
        for i in range(len(bounds)):
            bl, bu = bounds[i]
            if bl is not None and bu is None:
                self.ind[i] = 1
                self.bl[i] = bl
            elif bl is None and bu is not None:
                self.ind[i] = 2
                self.bu[i] = bu
            elif bl is not None and bu is not None:
                self.ind[i] = 3
                self.bl[i] = bl
                self.bu[i] = bu
            else:
                self.ind[i] = 4
    
    
    def solve(self, numpy.ndarray[numpy.float64_t,ndim=1] x0):
        """
        Using the initial guess `x0`, return the least-squares solution to the
        set of nonlinear algebraic equations defined by the :meth:`evaluate` 
        method of the derived class. This is the method that actually conducts
        the call to DQED. Returns the solution vector and a flag indicating
        the status of the solve. The possible output values of the flag are:
        
        ======== ===============================================================
        Value    Meaning
        ======== ===============================================================
        2        The norm of the residual is zero; the solution vector is a root of the system
        3        The bounds on the trust region are being encountered on each step; the solution vector may or may not be a local minimum
        4        The solution vector is a local minimum
        5        A significant amount of noise or uncertainty has been observed in the residual; the solution may or may not be a local minimum
        6        The solution vector is only changing by small absolute amounts; the solution may or may not be a local minimum
        7        The solution vector is only changing by small relative amounts; the solution may or may not be a local minimum
        8        The maximum number of iterations has been reached; the solution is the best found, but may or may not be a local minimum
        9-18     An error occurred during the solve operation; the solution is not a local minimum 
        ======== ===============================================================        
        
        """
        
        cdef int igo
        cdef numpy.ndarray[numpy.float64_t,ndim=1] x
        
        self.fjac = numpy.zeros((self.Neq+self.Ncons, self.Nvars+1), numpy.float64, order='F')
        self.fnorm = numpy.zeros((self.Neq), numpy.float64)
        
        # Make sure the length of the initial guess matches the expected number of variables
        if len(x0) != self.Nvars:
            raise DQEDError('Expected %i values of x0, got %i.' % (self.Nvars, len(x0)))
        self.x = x0
        
        # Call DQED
        igo = self.dqed(self.Neq, self.Nvars, self.Ncons, 
            self.ind, self.bl, self.bu, self.x, self.fnorm, self.fjac,
            self.iopt, self.ropt, self.iwork, self.rwork)
        
        # Return the result to the user
        x = self.x
        return x, igo
    
    cdef int dqed(self, int Neq, int Nvars, int Ncons, 
        numpy.ndarray[numpy.int32_t,ndim=1] ind,
        numpy.ndarray[numpy.float64_t,ndim=1] bl,
        numpy.ndarray[numpy.float64_t,ndim=1] bu,
        numpy.ndarray[numpy.float64_t,ndim=1] x,
        numpy.ndarray[numpy.float64_t,ndim=1] fnorm,
        numpy.ndarray[numpy.float64_t,ndim=2] fjac,
        numpy.ndarray[numpy.int32_t,ndim=1] iopt,
        numpy.ndarray[numpy.float64_t,ndim=1] ropt,
        numpy.ndarray[numpy.int32_t,ndim=1] iwork,
        numpy.ndarray[numpy.float64_t,ndim=1] rwork):
        
        cdef int ldfjac = self.Neq + self.Ncons
        cdef int igo = 0
        cdef void* dqedev = <void*> evaluate
        
        # Set the global DQED object to this object (so we can get back to
        # this object's residual and jacobian methods
        global dqedObject
        dqedObject = self
        
        dqed_(
            dqedev,
            &(Neq),
            &(Nvars),
            &(Ncons),
            <int*> ind.data,
            <numpy.float64_t*> bl.data,
            <numpy.float64_t*> bu.data,
            <numpy.float64_t*> x.data,
            <numpy.float64_t*> fjac.data,
            &(ldfjac),
            <numpy.float64_t*> fnorm.data,
            &(igo),
            <int*> iopt.data,
            <numpy.float64_t*> ropt.data,
            <int*> iwork.data,
            <numpy.float64_t*> rwork.data,
        )
        
        # Unset the global DQED object
        dqedObject = None
        
        return igo
    
    def evaluate(self, numpy.ndarray[numpy.float64_t,ndim=1] x):
        """
        Evaluate the nonlinear equations and constraints for this system, and
        the corresponding Jacobian matrices, at the given value of the solution
        vector `x`. Return a tuple containing three items:
        
        * A vector of the current values of the system of equations 
          :math:`\vector{f}(\vector{x})`.
          
        * A matrix of the current values of the Jacobian of the system of
          equations: :math:`J_{ij} = \frac{\partial f_i}{\partial x_j}`.
          
        * A matrix of the current values of the Jacobian of the (linear)
          constrains: :math:`J_{ij}^\prime = \frac{\partial g_i}{\partial x_j}`.
        
        """
        print('DQEDError: You must implement the evaluate() method in your derived class.')
        return None
        
################################################################################

# A module-level variable that contains the currently active DASSL object
# The residual and jacobian functions use this to call the object's residual
# and jacobian functions
cdef DQED dqedObject

cdef void evaluate(double* x, double* fj, int* ldfj, int* igo, int* iopt, double* ropt):
    
    cdef numpy.ndarray[numpy.float64_t,ndim=1] f, fcons
    cdef numpy.ndarray[numpy.float64_t,ndim=2] J, Jcons
    cdef int i, j, Neq, Nvars, Ncons
    
    f, J, fcons, Jcons = dqedObject.evaluate(dqedObject.x)
    
    Neq = J.shape[0]
    Nvars = J.shape[1]
    Ncons = Jcons.shape[0]
    for i in range(Neq):
        dqedObject.fjac[Ncons+i,Nvars] = f[i]
        for j in range(Nvars):
            dqedObject.fjac[Ncons+i,j] = J[i,j]
    for i in range(Ncons):
        dqedObject.fjac[i,Nvars] = fcons[i]
        for j in range(Nvars):
            dqedObject.fjac[i,j] = Jcons[i,j]
