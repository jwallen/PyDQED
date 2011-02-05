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
    cdef numpy.ndarray rwork, iwork, ropt, iopt, ind, bl, bu, x

    def __init__(self):
        self.Nvars = 0
        self.Ncons = 0
        self.Neq = 0
    
    cpdef initialize(self, int Neq, int Nvars, int Ncons):
        """
        Initialize the DQED solver.
        """
        
        self.Nvars = Nvars
        self.Ncons = Ncons
        self.Neq = Neq

        # Determine required sizes of work arrays
        Nt = 5                         # if not using option 15
        Na = Ncons + 2 * Nvars + Nt    # if not using option 14
        liwork = 3 * Ncons + 9 * Nvars + 4 * Nt + Na + 10 + 1
        lrwork = Na * (Na + 4 + 2) + Nvars * (Nt + 33) + (Neq + max(Neq,Nvars) + 14) * Nt + 9 * Ncons + 26 

        # Allocate work arrays
        self.rwork = numpy.zeros(lrwork, numpy.float64)
        self.iwork = numpy.zeros(liwork, numpy.int32)
        
        # Tell how much storage we gave the solver
        self.iwork[0] = lrwork
        self.iwork[1] = liwork

        # Allocate options arrays
        self.ropt = numpy.zeros(1, numpy.float64)
        self.iopt = numpy.zeros(1, numpy.int32)
        
        self.iopt[0] = 99       # No further options are changed

        # Set up bounds arrays
        self.ind = 4 * numpy.ones((Nvars), numpy.int32)
        self.bl = numpy.zeros((Nvars), numpy.float64)
        self.bu = numpy.zeros((Nvars), numpy.float64)    
    
    def solve(self, numpy.ndarray[numpy.float64_t,ndim=1] x0):
        """
        Using the initial guess `x0`, return the least-squares solution to the
        set of nonlinear algebraic equations defined by the :meth:`evaluate` 
        method of the derived class.        
        """
        
        # Make sure the length of the initial guess matches the expected number of variables
        if len(x0) != self.Nvars:
            raise DQEDError('Expected %i values of x0, got %i.' % (self.Nvars, len(x0)))
        self.x = x0
        
        # Set the global DQED object to this object (so we can get back to
        # this object's residual and jacobian methods
        global dqedObject
        dqedObject = self
        
        cdef int lrw = self.rwork.shape[0]
        cdef int liw = self.iwork.shape[0]
        cdef int ldfjac = self.Neq + self.Ncons
        cdef int igo = 0
        
        cdef numpy.ndarray[numpy.float64_t,ndim=1] rnorm, x
        cdef numpy.ndarray[numpy.float64_t,ndim=2] fjac
        fjac = numpy.zeros((self.Neq+self.Ncons, self.Nvars), numpy.float64)
        rnorm = numpy.zeros((self.Neq), numpy.float64)
        
        cdef void* dqedev = <void*> evaluate
        
        # Call DQED
        dqed_(
            dqedev,
            &(self.Neq),
            &(self.Nvars),
            &(self.Ncons),
            <int*> self.ind.data,
            <numpy.float64_t*> self.bl.data,
            <numpy.float64_t*> self.bu.data,
            <numpy.float64_t*> self.x.data,
            <numpy.float64_t*> fjac.data,
            &(ldfjac),
            <numpy.float64_t*> rnorm.data,
            &(igo),
            <int*> self.iopt.data,
            <numpy.float64_t*> self.ropt.data,
            <int*> self.iwork.data,
            <numpy.float64_t*> self.rwork.data,
        )
        
        # Unset the global DQED object
        dqedObject = None
        
        # Return the result to the user
        x = self.x
        return x, igo
        
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
    cdef numpy.ndarray[numpy.float64_t,ndim=1] res
    cdef numpy.ndarray[numpy.float64_t,ndim=2] J, Jcons
    cdef int i, j, Neq, Nvars, Ncons, M
    res, J, Jcons = dqedObject.evaluate(dqedObject.x)
    Neq = J.shape[0]
    Nvars = J.shape[1]
    Ncons = Jcons.shape[0]
    M = Neq + Ncons
    for i in range(Neq):
        fj[Nvars*M+Ncons+i] = res[i]
        for j in range(Nvars):
            fj[j*M+Ncons+i] = J[i,j]
    for i in range(Ncons):
        for j in range(Nvars):
            fj[j*M+i] = Jcons[i,j]
