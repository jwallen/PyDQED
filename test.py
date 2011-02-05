#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This file contains unit tests for PyDAS.
"""

import unittest

from pydqed import DQED
import math
import numpy

################################################################################

class Optimization1(DQED):
    """
    A simple optimization of the function f(x) = (x - 100)^4 with no 
    constraints.
    """
    
    def evaluate(self, x):
        Neq = self.Neq; Nvars = self.Nvars; Ncons = self.Ncons
        f = numpy.zeros((Neq), numpy.float64)
        J = numpy.zeros((Neq, Nvars), numpy.float64)
        fcons = numpy.zeros((Ncons), numpy.float64)
        Jcons = numpy.zeros((Ncons, Nvars), numpy.float64)
        
        f[0] = (x[0] - 100.)**4
        J[0,0] = 4 * (x[0] - 100.)**3
        
        return f, J, fcons, Jcons

class Optimization2(DQED):
    """
    An optimization of the parameters (a, b, c, d) of the equation
    
        f(t) = a*exp(b*t) + c*exp(d*t)

    and linear constraint
        
        0.05 <= b - d

    """
    
    def __init__(self, tdata, fdata):
        self.tdata = tdata
        self.fdata = fdata
    
    def evaluate(self, x):
        Neq = self.Neq; Nvars = self.Nvars; Ncons = self.Ncons
        f = numpy.zeros((Neq), numpy.float64)
        J = numpy.zeros((Neq, Nvars), numpy.float64)
        fcons = numpy.zeros((Ncons), numpy.float64)
        Jcons = numpy.zeros((Ncons, Nvars), numpy.float64)
        
        for i in range(Neq):
            f[i] = x[0] * math.exp(x[1]*self.tdata[i]) + x[2] * math.exp(x[3]*self.tdata[i]) - self.fdata[i]
        
        for i in range(Neq):
            J[i,0] = math.exp(x[1] * self.tdata[i])
            J[i,1] = x[0] * self.tdata[i] * math.exp(x[1] * self.tdata[i])
            J[i,2] = math.exp(x[3] * self.tdata[i])
            J[i,3] = x[2] * self.tdata[i] * math.exp(x[3] * self.tdata[i])
        
        fcons[0] = x[1] - x[3]
        Jcons[0,0] = 0.0
        Jcons[0,1] = 1.0
        Jcons[0,2] = 0.0
        Jcons[0,3] = -1.0
        
        return f, J, fcons, Jcons

################################################################################

class DQEDCheck(unittest.TestCase):
    """
    Contains unit tests of the DASSL wrapper.
    """
    
    def test1a(self):
        """
        Test the optimization of f(x) = (x - 100)^4 without bounds.
        """
        x0 = numpy.ones((1), numpy.float64)
        
        opt = Optimization1()
        opt.initialize(Nvars=1, Ncons=0, Neq=1, bounds=None, tolf=1e-16, told=1e-8, tolx=1e-8, maxIter=100)
        x, igo = opt.solve(x0)
        self.assertTrue(igo in [2,4,6,7], 'Unexpected return status %i from DQED' % igo)
        self.assertAlmostEqual(x[0] / 100.0, 1.0, 5) 
        
    def test1b(self):
        """
        Test the optimization of f(x) = (x - 100)^4 with an upper bound.
        """
        x0 = numpy.ones((1), numpy.float64)
        
        opt = Optimization1()
        opt.initialize(Nvars=1, Ncons=0, Neq=1, bounds=[(None,50)], tolf=1e-16, told=1e-8, tolx=1e-8, maxIter=100)
        x, igo = opt.solve(x0)
        self.assertTrue(igo in [2,4,6,7], 'Unexpected return status %i from DQED' % igo)
        self.assertAlmostEqual(x[0] / 50.0, 1.0, 5)
    
    def test1c(self):
        """
        Test the optimization of f(x) = (x - 100)^4 with a lower bound.
        """
        x0 = numpy.ones((1), numpy.float64)
        
        opt = Optimization1()
        opt.initialize(Nvars=1, Ncons=0, Neq=1, bounds=[(-50,None)], tolf=1e-16, told=1e-8, tolx=1e-8, maxIter=100)
        x, igo = opt.solve(x0)
        self.assertTrue(igo in [2,4,6,7], 'Unexpected return status %i from DQED' % igo)
        self.assertAlmostEqual(x[0] / 100.0, 1.0, 5)
    
    def test2(self):
        """
        An optimization of the parameters (a, b, c, d) of the equation
        
            f(t) = a*exp(b*t) + c*exp(d*t)

        given several pairs of values (t, f(t)) and subject to the bounds
        
            0 <= a
        -25.0 <= b <= 0
            0 <= c
        -25.0 <= d <= 0
        
        and linear constraint
            
            0.05 <= b - d

        """
        tdata = numpy.array([0.05, 0.1, 0.4, 0.5, 1.0], numpy.float64)
        fdata = numpy.array([2.206, 1.994, 1.350, 1.216, 0.7358], numpy.float64)
        bounds = [
            (0.0,None),
            (-25.0,0.0),
            (0.0,None),
            (-25.0,0.0),
            (0.05,None),
        ]
        x0 = numpy.zeros(4, numpy.float64)
        
        opt = Optimization2(tdata, fdata)
        opt.initialize(Nvars=4, Ncons=1, Neq=5, bounds=bounds, tolf=1e-5, told=1e-5, tolx=1e-5, maxIter=100)
        x, igo = opt.solve(x0)
        self.assertTrue(igo in [2,4,6,7], 'Unexpected return status %i from DQED' % igo)
        self.assertAlmostEqual(x[0] / 1.999475, 1.0, 4)
        self.assertAlmostEqual(x[1] / -0.999801, 1.0, 4)
        self.assertAlmostEqual(x[2] / 0.500057, 1.0, 4)
        self.assertAlmostEqual(x[3] / -9.953988, 1.0, 4)
        
################################################################################

if __name__ == '__main__':
    unittest.main( testRunner = unittest.TextTestRunner(verbosity=2) )
