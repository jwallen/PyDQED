************
Introduction
************

In the course of scientific computing, it is frequently desired to find the
solution of a set of nonlinear algebraic equations -- or at least the optimum
value -- subject to a number of constraints. Such a problem is often posed in
the following general way: minimize the sum of squares of the nonlinear
algebraic equations

.. math:: f_i(\vector{x}) = 0 \hspace{30pt} i = 1, 2, \ldots, N_\mathrm{eq}

for a vector :math:`\vector{x}` of :math:`N_\mathrm{vars}` variables, subject 
to the constraints

.. math:: g_j(\vector{x}) = 0 \hspace{30pt} j = 1, 2, \ldots, N_\mathrm{cons}

and optional lower and upper bounds

.. math:: l_k \le x_k \le u_k \hspace{30pt} k = 1, 2, \ldots, N_\mathrm{vars}

A variety of interesting, practical problems can be expressed in the above 
form.

The DQED solver provides an algorithm for solving the above problem for the case
where the constraints :math:`g_j(\vector{x})` are linear (or very nearly so).
For more on DQED, please see [Dongarra1979]_ [Schnabel1984]_ [Hanson1986]_ or
`this page <http://people.sc.fsu.edu/~jburkardt/f_src/dqed/dqed.html>`_.

.. [Dongarra1979] J. Dongarra, J. Bunch, C. Moler, and P. Stewart.
    "LINPACK User's Guide." *SIAM (Society for Industrial and Applied Mathematics)* Philadelphia, 1979.

.. [Schnabel1984] R. Schnabel and P. Frank. "Tensor Methods for Nonlinear Equations."
    *SIAM J. Numer. Analysis* **21**, p. 815-843 (1984).

.. [Hanson1986] R. Hanson. "Least Squares with Bounds and Linear Constraints."
    *SIAM J. Sci. Stat. Comp.* **7**, p. 826-834 (1986).

Why PyDQED?
===========

DQED is written in Fortran 90, which is not much fun to code in, especially for
novice programmers. (In particular, getting data into and out of a Fortran 
program via file input and output is quite difficult and awkward.) However, the 
strength of Fortran is that it produces code that is very efficient to execute, 
which is often important when solving large systems of nonlinear equations.

Meanwhile, the Python programming language is much easier to program in. In
particular, Python comes with a large library of free, open source packages 
that provide a wide range of functionality, limiting the amount of work the 
programmer needs to do from scratch. A number of packages, including 
`NumPy <http://numpy.scipy.org/>`_, `SciPy <http://www.scipy.org/>`_, and
`matplotlib <http://matplotlib.sourceforge.net/>`_, replace much of the
functionality of numerical computing environments such as MATLAB. However,
the optimization code within SciPy is still in need of some refinement.

PyDQED provides Python programmers with access to DQED via a Python wrapper.

Should I Use PyDQED?
====================

Depending on exactly the type of problem you are solving, PyDQED may be a useful
tool. By no means is PyDQED the best tool for every problem, but it provides
you with another option to choose from.

Read on to learn how to install and use PyDQED, and tips for getting the most
out of PyDQED.

