/** 
* @file dqed.h
*
* Contains the function prototype for exposing the Fortran bounded constrained
* nonlinear optimization/solver code DQED to C/C++.
*/

#ifdef __cplusplus
extern "C" {
#endif 

/**
 * Typedef for the function that DQED calls to obtains values of the nonlinear 
 * equations and constraints and corresponding analytical Jacobians. Generally,
 * the function uses the provided value of x to update the fj array with the
 * appropriate values.
 * @param x The current value of the unknown variables
 * @param fj An array to store the Jacobian of the constraints and the Jacobian and residual of the equations
 * @param ldfj The length of the leading dimension of fj
 * @param igo The status of the algorithm
 * @param iopt An array of options stored as integers; interpretation varies
 * @param ropt An array of options stored as floating-point numbers; interpretation varies
 */
typedef void (*dqed_function)(double* x, double* fj, int* ldfj, int* igo, int* iopt, double* ropt);

/** 
* Exposes the Fortran bounded constrained nonlinear optimization/solver code 
* DQED to C/C++.
*/
void dqed_(
    dqed_function dqedev,               /** The function that evaluates the nonlinear equations and constraints */
    int* mequa,                         /** The number of equations to be solved */
    int* nvars,                         /** The number of unknown variables */
    int* mcon,                          /** The number of general constraints (not including simple bounds) */
    int* ind,                           /** The type of simple bounds to use for each variable */
    double* bl,                         /** The lower bound to use for each variable (if applicable) */
    double* bu,                         /** The upper bound to use for each variable (if applicable) */
    double* x,                          /** The initial guess (on input) or final solution (on output) */
    double* fj,                         /** The final values of the Jacobian matrix for the constrants and equations */
    int* ldfj,                          /** The leading dimension of the Jacobian matrix */
    double* fnorm,                      /** The value of the Euclidean norm at the solution */
    int* igo,                           /** A flag containing the final status of the optimization/solve */
    int* iopt,                          /** Integer parameters modifying how DQED executes */
    double* ropt,                       /** Double-precision parameters modifying how DQED executes */
    int* iwork,                         /** Work space for integer values */
    double* rwork                       /** Work space for double-precision values */
);

#ifdef __cplusplus
}
#endif 
