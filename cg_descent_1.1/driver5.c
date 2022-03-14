/* Although there is a rigorous theory justifying a Wolfe line search,
   the performance of the Approximate Wolfe line search if often much
   better. Nonetheless, the user can completely turn off the Approximate
   Wolfe line search by setting AWolfe to false and AWolfeFac to 0.
   The output is the following:

   Termination status: 4
   Line search fails, too many secant steps
   Possible causes of this error message:
      - your tolerance may be too strict: grad_tol =  1.000000e-08
   absolute largest component of gradient: 1.590361e-07
   function value: -6.530787e+02
   cg iterations: 29
   function evaluations: 148
   gradient evaluations: 135

   Hence, due to numerical errors, it was not possible to achieve the specified
   1.e-8 error tolerance. By decreasing the error tolerance to 1.e-6 (the first
   argument of cg_descent), the problem is solved:

   Termination status: 0
   Convergence tolerance for gradient satisfied
   absolute largest component of gradient: 7.560665e-07
   function value: -6.530787e+02
   cg iterations: 26
   function evaluations: 49
   gradient evaluations: 33

   On the other hand, if we turn on the Approximate Wolfe line search by
   resetting AWolfe to true and AWolfeFac to its default value 1.e-3, we
   obtain slightly faster convergence:

   Termination status: 0
   Convergence tolerance for gradient satisfied
   absolute largest component of gradient: 6.938785e-07
   function value: -6.530787e+02
   cg iterations: 25
   function evaluations: 48
   gradient evaluations: 31 */

#include "cg_descent.h"
int mydim = 100 ;

double myvalue
(
    double   *x
) ;

void mygrad
(
    double    *g ,
    double    *x
) ;

void main (int argc, char *argv [])
{
    extern int mydim ;
    double *x, *work, step ;
    int i, status ;
    cg_stats Stats ;

    x = (double *) malloc (mydim*sizeof (double)) ;
    work = (double *) malloc (4*mydim*sizeof (double)) ;

/* starting guess */

    for (i = 0; i < mydim; i++) x [i] = 1. ;
    status = cg_descent (1.e-8, x, mydim, myvalue, mygrad, work, step, &Stats) ;
    free (x) ;
    free (work) ;
}

double myvalue
(
    double   *x
)
{
    extern int mydim ;
    double f, t ;
    int i ;
    f = 0. ;
    for (i = 0; i < mydim; i++)
    {
        t = i+1 ;
        t = sqrt (t) ;
        f += exp (x [i]) - t*x [i] ;
    }
    return (f) ;
}

void mygrad
(
    double    *g ,
    double    *x
)
{
    extern int mydim ;
    double t ;
    int i ;
    for (i = 0; i < mydim; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        g [i] = exp (x [i]) -  t ;
    }
    return ;
}
/*
.1        delta        (Wolfe line search parameter)
.9        sigma        (Wolfe line search parameter)
1.e-6     eps          (perturbation parameter for computing fpert)
.66       gamma        (required decay factor in interval)
5.        rho          (interval growth factor used to get bracketing interval)
.01       eta          (lower bound for cg's beta_k)
.01       psi0         (factor used in starting guess for iteration 1)
.1        psi1         (factor previous step multiplied by in QuadStep)
2.        psi2         (factor previous step is multipled by for startup)
1.e-12    QuadCutOff   (QuadStep if relative change in f > QuadCutOff)
0.e-12    StopFact     (factor multiplying starting |grad|_infty in StopRule)
0.e-3     AWolfeFac    (AWolfe = F => set AWolfe = T if |f-f0| < Awolfe_fac*Ck)
1.        restart_fac  (restart cg in restart_fac*n iterations)
500.      maxit_fac    (terminate in maxit_fac*n iterations)
0.        feps         (stop when value change <= feps*|f|)
.7        Qdecay       (used in Qk update: Qk = Qdecay*Qk + 1)
50        nexpand      (number of grow/shrink allowed in bracket)
50        nsecant      (number of secant steps allowed in line search)
1         PertRule     (F => eps, T => eps*Ck)
1         QuadStep     (use initial quad interpolation in line search)
0         PrintLevel   F (no print) T (intermediate results)
1         PrintFinal   F (no print) T (print error messages, final error)
1         StopRule     T (|grad|_infty <= max(tol,|grad|_0*StopFact) F (... <= tol*(1+|f|))
0         AWolfe       F (Wolfe) T (approx Wolfe)
0         Step         F (no initial line search guess) T (guess in step arg)
0         debug        F (no debugging) T (check for no increase in f)
*/
