/* This is the example appearing in the User's Guide. Default values are
   used for the parameters in cg_descent.parm (see listing below).
   The output is the following:

   Termination status: 0
   Convergence tolerance for gradient satisfied
   absolute largest component of gradient: 8.492665e-09
   function value: -6.530787e+02
   cg iterations: 31
   function evaluations: 54
   gradient evaluations: 43 */

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
    double *x, *work, step=0 ;
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
1.e-3     AWolfeFac    (AWolfe = F => set AWolfe = T if |f-f0| < Awolfe_fac*Ck)
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
1         AWolfe       F (Wolfe) T (approx Wolfe)
0         Step         F (no initial line search guess) T (guess in step arg)
0         debug        F (no debugging) T (check for no increase in f)
*/
