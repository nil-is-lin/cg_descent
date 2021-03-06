/* =========================================================================
   ============================ CG_DESCENT =================================
   =========================================================================
       ________________________________________________________________
      |      A conjugate gradient method with guaranteed descent       |
      |             C-code Version 1.1  (October 6, 2005)              |
      |                    Version 1.2  (November 14, 2005)            |
      |                    Version 2.0  (September 23, 2007)           |
      |                    Version 3.0  (May 18, 2008)                 |
      |                    Version 4.0  (March 28, 2011)               |
      |           William W. Hager    and   Hongchao Zhang             |
      |          hager@math.ufl.edu       hzhang@math.ufl.edu          |
      |                   Department of Mathematics                    |
      |                     University of Florida                      |
      |                 Gainesville, Florida 32611 USA                 |
      |                      352-392-0281 x 244                        |
      |                                                                |
      |                 Copyright by William W. Hager                  |
      |                                                                |
      |          http://www.math.ufl.edu/~hager/papers/CG              |
      |________________________________________________________________|
       ________________________________________________________________
      |This program is free software; you can redistribute it and/or   |
      |modify it under the terms of the GNU General Public License as  |
      |published by the Free Software Foundation; either version 2 of  |
      |the License, or (at your option) any later version.             |
      |This program is distributed in the hope that it will be useful, |
      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
      |GNU General Public License for more details.                    |
      |                                                                |
      |You should have received a copy of the GNU General Public       |
      |License along with this program; if not, write to the Free      |
      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
      |MA  02110-1301  USA                                             |
      |________________________________________________________________|

      References:
      1. W. W. Hager and H. Zhang, A new conjugate gradient method
         with guaranteed descent and an efficient line search,
         SIAM Journal on Optimization, 16 (2005), 170-192. 
      2. W. W. Hager and H. Zhang, Algorithm 851: CG_DESCENT,
         A conjugate gradient method with guaranteed descent,
         ACM Transactions on Mathematical Software, 32 (2006), 113-137. 
      3. W. W. Hager and H. Zhang, A survey of nonlinear conjugate gradient
         methods, Pacific Journal of Optimization, 2 (2006), pp. 35-58. */

#include "cg_user.h"
#include "cg_descent.h"
int cg_descent /*  return:
                      -2 (function value became nan)
                      -1 (starting function value is nan)
                       0 (convergence tolerance satisfied)
                       1 (change in func <= feps*|f|)
                       2 (total iterations exceeded maxit)
                       3 (slope always negative in line search)
                       4 (number secant iterations exceed nsecant)
                       5 (search direction not a descent direction)
                       6 (line search fails in initial interval)
                       7 (line search fails during bisection)
                       8 (line search fails during interval update)
                       9 (debugger is on and the function value increases)
                      10 (out of memory) */
(
    double            *x, /* input: starting guess, output: the solution */
    INT                n, /* problem dimension */
    cg_stats       *Stat, /* structure with statistics (can be NULL) */
    cg_parameter  *UParm, /* user parameters, NULL = use default parameters */
    double      grad_tol, /* StopRule = 1: |g|_infty <= max (grad_tol,
                                           StopFac*initial |g|_infty) [default]
                             StopRule = 0: |g|_infty <= grad_tol(1+|f|) */
    double      (*value) (double *, INT),  /* f = value (x, n) */
    void         (*grad) (double *, double *, INT), /* grad (g, x, n) */
    double    (*valgrad) (double *, double *, INT), /* f = valgrad (g, x, n),
                          NULL = compute value & gradient using value & grad */
    double         *Work  /* either size 4n work array or NULL */
)
{
    INT     n5, i, iter, maxit, nrestart, IterRestart ;
    int     IterQuad, status, StopRule ;
    double  delta2, Qk, Ck,
            f, ftemp, gnorm, xnorm, gnorm2, dnorm2, denom,
            t, t1, t2, t3, t4, t5, dphi, dphi0, alpha,
            yk, ykyk, ykgk, dkyk, yk1, yk2, yk3, yk4, yk5, beta, QuadTrust, tol,
           *d, *g, *xtemp, *gtemp, *work ;
    cg_parameter *Parm, ParmStruc ;
    cg_com Com ;

    /* initialize the parameters */
    if ( UParm == NULL )
    {
        Parm = &ParmStruc ;
        cg_default (Parm) ;
    }
    else Parm = UParm ;
    Com.Parm = Parm ;

    if ( Parm->PrintParms ) cg_printParms (Parm) ;

    /* allocate work arrays */
    if ( Work == NULL ) work = malloc (4*n*sizeof (double)) ;
    else                work = Work ;
    if ( work == NULL )
    {
        printf ("Insufficient memory for specified problem dimension %e\n",
                 (double) n) ;
        status = 10 ;
        return (status) ;
    }
    Com.x = x ;
    Com.d = d = work ;
    Com.g = g = d+n ;
    Com.xtemp = xtemp = g+n ;
    Com.gtemp = gtemp = xtemp+n ;
    Com.n = n ;          /* problem dimension */
    Com.nf = (INT) 0 ;   /* number of function evaluations */
    Com.ng = (INT) 0 ;   /* number of gradient evaluations */
    Com.AWolfe = Parm->AWolfe ; /* do not touch user's AWolfe */
    Com.cg_value = value ;
    Com.cg_grad = grad ;
    Com.cg_valgrad = valgrad ;
    StopRule = Parm->StopRule ;

    /* the conjugate gradient algorithm is restarted every nrestart iteration */
    nrestart = (INT) (((double) n)*Parm->restart_fac) ;

    /* abort when number of iterations reaches maxit */
    if ( Parm->maxit_fac == INF ) maxit = INT_INF ;
    else                          maxit = (INT) (((double) n)*Parm->maxit_fac) ;
    
    f = ZERO ;
    n5 = n % 5 ;

    Ck = ZERO ;
    Qk = ZERO ;

    /* initial function and gradient evaluations, initial direction */
    Com.alpha = ZERO ;
    cg_evaluate ("fg", "n", &Com) ;
    f = Com.f ;
    Com.f0 = f + f ;
    xnorm = ZERO ;
    for (i = 0; i < n5; i++) if ( xnorm < fabs (x [i]) ) xnorm = fabs (x [i]) ;
    for (; i < n; i += 5)
    {
         if ( xnorm < fabs (x [i]  ) ) xnorm = fabs (x [i]  ) ;
         if ( xnorm < fabs (x [i+1]) ) xnorm = fabs (x [i+1]) ;
         if ( xnorm < fabs (x [i+2]) ) xnorm = fabs (x [i+2]) ;
         if ( xnorm < fabs (x [i+3]) ) xnorm = fabs (x [i+3]) ;
         if ( xnorm < fabs (x [i+4]) ) xnorm = fabs (x [i+4]) ;
    }
    gnorm = ZERO ;
    gnorm2 = ZERO ;
    for (i = 0; i < n5; i++)
    {
        t = g [i] ;
        d [i] = -t ;
        gnorm2 += t*t ;
        if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
    }
    for (; i < n;)
    {
        t1 = g [i] ;
        d [i] = -t1 ;
        if ( gnorm < fabs (t1) ) gnorm = fabs (t1) ;
        i++ ;

        t2 = g [i] ;
        d [i] = -t2 ;
        if ( gnorm < fabs (t2) ) gnorm = fabs (t2) ;
        i++ ;

        t3 = g [i] ;
        d [i] = -t3 ;
        if ( gnorm < fabs (t3) ) gnorm = fabs (t3) ;
        i++ ;

        t4 = g [i] ;
        d [i] = -t4 ;
        if ( gnorm < fabs (t4) ) gnorm = fabs (t4) ;
        i++ ;

        t5 = g [i] ;
        d [i] = -t5 ;
        if ( gnorm < fabs (t5) ) gnorm = fabs (t5) ;
        i++ ;

        gnorm2 += t1*t1 + t2*t2 + t3*t3 + t4*t4 + t5*t5 ;
    }
    /* check if the starting function value is nan */
    if ( f != f )
    {
        status = -1 ;
        goto Exit ;
    }

    if ( Parm->StopRule ) tol = MAX (gnorm*Parm->StopFac, grad_tol) ;
    else                  tol = grad_tol ;

    if ( Parm->PrintLevel >= 1 )
    {
        printf ("iter: %5i f = %14.6e gnorm = %14.6e AWolfe = %2i\n",
          (int) 0, f, gnorm, Com.AWolfe) ;
    }

    if ( cg_tol (f, gnorm, StopRule, tol) )
    {
        iter = 0 ;
        status = 0 ;
        goto Exit ;
    }

    dphi0 = -gnorm2 ;
    delta2 = 2*Parm->delta - ONE ;
    alpha = Parm->step ;
    if ( alpha == 0. )
    {
        alpha = Parm->psi0*xnorm/gnorm ;
        if ( xnorm == ZERO )
        {
            if ( f != ZERO ) alpha = Parm->psi0*fabs (f)/gnorm2 ;
            else             alpha = ONE ;
        }
    }
    IterRestart = 0 ;  /* counts number of iterations since last restart */
    IterQuad = 0 ;     /* counts number of iterations that function change
                          is close to that of a quadratic */
 
    /* Start the conjugate gradient iteration.
       alpha starts as old step, ends as final step for current iteration
       f is function value for alpha = 0
       QuadOK = TRUE means that a quadratic step was taken */
 
    for (iter = 1; iter <= maxit; iter++)
    {
        Com.QuadOK = FALSE ;
        alpha = Parm->psi2*alpha ;
        if ( Parm->QuadStep )
        {
            if ( f != ZERO ) t = fabs ((f-Com.f0)/f) ;
            else             t = ONE ;
            /* test if quadratic interpolation step should be tried */
            if ( t > Parm->QuadCutOff )
            {
                Com.alpha = Parm->psi1*alpha ;
                status = cg_evaluate ("f", "y", &Com) ;
                if ( status ) goto Exit ;
                ftemp = Com.f ;

                if ( ftemp < f )              /* check if quadstep > 0 */
                {
                   denom = 2.*(((ftemp-f)/Com.alpha)-dphi0) ;
                   if ( denom > ZERO )        /* try a quadratic fit step */
                   {
                       Com.QuadOK = TRUE ;
                       alpha = -dphi0*Com.alpha/denom ;
                   }
                }
            }
        }
        Com.f0 = f ;                          /* f0 saved as prior value */

        if ( Parm->PrintLevel >= 1 )
        {
            printf ("QuadOK: %2i initial a: %14.6e f0: %14.6e dphi: %14.6e\n",
                    Com.QuadOK, alpha, Com.f0, dphi0) ;
        }

        /* parameters in Wolfe and approximate Wolfe conditions, and in update*/

        Qk = Parm->Qdecay*Qk + ONE ;
        Ck = Ck + (fabs (f) - Ck)/Qk ;        /* average cost magnitude */

        if ( Parm->PertRule ) Com.fpert = f + Parm->eps*Ck ;
        else                  Com.fpert = f + Parm->eps ;

        Com.wolfe_hi = Parm->delta*dphi0 ;
        Com.wolfe_lo = Parm->sigma*dphi0 ;
        Com.awolfe_hi = delta2*dphi0 ;
        Com.alpha = alpha ;        /* either prior step or quadratic fit step */
        Com.f = f ;

        if ( Com.AWolfe ) status = cg_line (dphi0, &Com) ; /* approx. Wolfe */
        else              status = cg_lineW (dphi0, &Com) ;/* ordinary Wolfe */
        if ( (status > 0) && !Com.AWolfe )/*try approximate Wolfe line search*/
        {
            if ( Parm->PrintLevel >= 1 )
            {
                 printf ("\nWOLFE LINE SEARCH FAILS\n") ;
            }
            Com.AWolfe = TRUE ;
            status = cg_line (dphi0, &Com) ;
        }

        alpha = Com.alpha ;
        f = Com.f ;
        dphi = Com.df ;

        if ( status ) goto Exit ;

        /* Test for convergence to within machine epsilon
           [set feps to zero to remove this test] */
 
        if ( -alpha*dphi0 <= Parm->feps*fabs (f) )
        {
            status = 1 ;
            goto Exit ;
        }

        /* test how close the cost function changes are to that of a quadratic
           QuadTrust = 0 means the function change matches that of a quadratic*/
        t = alpha*(dphi+dphi0) ;
        if ( fabs (t) <= Parm->qeps*MIN (Ck, ONE) ) QuadTrust = ZERO ;
        else QuadTrust = fabs((2.0*(f-Com.f0)/t)-ONE) ;
        if ( QuadTrust <= Parm->qrule) IterQuad++ ;
        else                           IterQuad = 0 ;

        IterRestart++ ;
        /* test if the CG algorithm should be restarted */
        if ( (IterRestart == nrestart) ||
             ((IterQuad == Parm->qrestart) && (IterQuad != IterRestart)) )
        {
            IterRestart = 0 ;
            IterQuad = 0 ;
            /* search direction d = -g */
            if ( Parm->PrintLevel >= 1 ) printf ("RESTART CG\n") ;
            gnorm = ZERO ;
            gnorm2 = ZERO ;
            cg_copy (x, xtemp, n) ;
            for (i = 0; i < n5; i++)
            {
                t = gtemp [i] ;
                g [i] = t ;
                d [i] = -t ;
                if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
                gnorm2 += t*t ;
            }
            for (; i < n; )
            {
                t1 = gtemp [i] ;
                g [i] = t1 ;
                d [i] = -t1 ;
                if ( gnorm < fabs (t1) ) gnorm = fabs (t1) ;
                i++ ;

                t2 = gtemp [i] ;
                g [i] = t2 ;
                d [i] = -t2 ;
                if ( gnorm < fabs (t2) ) gnorm = fabs (t2) ;
                i++ ;

                t3 = gtemp [i] ;
                g [i] = t3 ;
                d [i] = -t3 ;
                if ( gnorm < fabs (t3) ) gnorm = fabs (t3) ;
                i++ ;

                t4 = gtemp [i] ;
                g [i] = t4 ;
                d [i] = -t4 ;
                if ( gnorm < fabs (t4) ) gnorm = fabs (t4) ;
                i++ ;

                t5 = gtemp [i] ;
                g [i] = t5 ;
                d [i] = -t5 ;
                if ( gnorm < fabs (t5) ) gnorm = fabs (t5) ;
                i++ ;
                gnorm2 += t1*t1 + t2*t2 + t3*t3 + t4*t4 + t5*t5 ;
            }
            if ( cg_tol (f, gnorm, StopRule, tol) )
            {
                status = 0 ;
                goto Exit ;
            }
            dphi0 = -gnorm2 ;
        }
        else /* compute beta, yk2, gnorm, gnorm2, dnorm2, update x and g */
        {
            cg_copy (x, xtemp, n) ;
            dnorm2 = ZERO ;
            for (i = 0; i < n5; i++) dnorm2 = dnorm2 + d [i]*d [i] ;
            for (; i < n; i += 5)
            {
                dnorm2 = dnorm2 + d [i]*d [i] + d [i+1]*d [i+1]
                                              + d [i+2]*d [i+2]
                                              + d [i+3]*d [i+3]
                                              + d [i+4]*d [i+4] ;
            }
            gnorm = ZERO ;
            ykyk = ZERO ;
            ykgk = ZERO ;
            for (i = 0; i < n5; i++)
            {
                t = gtemp [i] ;
                if ( gnorm < fabs (t) ) gnorm = fabs (t) ;
                yk = t - g [i] ;
                g [i] = t ;
                ykgk += yk*t ;
                ykyk += yk*yk ;
            }
            for (; i < n; )
            {
                t1 = gtemp [i] ;
                yk1 = t1 - g [i] ;
                g [i] = t1 ;
                if ( gnorm < fabs (t1) ) gnorm = fabs (t1) ;
                i++ ;

                t2 = gtemp [i] ;
                yk2 = t2 - g [i] ;
                g [i] = t2 ;
                if ( gnorm < fabs (t2) ) gnorm = fabs (t2) ;
                i++ ;

                t3 = gtemp [i] ;
                yk3 = t3 - g [i] ;
                g [i] = t3 ;
                if ( gnorm < fabs (t3) ) gnorm = fabs (t3) ;
                i++ ;

                t4 = gtemp [i] ;
                yk4 = t4 - g [i] ;
                g [i] = t4 ;
                if ( gnorm < fabs (t4) ) gnorm = fabs (t4) ;
                i++ ;

                t5 = gtemp [i] ;
                yk5 = t5 - g [i] ;
                g [i] = t5 ;
                if ( gnorm < fabs (t5) ) gnorm = fabs (t5) ;

                i++ ;
                ykyk += yk1*yk1 + yk2*yk2 + yk3*yk3 + yk4*yk4 + yk5*yk5 ;
                ykgk += yk1*t1  + yk2*t2  + yk3*t3  + yk4*t4  + yk5*t5 ;
            }

            if ( cg_tol (f, gnorm, StopRule, tol) )
            {
                status = 0 ;
                goto Exit ;
            }
            dkyk = dphi - dphi0 ;
            if ( Parm->AdaptiveBeta ) t = 2. - ONE/(0.1*QuadTrust + ONE) ;
            else                      t = Parm->theta ;
            beta = (ykgk - t*dphi*ykyk/dkyk)/dkyk ;

/*  faster: initialize dnorm2 = gnorm2 at start, then
            dnorm2 = gnorm2 + beta**2*dnorm2 - 2.*beta*dphi
            gnorm2 = ||g_{k+1}||^2
            dnorm2 = ||d_{k+1}||^2
            dpi = g_{k+1}' d_k */

/* lower bound for beta is BetaLower*d_k'g_k/ ||d_k||^2 */
            beta = MAX (beta, Parm->BetaLower*dphi0/dnorm2) ;

/*    update search direction d = -g + beta*dold */

            gnorm2 = ZERO ;
            for (i = 0; i < n5; i++)
            {
                t = g [i] ;
                d [i] = -t + beta*d [i] ;
                gnorm2 += t*t ;
            }
            for (; i < n; )
            {
                t1 = g [i] ;
                d [i] = -t1 + beta*d [i] ;
                i++ ;

                t2 = g [i] ;
                d [i] = -t2 + beta*d [i] ;
                i++ ;

                t3 = g [i] ;
                d [i] = -t3 + beta*d [i] ;
                i++ ;

                t4 = g [i] ;
                d [i] = -t4 + beta*d [i] ;
                i++ ;

                t5 = g [i] ;
                d [i] = -t5 + beta*d [i] ;
                i++ ;

                gnorm2 += t1*t1 + t2*t2 + t3*t3 + t4*t4 + t5*t5 ;
            }
            dphi0 = -gnorm2 + beta*dphi ;
            if ( Parm->debug ) /* Check the dphi0 = d'g */
            {
                t = ZERO ;
                for (i = 0; i < n; i++)  t = t + d [i]*g [i] ;
                if ( fabs(t-dphi0) > Parm->debugtol*fabs(dphi0) )
                {
                    printf("Warning, dphi0 != d'g!\n");
                    printf("dphi0:%14.6e, d'g:%14.6e\n",dphi0, t) ;
                }
            }
        }
        if ( !Com.AWolfe )
        {
            if ( fabs (f-Com.f0) < Parm->AWolfeFac*Ck ) Com.AWolfe = TRUE ;
        }
    
        if ( Parm->PrintLevel >= 1 )
        {
            printf ("\niter: %5i f = %14.6e gnorm = %14.6e AWolfe = %2i\n",
               (int) iter, f, gnorm, Com.AWolfe) ;
        }

        if ( Parm->debug )
        {
            if ( f > Com.f0 + Parm->debugtol*Ck )
            {
                status = 9 ;
                goto Exit ;
            }
        }
                
        if ( dphi0 > ZERO )
        {
           status = 5 ;
           goto Exit ;
        }
    }
    status = 2 ;

Exit:
    if ( Stat != NULL )
    {
        Stat->f = f ;
        Stat->gnorm = gnorm ;
        Stat->nfunc = Com.nf ;
        Stat->ngrad = Com.ng ;
        Stat->iter = iter ;
    }
    if ( status > 2 )
    {
        gnorm = ZERO ;
        for (i = 0; i < n; i++)
        {
            x [i] = xtemp [i] ;
            g [i] = gtemp [i] ;
            t = fabs (g [i]) ;
            gnorm = MAX (gnorm, t) ;
        }
        if ( Stat != NULL ) Stat->gnorm = gnorm ;
    }
    if ( Parm->PrintFinal || Parm->PrintLevel >= 1 )
    {
        const char mess1 [] = "Possible causes of this error message:" ;
        const char mess2 [] = "   - your tolerance may be too strict: "
                              "grad_tol = " ;
        const char mess3 [] = "Line search fails" ;
        const char mess4 [] = "   - your gradient routine has an error" ;
        const char mess5 [] = "   - the parameter epsilon in cg_descent_c.parm "
                              "is too small" ;
        printf ("\nTermination status: %i\n", status) ;
        if ( status == -2 )
        {
            printf ("At iteration %10.0f function value became nan\n",
                    (double) iter) ;
        }
        else if ( status == -1 )
        {
            printf ("Objective function value is nan at starting point\n") ;
        }
        else if ( status == 0 )
        {
            printf ("Convergence tolerance for gradient satisfied\n") ;
        }
        else if ( status == 1 )
        {
            printf ("Terminating since change in function value "
                    "<= feps*|f|\n") ;
        }
        else if ( status == 2 )
        {
            printf ("Number of iterations exceed specified limit\n") ;
            printf ("Iterations: %10.0f maxit: %10.0f\n",
                    (double) iter, (double) maxit) ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, grad_tol) ;
        }
        else if ( status == 3 )
        {
            printf ("Slope always negative in line search\n") ;
            printf ("%s\n", mess1) ;
            printf ("   - your cost function has an error\n") ;
            printf ("%s\n", mess4) ;
        }
        else if ( status == 4 )
        {
            printf ("Line search fails, too many secant steps\n") ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, grad_tol) ;
        }
        else if ( status == 5 )
        {
            printf ("Search direction not a descent direction\n") ;
        }
        else if ( status == 6 ) /* line search fails */
        {
            printf ("%s\n", mess3) ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, grad_tol) ;
            printf ("%s\n", mess4) ;
            printf ("%s\n", mess5) ;
        }
        else if ( status == 7 ) /* line search fails */
        {
            printf ("%s\n", mess3) ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, grad_tol) ;
        }
        else if ( status == 8 ) /* line search fails */
        {
            printf ("%s\n", mess3) ;
            printf ("%s\n", mess1) ;
            printf ("%s %e\n", mess2, grad_tol) ;
            printf ("%s\n", mess4) ;
            printf ("%s\n", mess5) ;
        }
        else if ( status == 9 )
        {
            printf ("Debugger is on, function value does not improve\n") ;
            printf ("new value: %25.16e old value: %25.16e\n", f, Com.f0) ;
        }
        else if ( status == 10 )
        {
            printf ("Insufficient memory\n") ;
        }

        printf ("maximum norm for gradient: %13.6e\n", gnorm) ;
        printf ("function value:            %13.6e\n\n", f) ;
        printf ("cg  iterations:          %10.0f\n", (double) iter) ;
        printf ("function evaluations:    %10.0f\n", (double) Com.nf) ;
        printf ("gradient evaluations:    %10.0f\n", (double) Com.ng) ;
        printf ("===================================\n\n") ;
    }
    if ( Work == NULL ) free (work) ;
    return (status) ;
}

/* =========================================================================
   ==== cg_Wolfe ===========================================================
   =========================================================================
   Check whether the Wolfe or the approximate Wolfe conditions are satisfied
   ========================================================================= */
PRIVATE int cg_Wolfe
(
    double   alpha, /* stepsize */
    double       f, /* function value associated with stepsize alpha */
    double    dphi, /* derivative value associated with stepsize alpha */
    cg_com    *Com  /* cg com */
)
{
    if ( dphi >= Com->wolfe_lo )
    {

/* test original Wolfe conditions */

        if ( f - Com->f0 <= alpha*Com->wolfe_hi )
        {
            if ( Com->Parm->PrintLevel >= 2 )
            {
                printf ("wolfe f: %14.6e f0: %14.6e dphi: %14.6e\n",
                         f, Com->f0, dphi) ;
            }
            return (1) ;
        }
/* test approximate Wolfe conditions */
        else if ( Com->AWolfe )
        {
            if ( (f <= Com->fpert) && (dphi <= Com->awolfe_hi) )
            {
                if ( Com->Parm->PrintLevel >= 2 )
                {
                    printf ("f: %14.6e fpert: %14.6e dphi: %14.6e awolf_hi: "
                            "%14.6e\n", f, Com->fpert, dphi, Com->awolfe_hi) ;
                }
                return (1) ;
            }
        }
    }
    return (0) ;
}

/* =========================================================================
   ==== cg_tol =============================================================
   =========================================================================
   Check for convergence
   ========================================================================= */
PRIVATE int cg_tol
(
    double         f, /* function value associated with stepsize */
    double     gnorm, /* gradient sup-norm */
    int     StopRule, /* T => |grad|_infty <=max (tol, |grad|_infty*StopFact)
                         F => |grad|_infty <= tol*(1+|f|)) */
    double       tol  /* tolerance */
)
{
    if ( StopRule )
    {
        if ( gnorm <= tol ) return (1) ;
    }
    else if ( gnorm <= tol*(ONE + fabs (f)) ) return (1) ;
    return (0) ;
}

/* =========================================================================
   ==== cg_line ============================================================
   =========================================================================
   Approximate Wolfe line search routine
   Return:
      -2 (function nan)
       0 (convergence tolerance satisfied)
       3 (slope always negative in line search)
       4 (number secant iterations exceed nsecant)
       6 (line search fails in initial interval)
       7 (line search fails during bisection)
       8 (line search fails during interval update)
   ========================================================================= */
PRIVATE int cg_line
(
    double  dphi0, /* function derivative at starting point (alpha = 0) */
    cg_com   *Com  /* cg com structure */
)
{
    int iter, nshrink, ngrow, PrintLevel, status ;
    double a, dphia, b, dphib, c, alpha, phi, dphi,
           a0, da0, b0, db0, width, fquad, rho ;
    cg_parameter *Parm ;

    Parm = Com->Parm ;
    PrintLevel = Parm->PrintLevel ;
    if ( PrintLevel >= 1 ) printf ("Approximate Wolfe line search\n") ;
    status = cg_evaluate ("g", "y", Com) ;
    if ( status ) return (status) ; /* return if function is nan */
    alpha = Com->alpha ;
    dphi = Com->df ;
    rho = Com->rho ;
 
/*Find initial interval [a,b] such that dphia < 0, dphib >= 0,
         and phia <= phi0 + feps*fabs (phi0) */
 
    a = ZERO ;
    dphia = dphi0  ;
    ngrow = 0 ;
    nshrink = 0 ;
    while ( dphi < ZERO )
    {
        cg_evaluate ("f", "n", Com) ;
        phi = Com->f ;

/* if quadstep in effect and quadratic conditions hold, check wolfe condition*/

        if ( Com->QuadOK )
        {
            if ( ngrow == 0 ) fquad = MIN (phi, Com->f0) ;
            if ( phi <= fquad )
            {
                if ( PrintLevel >= 2 )
                {
                    printf ("alpha: %14.6e phi: %14.6e fquad: %14.6e\n",
                            alpha, phi, fquad) ;
                }
                if ( cg_Wolfe (alpha, phi, dphi, Com) ) return (0) ;
            }
        }
        if ( phi > Com->fpert )
        {
            /* contraction phase, only break at termination or Secant step */
            b = alpha ;
            while ( TRUE )
            {
                alpha = .5*(a+b) ;
                Com->alpha = alpha ;
                nshrink++ ;
                if ( nshrink > Parm->nexpand ) return (6) ;
                cg_evaluate ("g", "n", Com) ;
                dphi = Com->df ;
                if ( dphi >= ZERO ) goto Secant ;
                cg_evaluate ("f", "n", Com) ;
                phi = Com->f ;
                if ( PrintLevel >= 2 )
                {
                    printf ("contract, a: %14.6e b: %14.6e alpha: %14.6e phi: "
                            "%14.6e dphi: %14.6e\n", a, b, alpha, phi, dphi) ;
                }
                if ( Com->QuadOK && (phi <= fquad) )
                {
                    if ( cg_Wolfe (alpha, phi, dphi, Com) ) return (0) ;
                }
                if ( phi <= Com->fpert )
                {
                    a = alpha ;
                    dphia = dphi ;
                }
                else
                {
                    b = alpha ;
                }
            }
        }

/* expansion phase */

        a = alpha ;
        dphia = dphi ;
        ngrow++ ;
        if ( ngrow > Parm->nexpand ) return (3) ;
        alpha = rho*alpha ;
        Com->alpha = alpha ;
        cg_evaluate ("g", "n", Com) ;
        dphi = Com->df ;
        if ( PrintLevel >= 2 )
        {
            printf ("expand,   a: %14.6e alpha: %14.6e phi: "
                     "%14.6e dphi: %14.6e\n", a, alpha, phi, dphi) ;
        }
    }

Secant:
    b = alpha ;
    dphib = dphi ;
    if ( Com->QuadOK )
    {
        cg_evaluate ("f", "n", Com) ;
        phi = Com->f ;
        if ( ngrow + nshrink == 0 ) fquad = MIN (phi, Com->f0) ;
        if ( phi <= fquad )
        {
            if ( cg_Wolfe (alpha, phi, dphi, Com) ) return (0) ;
        }
    }
    for (iter = 1; iter <= Parm->nsecant; iter++)
    {
        if ( PrintLevel >= 2 )
        {
            printf ("secant, a: %14.6e b: %14.6e da: %14.6e db: %14.6e\n",
                     a, b, dphia, dphib) ;
        }
        width = Parm->gamma*(b - a) ;
        if ( -dphia <= dphib ) alpha = a - (a-b)*(dphia/(dphia-dphib)) ;
        else                   alpha = b - (a-b)*(dphib/(dphia-dphib)) ;
        Com->alpha = alpha ;
        a0 = a ;
        b0 = b ;
        da0 = dphia ;
        db0 = dphib ;
        status = cg_update (&a, &dphia, &b, &dphib, Com) ;
        if ( status >= 0 ) return (status) ;
        c = Com->alpha ;
        dphi = Com->df ;
        if ( status == -2 )
        {
            if ( c == a )
            {
                if ( dphi > da0 ) alpha = c - (c-a0)*(dphi/(dphi-da0)) ;
                else              alpha = a ;
            }
            else
            {
                if ( dphi < db0 ) alpha = c - (c-b0)*(dphi/(dphi-db0)) ;
                else              alpha = b ;
            }
            Com->alpha = alpha ;
            if ( (alpha > a) && (alpha < b) )
            {
                if ( PrintLevel >= 2 ) printf ("2nd secant\n") ;
                status = cg_update (&a, &dphia, &b, &dphib, Com) ;
                dphi = Com->df ;
                if ( status >= 0 ) return (status) ;
            }
        }

/* bisection iteration */

        if ( b-a >= width )
        {
            alpha = .5*(b+a) ;
            Com->alpha = alpha ;
            if ( PrintLevel >= 2 ) printf ("bisection\n") ;
            status = cg_update (&a, &dphia, &b, &dphib, Com) ;
            dphi = Com->df ;
            if ( status >= 0 ) return (status) ;
        }
        else if ( b <= a ) return (7) ;
    }
    return (4) ;
}

/* =========================================================================
   ==== cg_update ==========================================================
   =========================================================================
   update returns: 8 if too many iterations
                   0 if Wolfe condition is satisfied
                  -1 if interval is updated and a search is done
                  -2 if the interval updated successfully
   ========================================================================= */
PRIVATE int cg_update
(
    double        *a , /* left side of bracketing interval */
    double    *dphia , /* derivative at a */
    double        *b , /* right side of bracketing interval */
    double    *dphib , /* derivative at b */
    cg_com      *Com   /* cg com structure */
)
{
    int nshrink, status ;
    cg_parameter *Parm ;

    Parm = Com->Parm ;
    cg_evaluate ("fg", "n", Com) ;
    if ( Parm->PrintLevel >= 2 )
    {
        printf ("update alpha: %14.6e phi: %14.6e dphi: %14.6e\n",
                 Com->alpha, Com->f, Com->df) ;
    }
    if ( cg_Wolfe (Com->alpha, Com->f, Com->df, Com) )
    {
        status = 0 ;
        goto Exit ;
    }
    status = -2 ;
    if ( Com->df >= ZERO )
    {
        *b = Com->alpha ;
        *dphib = Com->df ;
        goto Exit ;
    }
    else if ( Com->f <= Com->fpert )
    {
        *a = Com->alpha ;
        *dphia = Com->df ;
        goto Exit ;
    }
    nshrink = 0 ;
    *b = Com->alpha ;
    while ( TRUE )
    {
        Com->alpha = .5*(*a + *b) ;
        nshrink++ ;
        if ( nshrink > Parm->nexpand )
        {
            status = 8 ;
            goto Exit ;
        }
        cg_evaluate ("fg", "n", Com) ;
        if ( Parm->PrintLevel >= 2 )
        {
            printf ("contract, a: %14.6e alpha: %14.6e phi: "
                    "%14.6e dphi: %14.6e\n", *a, Com->alpha, Com->f, Com->df) ;
        }
        if ( cg_Wolfe (Com->alpha, Com->f, Com->df, Com) )
        {
            status = 0 ;
            goto Exit ;
        }
        if ( Com->df >= ZERO )
        {
            *b = Com->alpha ;
            *dphib = Com->df ;
            status = -1 ;
            goto Exit ;
        }
        if ( Com->f <= Com->fpert )
        {
            if ( Parm->PrintLevel >= 2 )
            {
                printf ("update a: %14.6e dphia: %14.6e\n",
                         Com->alpha, Com->df) ;
            }
            *a = Com->alpha ;
            *dphia = Com->df ;
        }
        else *b = Com->alpha ;
    }
Exit:
    if ( Parm->PrintLevel >= 2 )
    {
        printf ("UP a: %14.6e b: %14.6e da: %14.6e db: %14.6e status: %i\n",
                 *a, *b, *dphia, *dphib, status) ;
    }
    return (status) ;
}

/* =========================================================================
   ==== cg_lineW ===========================================================
   =========================================================================
   Ordinary Wolfe line search routine.
   This routine is identical to cg_line except that the function
   psi [a] = phi [a] - phi [0] - a*delta*dphi [0] is minimized instead of
   the function phi
   Return:
      -2 (function nan)
       0 (convergence tolerance satisfied)
       3 (slope always negative in line search)
       4 (number secant iterations exceed nsecant)
       6 (line search fails in initial interval)
       7 (line search fails during bisection)
       8 (line search fails during interval update)
   ========================================================================= */
PRIVATE int cg_lineW
(
    double  dphi0, /* function derivative at starting point (alpha = 0) */
    cg_com   *Com  /* cg com structure */
)
{
    int iter, nshrink, ngrow, PrintLevel, status ;
    double a, dpsia, b, dpsib, c, alpha, phi, dphi,
           a0, da0, b0, db0, width, fquad, rho, psi, dpsi ;
    cg_parameter *Parm ;

    Parm = Com->Parm ;
    PrintLevel = Parm->PrintLevel ;
    if ( PrintLevel >= 1 ) printf ("Wolfe line search\n") ;
    status = cg_evaluate ("g", "y", Com) ;
    if ( status ) return (status) ;
    alpha = Com->alpha ;
    dphi = Com->df ;
    rho = Com->rho ;
    dpsi = dphi - Com->wolfe_hi ;
 
/*Find initial interval [a,b] such that dphia < 0, dphib >= 0,
         and phia <= phi0 + feps*fabs (phi0) */
 
    a = ZERO ;
    dpsia = dphi0 - Com->wolfe_hi ;
    ngrow = 0 ;
    nshrink = 0 ;
    while ( dpsi < ZERO )
    {
        cg_evaluate ("f", "n", Com) ;
        phi = Com->f ;
        psi = phi - alpha*Com->wolfe_hi ;

/* if quadstep in effect and quadratic conditions hold, check Wolfe condition*/

        if ( Com->QuadOK )
        {
            if ( ngrow == 0 ) fquad = MIN (phi, Com->f0) ;
            if ( phi <= fquad )
            {
                if ( PrintLevel >= 2 )
                {
                    printf ("alpha: %14.6e phi: %14.6e fquad: %14.6e\n",
                            alpha, phi, fquad) ;
                }
                if ( cg_Wolfe (alpha, phi, dphi, Com) ) return (0) ;
            }
        }
        if ( psi <= Com->fpert )
        {
            a = alpha ;
            dpsia = dphi ;
        }
        else
        {
            /* contraction phase, only break at termination or Secant step */
            b = alpha ;
            while ( TRUE )
            {
                alpha = .5*(a+b) ;
                Com->alpha = alpha ;
                nshrink++ ;
                if ( nshrink > Parm->nexpand ) return (6) ;
                cg_evaluate ("g", "n", Com) ;
                dphi = Com->df ;
                dpsi = dphi - Com->wolfe_hi ;
                if ( dpsi >= ZERO ) goto Secant ;
                cg_evaluate ("f", "n", Com) ;
                phi = Com->f ;
                psi = phi - alpha*Com->wolfe_hi ;
                if ( PrintLevel >= 2 )
                {
                    printf ("contract, a: %14.6e b: %14.6e alpha: %14.6e phi: "
                            "%14.6e dphi: %14.6e\n", a, b, alpha, phi, dphi) ;
                }
                if ( Com->QuadOK && (phi <= fquad) )
                {
                    if ( cg_Wolfe (alpha, phi, dphi, Com) ) return (0) ;
                }
                if ( psi <= Com->fpert )
                {
                    a = alpha ;
                    dpsia = dpsi ;
                }
                else
                {
                    b = alpha ;
                }
            }
        }

/* expansion phase */

        ngrow++ ;
        if ( ngrow > Parm->nexpand ) return (3) ;
        alpha *= rho ;
        Com->alpha = alpha ;
        cg_evaluate ("g", "n", Com) ;
        dphi = Com->df ;
        dpsi = dphi - Com->wolfe_hi ;
        if ( PrintLevel >= 2 )
        {
            printf ("expand,   a: %14.6e alpha: %14.6e phi: "
                     "%14.6e dphi: %14.6e\n", a, alpha, phi, dphi) ;
        }
    }

Secant:
    b = alpha ;
    dpsib = dpsi ;
    if ( Com->QuadOK )
    {
        cg_evaluate ("f", "n", Com) ;
        phi = Com->f ;
        if ( ngrow + nshrink == 0 ) fquad = MIN (phi, Com->f0) ;
        if ( phi <= fquad )
        {
            if ( cg_Wolfe (alpha, phi, dphi, Com) ) return (0) ;
        }
    }
    for (iter = 1; iter <= Parm->nsecant; iter++)
    {
        if ( PrintLevel >= 2 )
        {
            printf ("secant, a: %14.6e b: %14.6e da: %14.6e db: %14.6e\n",
                     a, b, dpsia, dpsib) ;
        }
        width = Parm->gamma*(b - a) ;
        if ( -dpsia <= dpsib ) alpha = a - (a-b)*(dpsia/(dpsia-dpsib)) ;
        else                   alpha = b - (a-b)*(dpsib/(dpsia-dpsib)) ;
        Com->alpha = alpha ;
        a0 = a ;
        b0 = b ;
        da0 = dpsia ;
        db0 = dpsib ;
        status = cg_updateW (&a, &dpsia, &b, &dpsib, &dpsi, Com) ;
        if ( status >= 0 ) return (status) ;
        c = Com->alpha ;
        if ( status == -2 )
        {
            if ( c == a )
            {
                if ( dpsi > da0 ) alpha = c - (c-a0)*(dpsi/(dpsi-da0)) ;
                else              alpha = a ;
            }
            else
            {
                if ( dpsi < db0 ) alpha = c - (c-b0)*(dpsi/(dpsi-db0)) ;
                else              alpha = b ;
            }
            Com->alpha = alpha ;
            if ( (alpha > a) && (alpha < b) )
            {
                if ( PrintLevel >= 2 ) printf ("2nd secant\n") ;
                status = cg_updateW (&a, &dpsia, &b, &dpsib, &dpsi, Com) ;
                if ( status >= 0 ) return (status) ;
            }
        }

/* bisection iteration */

        if ( b-a >= width )
        {
            alpha = .5*(b+a) ;
            Com->alpha = alpha ;
            if ( PrintLevel >= 2 ) printf ("bisection\n") ;
            status = cg_updateW (&a, &dpsia, &b, &dpsib, &dpsi, Com) ;
            if ( status >= 0 ) return (status) ;
        }
        else if ( b <= a ) return (7) ;
    }
    return (4) ;
}

/* =========================================================================
   ==== cg_updateW =========================================================
   =========================================================================
   This routine is identical to cg_update except that the function
   psi [a] = phi [a] - phi [0] - a*delta*dphi [0] is minimized instead of
   the function phi. The return int has the following meaning:
                   8 if too many iterations
                   0 if Wolfe condition is satisfied
                  -1 if interval is updated and a search is done
                  -2 if the interval updated successfully
   ========================================================================= */
PRIVATE int cg_updateW
(
    double        *a , /* left side of bracketing interval */
    double    *dpsia , /* derivative at a */
    double        *b , /* right side of bracketing interval */
    double    *dpsib , /* derivative at b */
    double     *dpsi , /* derivative of psi at alpha (returned) */
    cg_com      *Com   /* cg com structure */
)
{
    int nshrink, status ;
    double psi ;
    cg_parameter *Parm ;

    Parm = Com->Parm ;
    cg_evaluate ("fg", "n", Com) ;
    psi = Com->f - Com->alpha*Com->wolfe_hi ;
    *dpsi = Com->df - Com->wolfe_hi ;
    if ( Parm->PrintLevel >= 2 )
    {
        printf ("update alpha: %14.6e psi: %14.6e dpsi: %14.6e\n",
                 Com->alpha, psi, *dpsi) ;
    }
    if ( cg_Wolfe (Com->alpha, Com->f, Com->df, Com) )
    {
        status = 0 ;
        goto Exit ;
    }
    status = -2 ;
    if ( *dpsi >= ZERO )
    {
        *b = Com->alpha ;
        *dpsib = *dpsi ;
        goto Exit ;
    }
    else
    {
        if ( psi <= Com->fpert )
        {
            *a = Com->alpha ;
            *dpsia = *dpsi ;
            goto Exit ;
        }
    }
    nshrink = 0 ;
    *b = Com->alpha ;
    while ( TRUE )
    {
        Com->alpha = .5*(*a + *b) ;
        nshrink++ ;
        if ( nshrink > Parm->nexpand )
        {
            status = 8 ;
            goto Exit ;
        }
        cg_evaluate ("fg", "n", Com) ;
        *dpsi = Com->df - Com->wolfe_hi ;
        psi = Com->f - Com->alpha*Com->wolfe_hi ;
        if ( Parm->PrintLevel >= 2 )
        {
            printf ("contract, a: %14.6e alpha: %14.6e phi: "
                    "%14.6e dphi: %14.6e\n", *a, Com->alpha, Com->f, Com->df) ;
        }
        if ( cg_Wolfe (Com->alpha, Com->f, Com->df, Com) )
        {
            status = 0 ;
            goto Exit ;
        }
        if ( *dpsi >= ZERO )
        {
            *b = Com->alpha ;
            *dpsib = *dpsi ;
            status = -1 ;
            goto Exit ;
        }
        if ( psi <= Com->fpert )
        {
            if ( Parm->PrintLevel >= 2 )
            {
                printf ("update a: %14.6e dpsia: %14.6e\n", Com->alpha, *dpsi) ;
            }
            *a = Com->alpha ;
            *dpsia = *dpsi ;
        }
        else *b = Com->alpha ;
    }
Exit:
    if ( Parm->PrintLevel >= 2 )
    {
        printf ("UP a: %14.6e b: %14.6e da: %14.6e db: %14.6e status: %i\n",
                 *a, *b, *dpsia, *dpsib, status) ;
    }
    return (status) ;
}

/* =========================================================================
   ==== cg_fg_evaluate =====================================================
   Evaluate the function and/or gradient.  Also, possibly check if either is nan
   and if so, then reduce the stepsize. Only used at the start of an iteration.
   =========================================================================*/

PRIVATE int cg_evaluate
(
    char    *what, /* fg = evaluate func and grad, g = grad only,f = func only*/
    char     *nan, /* y means check function/derivative values for nan */
    cg_com   *Com
)
{
    INT n ;
    int i ;
    double alpha, *d, *gtemp, *x, *xtemp ;
    cg_parameter *Parm ;
    Parm = Com->Parm ;
    n = Com->n ;
    x = Com->x ;
    d = Com->d ;
    xtemp = Com->xtemp ;
    gtemp = Com->gtemp ;
    alpha = Com->alpha ;
    /* check to see if values are nan */
    if ( !strcmp (nan, "y") )
    {
        if ( !strcmp (what, "f") ) /* compute function */
        {
            cg_step (xtemp, x, d, alpha, n) ;
            /* provisional function value */
            Com->f = Com->cg_value (xtemp, n) ;
            Com->nf++ ;
    
            /* reduce stepsize if function value is nan */
            if ( Com->f != Com->f )
            {
                for (i = 0; i < Parm->nexpand; i++)
                {
                    alpha *= Parm->nan_decay ;
                    cg_step (xtemp, x, d, alpha, n) ;
                    Com->f = Com->cg_value (xtemp, n) ;
                    Com->nf++ ;
                    if ( Com->f == Com->f ) break ;
                }
                if ( i == Parm->nexpand ) return (-2) ;
            }
            Com->alpha = alpha ;
        }
        else                            /* compute gradient */
        {
            cg_step (xtemp, x, d, alpha, n) ;
            Com->cg_grad (gtemp, xtemp, n) ;
            Com->ng++ ;
            Com->df = cg_dot (gtemp, d, n) ;
            /* reduce stepsize if derivative is nan */
            if ( Com->df != Com->df )
            {
                for (i = 0; i < Parm->nexpand; i++)
                {
                    alpha *= Parm->nan_decay ;
                    cg_step (xtemp, x, d, alpha, n) ;
                    Com->cg_grad (gtemp, xtemp, n) ;
                    Com->ng++ ;
                    Com->df = cg_dot (gtemp, d, n) ;
                    if ( Com->df == Com->df ) break ;
                }
                if ( i == Parm->nexpand ) return (-2) ;
                Com->rho = Parm->nan_rho ;
            }
            else Com->rho = Parm->rho ;
            Com->alpha = alpha ;
        }
    }
    else                                /* evaluate without nan checking */
    {
        if ( !strcmp (what, "fg") )     /* compute function and gradient */
        {
            if ( alpha == ZERO )        /* evaluate at x */
            {
                if ( Com->cg_valgrad != NULL )
                {
                    Com->f = Com->cg_valgrad (Com->g, x, n) ;
                }
                else
                {
                    Com->cg_grad (Com->g, x, n) ;
                    Com->f = Com->cg_value (x, n) ;
                }
            }
            else
            {
                cg_step (xtemp, x, d, alpha, n) ;
                if ( Com->cg_valgrad != NULL )
                {
                    Com->f = Com->cg_valgrad (gtemp, xtemp, n) ;
                }
                else
                {
                    Com->cg_grad (gtemp, xtemp, n) ;
                    Com->f = Com->cg_value (xtemp, n) ;
                }
                Com->df = cg_dot (gtemp, d, n) ;
            }
            Com->nf++ ;
            Com->ng++ ;
        }
        else if ( !strcmp (what, "f") ) /* compute function */
        {
            cg_step (xtemp, x, d, alpha, n) ;
            Com->f = Com->cg_value (xtemp, n) ;
            Com->nf++ ;
        }
        else
        {
            cg_step (xtemp, x, d, alpha, n) ;
            Com->cg_grad (gtemp, xtemp, n) ;
            Com->df = cg_dot (gtemp, d, n) ;
            Com->ng++ ;
        }
    }
    return (0) ;
}

/* =========================================================================
   ==== cg_dot =============================================================
   =========================================================================
   Compute dot product of x and y, vectors of length n
   ========================================================================= */
PRIVATE double cg_dot
(
    double *x, /* first vector */
    double *y, /* second vector */
    INT     n /* length of vectors */
)
{
    INT i, n5 ;
    double t ;
    t = 0. ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) t += x [i]*y [i] ;
    for (; i < n; i += 5)
    {
        t += x [i]*y[i] + x [i+1]*y [i+1] + x [i+2]*y [i+2]
                        + x [i+3]*y [i+3] + x [i+4]*y [i+4] ;
    }
    return (t) ;
}

/* =========================================================================
   === cg_copy =============================================================
   =========================================================================
   Copy vector x into vector y
   ========================================================================= */
PRIVATE void cg_copy
(
    double *y, /* output of copy */
    double *x, /* input of copy */
    int     n  /* length of vectors */
)
{
    int j, n10 ;
    n10 = n % 10 ;
    for (j = 0; j < n10; j++) y [j] = x [j] ;
    for (; j < n; j += 10)
    {
        y [j] = x [j] ;
        y [j+1] = x [j+1] ;
        y [j+2] = x [j+2] ;
        y [j+3] = x [j+3] ;
        y [j+4] = x [j+4] ;
        y [j+5] = x [j+5] ;
        y [j+6] = x [j+6] ;
        y [j+7] = x [j+7] ;
        y [j+8] = x [j+8] ;
        y [j+9] = x [j+9] ;
    }
}

/* =========================================================================
   ==== cg_step ============================================================
   =========================================================================
   Compute xtemp = x + alpha d
   ========================================================================= */
PRIVATE void cg_step
(
    double *xtemp, /*output vector */
    double     *x, /* initial vector */
    double     *d, /* search direction */
    double  alpha, /* stepsize */
    INT         n  /* length of the vectors */
)
{
    INT n5, i ;
    n5 = n % 5 ;
    for (i = 0; i < n5; i++) xtemp [i] = x[i] + alpha*d[i] ;
    for (; i < n; i += 5)
    { 
        xtemp [i]   = x [i]   + alpha*d [i] ;
        xtemp [i+1] = x [i+1] + alpha*d [i+1] ;
        xtemp [i+2] = x [i+2] + alpha*d [i+2] ;
        xtemp [i+3] = x [i+3] + alpha*d [i+3] ;
        xtemp [i+4] = x [i+4] + alpha*d [i+4] ;
    }
}

/* =========================================================================
   === cg_default ==========================================================
   =========================================================================
   Set default conjugate gradient parameter values. If the parameter argument
   of cg_descent is NULL, this routine is called by cg_descent automatically.
   If the user wishes to set parameter values, then the cg_parameter structure
   should be allocated in the main program. The user could call cg_default
   to initialize the structure, and then individual elements in the structure
   could be changed, before passing the structure to cg_descent.
   =========================================================================*/
void cg_default
(
    cg_parameter   *Parm
)
{
    /* T => print final function value
       F => no printout of final function value */
    Parm->PrintFinal = TRUE ;

   /* Level 0 = no printing, ... , Level 3 = maximum printing */
    Parm->PrintLevel = 0 ;

    /* T => print parameters values
       F => do not display parmeter values */
    Parm->PrintParms = FALSE ;

    /* T => use approximate Wolfe line search
       F => use ordinary Wolfe line search, switch to approximate Wolfe when
                |f_k+1-f_k| < AWolfeFac*C_k, C_k = average size of cost */
    Parm->AWolfe = FALSE ;
    Parm->AWolfeFac = 1.e-3 ;

    /* factor in [0, 1] used to compute average cost magnitude C_k as follows:
       Q_k = 1 + (Qdecay)Q_k-1, Q_0 = 0,  C_k = C_k-1 + (|f_k| - C_k-1)/Q_k */
    Parm->Qdecay = .7 ;

    /* Stop Rules:
       T => ||grad||_infty <= max(grad_tol, initial |grad|_infty*StopFact)
       F => ||grad||_infty <= grad_tol*(1 + |f_k|) */
    Parm->StopRule = TRUE ;
    Parm->StopFac = 0.e-12 ;

    /* T => estimated error in function value is eps*Ck,
       F => estimated error in function value is eps */
    Parm->PertRule = TRUE ;
    Parm->eps = 1.e-6 ;

    /* T => attempt quadratic interpolation in line search when
                |f_k+1 - f_k|/f_k > QuadCutOff
       F => no quadratic interpolation step */
    Parm->QuadStep = TRUE ;
    Parm->QuadCutOff = 1.e-12 ;

    /* T => check that f_k+1 - f_k <= debugtol*C_k
       F => no checking of function values */
    Parm->debug = FALSE ;
    Parm->debugtol = 1.e-10 ;

    /* if step is nonzero, it is the initial step of the initial line search */
    Parm->step = ZERO ;

    /* abort cg after maxit_fac*n iterations */
    Parm->maxit_fac = INF ;

    /* maximum number of times the bracketing interval grows or shrinks
       in the line search is nexpand */
    Parm->nexpand = (int) 50 ;

    /* maximum number of secant iterations in line search is nsecant */
    Parm->nsecant = (int) 50 ;

    /* conjugate gradient method restarts after (n*restart_fac) iterations */
    Parm->restart_fac = 6.0 ;

    /* stop when -alpha*dphi0 (estimated change in function value) <= feps*|f|*/
    Parm->feps = ZERO ;

    /* after encountering nan, growth factor when searching for
       a bracketing interval */
    Parm->nan_rho = 1.3 ;

    /* after encountering nan, decay factor for stepsize */
    Parm->nan_decay = 0.1 ;

    /* Wolfe line search parameter, range [0, .5]
       phi (a) - phi (0) <= delta phi'(0) */
    Parm->delta = .1 ;

    /* Wolfe line search parameter, range [delta, 1]
       phi' (a) >= sigma phi' (0) */
    Parm->sigma = .9 ;

    /* decay factor for bracket interval width in line search, range (0, 1) */
    Parm->gamma = .66 ;

    /* growth factor in search for initial bracket interval */
    Parm->rho = 5. ;

    /* starting guess for line search =
         psi0 ||x_0||_infty over ||g_0||_infty if x_0 != 0
         psi0 |f(x_0)|/||g_0||_2               otherwise */
    Parm->psi0 = .01 ;      /* factor used in starting guess for iteration 1 */

    /* for a QuadStep, function evalutated at psi1*previous step */
    Parm->psi1 = .1 ;

    /* when starting a new cg iteration, our initial guess for the line
       search stepsize is psi2*previous step */
    Parm->psi2 = 2. ;

    /* choose theta adaptively if AdaptiveBeta = T */
    Parm->AdaptiveBeta = FALSE ;

    /* lower bound for beta is BetaLower*d_k'g_k/ ||d_k||^2 */
    Parm->BetaLower = 0.4 ;

    /* value of the parameter theta in the cg_descent update formula:
       W. W. Hager and H. Zhang, A survey of nonlinear conjugate gradient
       methods, Pacific Journal of Optimization, 2 (2006), pp. 35-58. */
    Parm->theta = 1.0 ;

    /* parameter used in cost error estimate for quadratic restart criterion */
    Parm->qeps = 1.e-12 ;

    /* number of iterations the function is nearly quadratic before a restart */
    Parm->qrestart = 3 ;

    /* treat cost as quadratic if
       |1 - (cost change)/(quadratic cost change)| <= qrule */
    Parm->qrule = 1.e-8 ;
}

/* =========================================================================
   ==== cg_printParms ======================================================
   =========================================================================
   Print the contents of the cg_parameter structure
   ========================================================================= */
PRIVATE void cg_printParms
(
    cg_parameter  *Parm
)
{
    printf ("PARAMETERS:\n") ;
    printf ("\n") ;
    printf ("Wolfe line search parameter ..................... delta: %e\n",
             Parm->delta) ;
    printf ("Wolfe line search parameter ..................... sigma: %e\n",
             Parm->sigma) ;
    printf ("decay factor for bracketing interval ............ gamma: %e\n",
             Parm->gamma) ;
    printf ("growth factor for bracket interval ................ rho: %e\n",
             Parm->rho) ;
    printf ("growth factor for bracket interval after nan .. nan_rho: %e\n",
             Parm->nan_rho) ;
    printf ("decay factor for stepsize after nan ......... nan_decay: %e\n",
             Parm->nan_decay) ;
    printf ("parameter in lower bound for beta ........... BetaLower: %e\n",
             Parm->BetaLower) ;
    printf ("parameter describing cg_descent family .......... theta: %e\n",
             Parm->theta) ;
    printf ("perturbation parameter for function value ......... eps: %e\n",
             Parm->eps) ;
    printf ("factor for computing average cost .............. Qdecay: %e\n",
             Parm->Qdecay) ;
    printf ("relative change in cost to stop quadstep ... QuadCufOff: %e\n",
             Parm->QuadCutOff) ;
    printf ("factor multiplying gradient in stop condition . StopFac: %e\n",
             Parm->StopFac) ;
    printf ("cost change factor, approx Wolfe transition . AWolfeFac: %e\n",
             Parm->AWolfeFac) ;
    printf ("restart cg every restart_fac*n iterations . restart_fac: %e\n",
             Parm->restart_fac) ;
    printf ("cost error in quadratic restart is qeps*cost ..... qeps: %e\n",
             Parm->qeps) ;
    printf ("number of quadratic iterations before restart  qrestart: %i\n",
             Parm->qrestart) ;
    printf ("parameter used to decide if cost is quadratic ... qrule: %e\n",
             Parm->qrule) ;
    printf ("stop when cost change <= feps*|f| ................. eps: %e\n",
             Parm->eps) ;
    printf ("starting guess parameter in first iteration ...... psi0: %e\n",
             Parm->psi0) ;
    printf ("starting step in first iteration if nonzero ...... step: %e\n",
             Parm->step) ;
    printf ("factor multiply starting guess in quad step ...... psi1: %e\n",
             Parm->psi1) ;
    printf ("initial guess factor for general iteration ....... psi2: %e\n",
             Parm->psi2) ;
    printf ("max iterations is n*maxit_fac ............... maxit_fac: %e\n",
             Parm->maxit_fac) ;
    printf ("max expansions in line search ................. nexpand: %i\n",
             Parm->nexpand) ;
    printf ("max secant iterations in line search .......... nsecant: %i\n",
             Parm->nsecant) ;
    printf ("print level (0 = none, 2 = maximum) ........ PrintLevel: %i\n",
             Parm->PrintLevel) ;
    printf ("Logical parameters:\n") ;
    if ( Parm->PertRule )
        printf ("    Error estimate for function value is eps\n") ;
    else
        printf ("    Error estimate for function value is eps*Ck\n") ;
    if ( Parm->QuadStep )
        printf ("    Use quadratic interpolation step\n") ;
    else
        printf ("    No quadratic interpolation step\n") ;
    if ( Parm->AdaptiveBeta )
        printf ("    Adaptively adjust direction update parameter beta\n") ;
    else
        printf ("    Use fixed parameter theta in direction update\n") ;
    if ( Parm->PrintFinal )
        printf ("    Print final cost and statistics\n") ;
    else
        printf ("    Do not print final cost and statistics\n") ;
    if ( Parm->PrintParms )
        printf ("    Print the parameter structure\n") ;
    else
        printf ("    Do not print parameter structure\n") ;
    if ( Parm->AWolfe)
        printf ("    Approximate Wolfe line search\n") ;
    else
        printf ("    Wolfe line search") ;
        if ( Parm->AWolfeFac > 0. )
            printf (" ... switching to approximate Wolfe\n") ;
        else
            printf ("\n") ;
    if ( Parm->StopRule )
        printf ("    Stopping condition uses initial grad tolerance\n") ;
    else
        printf ("    Stopping condition weighted by absolute cost\n") ;
    if ( Parm->debug)
        printf ("    Check for decay of cost, debugger is on\n") ;
    else
        printf ("    Do not check for decay of cost, debugger is off\n") ;
}

/*
Version 1.2 Change:
  1. The variable dpsi needs to be included in the argument list for
     subroutine cg_updateW (update of a Wolfe line search)

Version 2.0 Changes:
    The user interface was redesigned. The parameters no longer need to
    be read from a file. For compatibility with earlier versions of the
    code, we include the routine cg_readParms to read parameters.
    In the simplest case, the user can use NULL for the
    parameter argument of cg_descent, and the code sets the default
    parameter values. If the user wishes to modify the parameters, call
    cg_default in the main program to initialize a cg_parameter
    structure. Individual elements of the structure could be modified.
    The header file cg_user.h contains the structures and prototypes
    that the user may need to reference or modify, while cg_descent.h
    contains header elements that only cg_descent will access.  Note
    that the arguments of cg_descent have changed.

Version 3.0 Changes:
    Major overhaul

Version 4.0 Changes:
  1. Set theta = 1.0 by default in the cg_descent rule for beta_k
  2. Increase the default value of restart_fac to 6 (a value larger than 1 is
     more efficient when the problem dimension is small)
  3. Restart the CG iteration if the objective function is nearly quadratic
     for several iterations (qrestart). This type of restart is described
     in the paper "A nonlinear conjugate gradient algorithm with an
     optimal property and an improved Wolfe line search by Yu-Hong Dai and
     Cai-Xia Kou.
  4. New lower bound for beta: BetaLower*d_k'g_k/ ||d_k||^2
  5. Evaluation of the objective function and gradient is now handled by
     the routine cg_evaluate.
*/
