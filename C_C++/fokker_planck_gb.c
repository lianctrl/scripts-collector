/* do_bayesian_analysis.c */

#include "mkl_types.h"
#include "mkl_cblas.h"
#include "mkl_lapack.h"
#include "general.h"
#include "do_bayesian_analysis.h"
#include "sa.h"
#include "cg.h"
#include "linear_algebra.h"

/* Global static variables to be used in SA func */
static  int  SA_func_ndof;
static  real *SAfunc_parm1;
static  real SAfunc_parm2;
static  real **SAfunc_parm3;


/********************************************************************/ 
/*                                                                  */
/*                  ROUTINE: main                                   */
/*                                                                  */
/********************************************************************/

int main(int argc, char **argv)
{
/*==================================================================*/
/*            Local variable declarations                           */ 

    real gamma, start, end, delta, au2eV, ix_delta, xj, pij;
    real *x, *bin, *P, Ptot, dt, twopi, pi, minbin, maxbin;
    int long ranseed = -4321234;
    real **Rij, *Ri, *P_t, Pi_t, xc, *Ptheo, **Rij_theo, **Rij_theo2, Rij_sum, *Dtheo, *Dfit;
    real **M_t0, **M_t1, **M_t2, **M_t3, **M_t4, **dM_t0, **M_t0_inv, **M_t0_theo, **M_t0_theo_sym;
    real *Ddyn, dist, sign;
    int  nbin, npoints, lag, nskip, iter, ndt, nlagtime;
    int  info, *ipiv, lda, lwork;
    real *work, temp, temp2, lagtime;
    real **V, **VT, **V1, **V2, **Diag, **Diag2, *Lambda;
    int  nv;
    int  n, m, i, j, k, ix, jx;
    FILE *fin1, *fin2, *fout, *fout2;

/*==================================================================*/
/*  I)        Read input data                                       */

    Usage(argc, argv);


    if ((fin1 = fopen(argv[1],"r"))==NULL) {
      printf("Cannot open input data file!\n");
      exit(1);
    }

    npoints = atol(argv[2]);
    nbin    = atol(argv[3]);
    nskip   = atol(argv[4]);
    lag     = atol(argv[5]);

    if ((fout = fopen(argv[6],"w"))==NULL) {
      printf("Cannot create output file!\n");
      exit(1);
    }

/*==================================================================*/
/*  II)       Read input data                                       */

    pi = 3.141592653589793;
    twopi = 2.0 * 3.141592653589793;
 
    dt  = 0.004;
    lagtime = dt * lag;

    x   = AllocVecR(npoints);
    Rij = AllocMatR(nbin,nbin);
    Ri  = AllocVecR(nbin);
    P   = AllocVecR(nbin);
    P_t = AllocVecR(nbin);
    Ptheo = AllocVecR(nbin);
    Dtheo = AllocVecR(nbin);
    Dfit  = AllocVecR(nbin);
    Ddyn  = AllocVecR(nbin);
    bin = AllocVecR(nbin+1);
    delta = twopi / ((real) nbin);

    /* MAYBE MODIFY THIS */
    ndt      = 5;
    M_t0     = AllocMatR(nbin,nbin);
    M_t1     = AllocMatR(nbin,nbin);
    M_t2     = AllocMatR(nbin,nbin);
    M_t3     = AllocMatR(nbin,nbin);
    M_t4     = AllocMatR(nbin,nbin);
    dM_t0    = AllocMatR(nbin,nbin);
    M_t0_inv = AllocMatR(nbin,nbin);

    Rij_theo  = AllocMatR(nbin,nbin);
    Rij_theo2 = AllocMatR(nbin,nbin);
    M_t0_theo = AllocMatR(nbin,nbin);
    M_t0_theo_sym = AllocMatR(nbin,nbin);

    V1       = AllocMatR(nbin,nbin);
    V2       = AllocMatR(nbin,nbin);
    V        = AllocMatR(nbin,nbin);
    VT       = AllocMatR(nbin,nbin);
    Diag     = AllocMatR(nbin,nbin);
    Diag2    = AllocMatR(nbin,nbin);
    Lambda   = AllocVecR(nbin);

/*
    mu     = AllocMatR(npoints,3);
    freq   = AllocVecR(npoints);
    osc    = AllocVecR(npoints);
    alphaR = AllocMatR(3,3);
    alphaI = AllocMatR(3,3);
*/

    minbin = 6.6;
    maxbin = 9.1;

    bin[1] = minbin;
    for (i = 1; i <= nbin; i++) {
      Ri[i] = 0.0;
      P[i]  = 0.0;
      P_t[i]  = 0.0;
      Ptheo[i]  = 0.0;
      Dtheo[i] = 0.0;
      Dfit[i] = 0.0;
      Ddyn[i] = 0.0;
      bin[i + 1] = i *  (maxbin - minbin) / nbin + minbin;
      for (j = 1; j <= nbin; j++) {
        Rij[i][j] = 0.0;
        M_t0[i][j] = 0.0;
        M_t1[i][j] = 0.0;
        M_t2[i][j] = 0.0;
        M_t3[i][j] = 0.0;
        M_t4[i][j] = 0.0;
        dM_t0[i][j] = 0.0;
        M_t0_inv[i][j] = 0.0;
        Rij_theo[i][j] = 0.0;
        Rij_theo2[i][j] = 0.0;
        M_t0_theo[i][j] = 0.0;
        M_t0_theo_sym[i][j] = 0.0;

        V[i][j]     = 0.0;
        VT[i][j]    = 0.0;
        V1[i][j]    = 0.0;
        V2[i][j]    = 0.0;
        Diag[i][j]  = 0.0;
        Diag2[i][j] = 0.0;
        
      }
    }

#ifdef out
    V1[1][1] =  1.96;   V1[1][2] =  0.00;   V1[1][3] =  0.00;   V1[1][4] = 0.00;   V1[1][5] =  0.00;
    V1[2][1] = -6.49;   V1[2][2] =  3.80;   V1[2][3] =  0.00;   V1[2][4] = 0.00;   V1[2][5] =  0.00;
    V1[3][1] = -0.47;   V1[3][2] = -6.39;   V1[3][3] =  4.17;   V1[3][4] = 0.00;   V1[3][5] =  0.00;
    V1[4][1] = -7.20;   V1[4][2] =  1.50;   V1[4][3] = -1.51;   V1[4][4] = 5.70;   V1[4][5] =  0.00;
    V1[5][1] = -0.65;   V1[5][2] = -6.34;   V1[5][3] =  2.67;   V1[5][4] = 1.80;   V1[5][5] = -7.10;

/*
    V1[1][1] =  1.96;   V1[2][1] =  0.00;   V1[3][1] =  0.00;   V1[4][1] = 0.00;   V1[5][1] =  0.00;
    V1[1][2] = -6.49;   V1[2][2] =  3.80;   V1[3][2] =  0.00;   V1[4][2] = 0.00;   V1[5][2] =  0.00;
    V1[1][3] = -0.47;   V1[2][3] = -6.39;   V1[3][3] =  4.17;   V1[4][3] = 0.00;   V1[5][3] =  0.00;
    V1[1][4] = -7.20;   V1[2][4] =  1.50;   V1[3][4] = -1.51;   V1[4][4] = 5.70;   V1[5][4] =  0.00;
    V1[1][5] = -0.65;   V1[2][5] = -6.34;   V1[3][5] =  2.67;   V1[4][5] = 1.80;   V1[5][5] = -7.10;
*/
    for (i = 1; i <= 5; i++)
      for (j = i + 1; j <= 5; j++)
        V1[i][j] = V1[j][i];        

    k = 1;
    for (i = 1; i <= 5; i++) {
      Lambda[i]  = 0.0;
      for (j = 1; j <= 5; j++) {
        V[i][j] = 0.0;
        VT[i][j] = 0.0;
        Diag[i][j] = 0.0;
        V2[i][j] = 0.0;
      }
    }

    nv = 5; 
    fprintf(stdout," TUTTO OK!!!\n\n"); fflush(stdout);
 
    DiagonalMatrix(nv, V1, VT, Lambda);

    TransposeMatrix(nv, VT, V);
    ProductMatrix(nv, VT, V1, V2);
    ProductMatrix(nv, V2, V, Diag);


    fprintf(stdout," TUTTO OK!!!\n\n"); fflush(stdout);

    fprintf(fout,"\n\n\n");
    fprintf(fout," Diagonal\n");
    for (i = 1; i <= nv; i++) {
      for (j = 1; j <= nv; j++) {
         fprintf(fout,"%12.8lf", V1[i][j]);
      }
      fprintf(fout,"\n");
    }
    fprintf(fout," Diagonal Eigenvector\n");
    for (i = 1; i <= nv; i++) {
      for (j = 1; j <= nv; j++) {
         fprintf(fout,"%12.8lf", V[i][j]);
      }
      fprintf(fout,"\n");
    }
    fprintf(fout," Diagonal Eigenvalue\n");
    for (i = 1; i <= nv; i++) {
      fprintf(fout,"%12.8lf\n", Lambda[i]);
    }

    fprintf(fout," Diag\n");
    for (i = 1; i <= nv; i++) {
      for (j = 1; j <= nv; j++) {
         fprintf(fout,"%12.8lf", Diag[i][j]);
      }
      fprintf(fout,"\n");
    }
    
    exit(0);
#endif


    for (i = 1; i <= nbin; i++) {
      fprintf(stdout," Bin %5d   %8.5lf  %8.5lf\n", i, bin[i], bin[i+1]);
      fflush(stdout);

    }

    /* Read in data */
    for (n = 1; n <= npoints; n++) {
      if (fscanf(fin1,"%lf", &x[n]) == EOF)  break;
    }

    /* Compute probability distribution */
    for (n = 1; n <= npoints; n++) {
      for (i = 1; i <= nbin; i++) {
        if ((x[n] >= bin[i]) && (x[n] < bin[i+1])) {
          P[i] = P[i] + 1.0;
        }
      }
    }
 
    Ptot = 0.0;
    for (i = 1; i <= nbin; i++)
      Ptot = Ptot + P[i];

    for (i = 1; i <= nbin; i++)
      P[i] = P[i] / Ptot;

/*==================================================================*/
/*  II)       Compute diffusion from dynamics                       */

    
    for (n = 1; n <= npoints; n = n + nskip) {
      if ((n + lag) > npoints) continue;

      temp = x[n+ndt] - x[n];
      // dist = fabs(temp);
      // sign = temp / dist;

      // if (dist > pi) temp = temp - sign*twopi;
            
      for (i = 1; i <= nbin; i++) {
        if ((x[n] >= bin[i]) && (x[n] < bin[i+1])) {
          P_t[i]  = P_t[i] + 1.0;
          Ddyn[i] += temp * temp / (2.0 * ndt * dt);
        }
      }

    }

    for (i = 1; i <= nbin; i++) {
      Ddyn[i] = Ddyn[i] / P_t[i];
      P_t[i]  = 0;
    }

/*==================================================================*/
/*  II)       Compute exp(tR) at times t=t', t=t'-dt and t=t'+dt    */


    for (iter = -2; iter <= 2; iter++) {

      nlagtime = lag + iter * ndt;

      /* Initialize matrix */
      for (i = 1; i <= nbin; i++) {
        Ri[i] = 0.0;
        for (j = 1; j <= nbin; j++) {
          Rij[i][j] = 0.0;
        }
      }

    for (n = 1; n <= npoints; n = n + nskip) {
      ix = 0;
      jx = 0;

      for (i = 1; i <= nbin; i++) {
        if ((x[n] >= bin[i]) && (x[n] < bin[i+1])) {
          ix = i;
          ix_delta = x[n] - (bin[i+1] + bin[i]) / 2.0;
        }
      }

      if ((n + nlagtime) > npoints) continue;

      xj = x[n+nlagtime] - ix_delta;

      // Apply PBC
/*
      if (xj > pi) xj = xj - twopi;
      if (xj < -1.0 * pi) xj = xj + twopi;
*/

      // xj = x[n+nlagtime];
      for (j = 1; j <= nbin; j++) {
        if ((xj >= bin[j]) && (xj < bin[j+1])) {
          jx = j;
        }
      }


      // fprintf(stdout," xi %10.5lf i %3d    xj %10.5lf %10.5lf  j %3d\n", x[n], ix, x[n+nlagtime], xj, jx);

      /* Note: indexes inverted so that (j,t=0) (i,t=t')*/
      if ((ix == 0) || (jx == 0)) 
        fprintf(stdout," Problem!!!! %d\n", n);
      fflush(stdout);

      if ((ix != 0) && (jx != 0))
        Rij[jx][ix] = Rij[jx][ix] + 1.0;
    }

    /* Symmetrize matrix to enforce detailed balance */
/*
    for (i = 1; i <= nbin - 1; i++) {
      for (j = i + 1; j <= nbin; j++) {
        pij = (Rij[i][j] + Rij[j][i]) / 2.0;
        Rij[i][j] = pij;
        Rij[j][i] = pij;
      }
    }
*/

      for (i = 1; i <= nbin; i++) {
        for (j = 1; j <= nbin; j++) {
          Ri[i] = Ri[i] + Rij[j][i];
        }
      }

        for (i = 1; i <= nbin; i++) {
          for (j = 1; j <= nbin; j++) {
            if (iter == -1) {
              M_t1[i][j] = ((real) Rij[i][j]) / ((real) Ri[j]);
            } else if (iter == 0) {
              M_t0[i][j] = ((real) Rij[i][j]) / ((real) Ri[j]);
            } else if (iter == 1) {
              M_t2[i][j] = ((real) Rij[i][j]) / ((real) Ri[j]);
            } else if (iter == -2) {
              M_t3[i][j] = ((real) Rij[i][j]) / ((real) Ri[j]);
            } else if (iter == 2) {
              M_t4[i][j] = ((real) Rij[i][j]) / ((real) Ri[j]);
            }
          }
        }


/*
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        if (Rij[i][j] > 0) {
          pij = ((real) Rij[i][j]) / ((real) Ri[j]);
          fprintf(fout," %3d %3d Nij %8.1lf Nj %8.5lf  Mij %8.5lf  MijPj %8.5lf  Pj %8.5lf\n", i, j, Rij[i][j], Ri[j], pij, pij * P[j], P[j]);
          fflush(fout);
        }
      }
    }
*/

    } /*endfor iter*/

/*==================================================================*/
/*  II)       Print out exp(tR) at times t=t', t=t'-dt and t=t'+dt  */

    /* Theoretical values */
    for (i = 1; i <= nbin; i++) {
      xc = (bin[i+1] + bin[i]) / 2.0;
      Ptheo[i] = delta * exp(-1.0*(-cos(2*xc)+2.0));
    }

    Ptot = 0.0;
    for (i = 1; i <= nbin; i++)
      Ptot = Ptot + Ptheo[i];

    for (i = 1; i <= nbin; i++)
      Ptheo[i] = Ptheo[i] / Ptot;
    /* End theoretical values */

    Pi_t = 0.0;
    for (i = 1; i <= nbin; i++) {
      temp = 0.0;
      for (j = 1; j <= nbin; j++) {
        temp = temp + M_t0[i][j] * P[j];
      }
      for (j = 1; j <= nbin; j++) {
        if (j == i) continue;
        temp = temp - M_t0[j][i] * P[i];
      }
      Pi_t += temp;
      P_t[i] = temp;
    }

    for (i = 1; i <= nbin; i++) {
      fprintf(fout," %3d  Xbin  %12.5lf -  %12.5lf  Pi_t %8.5lf  Pi %12.8lf  Ptheo %8.5lf\n", i, bin[i], bin[i+1], P_t[i]/Pi_t, P[i], Ptheo[i]);
      fflush(fout);
    }

    fprintf(fout,"\n\n\n");
    fprintf(fout," M_t0\n");
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        fprintf(fout,"%12.8lf", M_t0[i][j]);
      }
      fprintf(fout,"\n");
    }

    fprintf(fout,"\n\n\n");
    fprintf(fout," M_t1\n");
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        fprintf(fout,"%12.8lf", M_t1[i][j]);
      }
      fprintf(fout,"\n");
    }

    fprintf(fout,"\n\n\n");
    fprintf(fout," M_t2\n");
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        fprintf(fout,"%12.8lf", M_t2[i][j]);
      }
      fprintf(fout,"\n");
    }


/*==================================================================*/
/*  II)       Compute the inverse of exp(tR)                        */

    info  = 0;
    lda   = nbin;
    work  = AllocVecR(64 * nbin);
    lwork = 64 * nbin;
    ipiv  = AllocVecI(nbin);

    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        M_t0_inv[i][j] = M_t0[i][j];
      }
    }

    dgetrf(&nbin, &nbin, &(M_t0_inv[1][1]), &nbin, &ipiv[1], &info);
    if (info != 0)  {
      fprintf(stdout,"     Error in routine dgetrf (info: %3d)\n", info);
      exit(1);
    }
    dgetri(&nbin, &(M_t0_inv[1][1]), &nbin, &ipiv[1], &work[1], &lwork, &info);
    if (info != 0)  {
      fprintf(stdout,"     Error in routine dgetri (info: %3d)\n", info);
      exit(1);
    }

    fprintf(fout,"\n\n\n");
    fprintf(fout," M_t0_inv\n");
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        fprintf(fout,"%12.8lf", M_t0_inv[i][j]);
      }
      fprintf(fout,"\n");
    }

/*==================================================================*/
/*  II)       Compute d exp(tR) / dt and R                          */


    fprintf(fout,"\n\n\n");
    fprintf(fout," dM 2nd order dt=%.5lf\n", dt * ndt);
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        dM_t0[i][j] = ((1.0/2.0)*M_t2[i][j] + (-1.0/2.0)*M_t1[i][j]) / (dt * ndt);
        fprintf(fout,"%12.8lf", dM_t0[i][j]);
      }   
      fprintf(fout,"\n");
    }

    fprintf(fout,"\n\n\n");
    fprintf(fout," dM 4th order dt=%.5f\n", dt * ndt);
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        dM_t0[i][j] = ((-1.0/12.0)*M_t4[i][j] + (2.0/3.0)*M_t2[i][j] + (-2.0/3.0)*M_t1[i][j] + (1.0/12.0)*M_t3[i][j]) / (dt * ndt);
        fprintf(fout,"%12.8lf", dM_t0[i][j]);
      }
      fprintf(fout,"\n");
    }

    fprintf(fout,"\n\n\n");
    fprintf(fout," R\n");
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        Rij[i][j] = 0.0;
        for (k = 1; k <= nbin; k++)
          Rij[i][j] = Rij[i][j] + dM_t0[i][k] * M_t0_inv[k][j];

        fprintf(fout,"%12.8lf", Rij[i][j]);
      }
      fprintf(fout,"\n");
    }

    fprintf(fout,"\n\n\n");
    fprintf(fout," Rii Sum Rji\n");
    for (i = 1; i <= nbin; i++) {
      Rij_sum = 0.0;
      for (j = 1; j <= nbin; j++) {
        if (j != i)
          Rij_sum = Rij_sum - Rij[j][i];
      }

      fprintf(fout,"%3d  %12.8lf   %12.8lf\n", i, Rij[i][i], Rij_sum);
    }

    fprintf(fout,"\n\n\n");
    fprintf(fout," R_i+1,i\n");
    for (i = 1; i <= nbin - 1; i++) {


      xc = (bin[i+1] + bin[i]) / 2.0;

//      xc = bin[i+1];
      Rij_theo[i+1][i] = (0.1*(2.0 + sin(xc))) * sqrt(Ptheo[i+1]/Ptheo[i]) / (delta*delta);

      fprintf(fout,"%3d  %3d   %12.8lf   %12.8lf   %12.8lf\n", i+1, i, Rij[i+1][i], Rij_theo[i+1][i], 
      Rij[i][i+1]);
    }
    i = 1;
    Rij_theo[1][nbin] = (0.1*(2.0 + sin(pi))) * sqrt(Ptheo[i]/Ptheo[nbin]) / (delta*delta);
    fprintf(fout,"%3d  %3d   %12.8lf   %12.8lf  %12.8lf\n", i, nbin, Rij[1][nbin], Rij_theo[1][nbin],
    Rij[nbin][1]);


    /* Further check */
    for (i = 1; i <= nbin - 1; i++) {
      Rij_theo[i][i+1] = Rij_theo[i+1][i] * Ptheo[i] / Ptheo[i+1];
    }
    Rij_theo[nbin][1] = Rij_theo[1][nbin] * Ptheo[nbin] / Ptheo[1];

    for (i = 1; i <= nbin; i++) {
      temp = 0.0;
      for (j = 1; j <= nbin; j++) {
        if (j != i)
          temp = temp - Rij_theo[j][i];
      }
      Rij_theo[i][i] = temp;
    }

    fprintf(fout,"\n\n\n");
    fprintf(fout," R theo\n");
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        fprintf(fout,"%12.8lf", Rij_theo[i][j]);
      }
      fprintf(fout,"\n");
    }

    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        V1[i][j] = (1.0/sqrt(Ptheo[i])) * Rij_theo[i][j] * sqrt(Ptheo[j]);
      }
    }

    fprintf(fout,"\n\n\n");
    fprintf(fout," R theo sym\n");
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        fprintf(fout,"%12.8lf", V1[i][j]);
      }
      fprintf(fout,"\n");
    }

    DiagonalMatrix(nbin, V1, VT, Lambda);
    TransposeMatrix(nbin, VT, V);

    for (i = 1; i <= nbin; i++) {
        Diag[i][i] = exp(dt * lag * Lambda[i]);
    }

    
    ProductMatrix(nbin, V, Diag, V2);
    ProductMatrix(nbin, V2, VT, M_t0_theo_sym);

    for (i = 1; i <= nbin; i++) {
      Diag[i][i]  = sqrt(Ptheo[i]);
      Diag2[i][i] = 1.0 / sqrt(Ptheo[i]);
    }

    ProductMatrix(nbin, Diag, M_t0_theo_sym, V2);
    ProductMatrix(nbin, V2, Diag2, M_t0_theo);

    fprintf(fout,"\n\n\n");
    fprintf(fout," M_t0 theo from eigenvalue method\n");
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        fprintf(fout,"%12.8lf", M_t0_theo[i][j]);
      }
      fprintf(fout,"\n");
    }

    fprintf(fout,"\n\n\n");
    fprintf(fout," M_t0 theo from taylor expansion\n");
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        V1[i][j] = dt * lag * Rij_theo[i][j];
      }
    }    
    ExponentialMatrix(nbin, V1, M_t0_theo, 60);
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        fprintf(fout,"%12.8lf", M_t0_theo[i][j]);
      }
      fprintf(fout,"\n");
    }

    fprintf(fout,"\n\n\n");
    fprintf(fout," R theo2\n");

    BuildRtheo(nbin, Rij_theo2, bin);
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        fprintf(fout,"%12.8lf", Rij_theo2[i][j]);
      }
      fprintf(fout,"\n");
    }

    fprintf(fout,"\n\n\n");
    fprintf(fout," M_t0 theo2 from taylor expansion\n");
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        V1[i][j] = dt * lag * Rij_theo2[i][j];
      }
    }
    ExponentialMatrix(nbin, V1, M_t0_theo, 60);
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        fprintf(fout,"%12.8lf", M_t0_theo[i][j]);
      }
      fprintf(fout,"\n");
    }


/*
    SimulatedAnnealing (nbin, lagtime, Ptheo, Rij_theo, M_t0);
*/

//    ConjugateGradients (nbin, lagtime, Ptheo, Rij_theo, M_t0);


    ConjugateGradients (nbin, lagtime, P, Rij_theo, M_t0);


    fprintf(fout,"\n\n\n");
    fprintf(fout," R fit sim. anneal.\n");
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        fprintf(fout,"%12.8lf", Rij_theo[i][j]);
      }
      fprintf(fout,"\n");
    }


    fprintf(fout,"\n\n\n");
    fprintf(fout,"            Dtheo      Dfit     Ddyn\n");
    for (i = 1; i <= nbin-1; i++) {

      xc = (bin[i+1] + bin[i]) / 2.0;
//      xc = bin[i+1];

      Dtheo[i] = Rij[i+1][i] * (delta*delta) * sqrt(P[i]/P[i+1]);
      temp2    = Rij[i][i+1] * (delta*delta) * sqrt(P[i+1]/P[i]);

      Dfit[i]  = Rij_theo[i+1][i] * (delta*delta) * sqrt(P[i]/P[i+1]);
      temp = Rij_theo[i][i+1] * (delta*delta) * sqrt(P[i+1]/P[i]);

      fprintf(fout,"%3d  %3d   %12.8lf |  %10.5lf  %10.5lf |  %10.5lf   %10.5lf |  %10.5lf\n", i+1, i, xc, Dtheo[i], temp2, Dfit[i], temp, Ddyn[i]);
    }

/*
    i = 1;
    xc = (bin[nbin+1] + bin[nbin]) / 2.0;

    Dtheo[nbin] = (0.1*(2.0 + sin(xc)));
    Dfit[nbin]  = Rij_theo[1][nbin] * (delta*delta) *  sqrt(Ptheo[nbin]/Ptheo[i]);
    fprintf(fout,"%3d  %3d   %12.8lf   %10.5lf   %10.5lf   %10.5lf\n", i, nbin, xc, Dtheo[nbin], Dfit[nbin], Ddyn[i]);
*/


    fclose(fin1);
    fclose(fout);

    return 0;
}

/********************************************************************/
/*                                                                  */
/*                  ROUTINE: Usage                                  */
/*                                                                  */
/********************************************************************/

int Propagate(real x0, real *x1, real dt, long *ranseed)
{
/*==================================================================*/
/*            Local variable declarations                           */
 
    real D, D1, F1;
    int n;

/*==================================================================*/
/*  I)        Compute next step                                     */

    D  = 0.1 * (2.0 + sin(x0));
    D1 = 0.1 * cos(x0);
    F1 = 2.0 * sin(2.0*x0);

    *x1 = x0  + (D1 - D*F1)*dt + gasdev(ranseed) * sqrt(2*D*dt);

    return 0;
}



/********************************************************************/ 
/*                                                                  */
/*                  ROUTINE: Usage                                  */
/*                                                                  */
/********************************************************************/

int Usage(int argc, char **argv)
{
/*==================================================================*/
/*            Local variable declarations                           */ 
    int n;

/*==================================================================*/
/*  I)        Usage                                                 */

    fprintf(stdout, "\n     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fprintf(stdout,   "     @                                                @\n");
    fprintf(stdout,   "     @              DO_BAYESIAN_ANALYSIS              @\n");
    fprintf(stdout,   "     @                                                @\n");
    fprintf(stdout,   "     @                                                @\n");
    fprintf(stdout,   "     @             by  Giuseppe Brancato              @\n");
    fprintf(stdout,   "     @                                                @\n");
    fprintf(stdout,   "     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n");
    fprintf(stdout,   "  Usage:                                               \n");
    fprintf(stdout,   "  do_bayesian_analysis input npoints nbin nskip lag output\n\n");
    fflush(stdout);

    if (argc != 7)
       exit(1);

    return 0;
}

float gasdev(long *idum)
{
        float ran1(long *idum);
        static int iset=0;
        static float gset;
        float fac,rsq,v1,v2;

        if  (iset == 0) {
                do {
                        v1=2.0*ran1(idum)-1.0;
                        v2=2.0*ran1(idum)-1.0;
                        rsq=v1*v1+v2*v2;
                } while (rsq >= 1.0 || rsq == 0.0);
                fac=sqrt(-2.0*log(rsq)/rsq);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
        int j;
        long k;
        static long iy=0;
        static long iv[NTAB];
        float temp;

        if (*idum <= 0 || !iy) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ;
                        *idum=IA*(*idum-k*IQ)-IR*k;
                        if (*idum < 0) *idum += IM;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ;
        *idum=IA*(*idum-k*IQ)-IR*k;
        if (*idum < 0) *idum += IM;
        j=iy/NDIV;
        iy=iv[j];
        iv[j] = *idum;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/********************************************************************/
/*                                                                  */
/*                  ROUTINE: ConjugateGradients                     */
/*                                                                  */
/********************************************************************/

int ConjugateGradients (int  nbin,
                        real lagtime,
                        real *Ptheo,
                        real **Rij,
                        real **M_t0)

{
/*==================================================================*/
/*            Local variable declarations                           */

    real  *p, fval;
    real  fret;
    real  temp;
    int   i, j, n, iter;

/*==================================================================*/
/*  I)        Initialize local variables                            */

    SA_func_ndof = nbin;
    SAfunc_parm1 = Ptheo;
    SAfunc_parm2 = lagtime;
    SAfunc_parm3 = M_t0;

    p = AllocVecR(SA_func_ndof);

/*==================================================================*/
/*  II)       Initialize coordinate vector, compute initial energy  */
/*            and write the initial configuration                   */

    for (i = 1; i <= SA_func_ndof - 1; i++) {
      p[i] = Rij[i+1][i];
    }
    p[SA_func_ndof] = Rij[1][SA_func_ndof];

    fval = func (p);

/*==================================================================*/
/*  III)      Print on logfile                                      */

    PrintSepIni(stdout);
    fprintf(stdout,"  Conjugate Gradients Minimization\n");
    PrintSep(stdout);
    fprintf(stdout," Nr. Degree of Freedom: %10d  \n",  SA_func_ndof);
    fprintf(stdout," Starting Energy      : %10.5f\n",  fval);
    fflush(stdout);
/*==================================================================*/
/*  IV)       Performe Conjugate Gradients Minimization             */

    frprmn(p, SA_func_ndof, 1.0e-14, &iter, &fret, func, dfunc);

/*==================================================================*/
/*  V)        Print on logfile and write final configuration        */
    fprintf(stdout," Final Energy         : %10.5f\n",  fret);
    fprintf(stdout," Number of iterations : %10d\n", iter);
    PrintSepEnd(stdout);



    for (i = 1; i <= SA_func_ndof; i++) {
      for (j = 1; j <= SA_func_ndof; j++) {
        Rij[i][j] = 0.0;
      }
    }
    for (i = 1; i <= SA_func_ndof - 1; i++) {
      Rij[i+1][i] = p[i];
    }
    Rij[1][SA_func_ndof] = p[SA_func_ndof];

    /* Further check */
    for (i = 1; i <= SA_func_ndof - 1; i++) {
      Rij[i][i+1] = Rij[i+1][i] * SAfunc_parm1[i] / SAfunc_parm1[i+1];
    }
    Rij[SA_func_ndof][1] = Rij[1][SA_func_ndof] * SAfunc_parm1[SA_func_ndof] / SAfunc_parm1[1];

    for (i = 1; i <= SA_func_ndof; i++) {
      temp = 0.0;
      for (j = 1; j <= SA_func_ndof; j++) {
        if (j != i)
          temp = temp - Rij[j][i];
      }
      Rij[i][i] = temp;
    }

    FreeVecR(p);

    return 0;
}
/*----------------- END ROUTINE: ConjugateGradients ----------------*/



/********************************************************************/
/*                                                                  */
/*                  ROUTINE: SimulatedAnnealing                     */
/*                                                                  */
/********************************************************************/

/* The value that causes the "set/query" functions to just query */
#define NO_VALUE        -1 
#define NO_VALUE_INT    -1 
#define NO_VALUE_REAL   -1.0 

/* #define MELT_ONLY  */
/* #define EQUIL_ONLY   */
/* #define MAXIT 1000 */
/* real lambda = 0.0; */            /* derivative weight factor */

int SimulatedAnnealing (int  nbin,
                        real lagtime,
                        real *Ptheo,
                        real **Rij,
                        real **M_t0)
{
/*==================================================================*/
/*            Local variable declarations                           */

    real    TempIni       = 0.05;
    real    Tscale        = 0.01;
    real    Boltzmann     = 1.0;
    real    Learning_rate = 0.500;
    real    Dwell         = 1000;
    real    Range         = PI / 4.0;
    int     Niter         = 2000;

    real    TempEnd, TempEnd2;
    real    *p, fval;
    real    temp;
    int     i, j, n;

/*==================================================================*/
/*  I)        Initialize local variables                            */

    SA_func_ndof = nbin;
    SAfunc_parm1 = Ptheo;
    SAfunc_parm2 = lagtime;
    SAfunc_parm3 = M_t0;

    p = AllocVecR(SA_func_ndof);

/*==================================================================*/
/*  II)       Initialize coordinate vector, compute initial energy  */
/*            and write the initial configuration                   */

    for (i = 1; i <= SA_func_ndof - 1; i++) {
      p[i] = Rij[i+1][i];
    }
    p[SA_func_ndof] = Rij[1][SA_func_ndof];

    fval = func (p);
/*
    WriteXVF (step, Time, FileName, x, v, f, Gen);
*/

/*==================================================================*/
/*  III)      Initialize Monte Carlo parameters                     */

    if ( SAInit(costf, SA_func_ndof) == 0 ) {
      fprintf(stderr,"Trouble initializing in SAInit\n");
      exit(1);
    }

    SAinitial (p, fval);                /* set initial condition */
    SAtemperature (TempIni);
    SABoltzmann (Boltzmann);
    SAlearning_rate (Learning_rate); /* SAlearning_rate (0.0005); */
    SAtscale (Tscale);
    SAdwell  (Dwell);
    SArange  (Range); /* SArange  (PI/4.0); */
    SAiterations (Niter);

    /* Estimated Final Temperature */
    TempEnd = TempIni / (1.0 + (Niter - 1) * SAtscale(NO_VALUE_REAL));

/*==================================================================*/
/*  IV)       Print on logfile                                      */

    PrintSepIni(stdout);
    fprintf(stdout,"  Monte Carlo Simulated Annealing \n");
    PrintSep(stdout);
    fprintf(stdout," Nr. Iterations       : %10d  \n",  SAiterations(NO_VALUE));
    fprintf(stdout," Nr. Equi. Iterations : %10d  \n",  SAdwell(NO_VALUE));
    fprintf(stdout," Starting Temperature : %10.3f  K\n",  SAtemperature(NO_VALUE_REAL));
    fprintf(stdout," Final Temperature    : %10.3f  K\n",  TempEnd);
    fprintf(stdout," Temper. Scale Factor : %10.3f\n",  SAtscale(NO_VALUE_REAL));
    fprintf(stdout," Nr. Degree of Freedom: %10d  \n",  SA_func_ndof);
    PrintSep(stdout);
    fprintf(stdout," Starting Energy      : %15.8f\n",  fval);
    PrintSep(stdout);
    fflush(stdout);

/*==================================================================*/
/*  V)        Performe Monte Carlo Simulated Annealing              */

    /* SAmelt(NO_VALUE); */       /* melt the system  */
#ifndef MELT_ONLY 
        /* make it a bit warmer than melting temperature */
        /* TempIni = 1.2 * SAtemperature(NO_VALUE); TempIni = 300.0;  
           SAtemperature(TempIni); */
#ifndef EQUIL_ONLY       
    TempEnd2 = SAanneal(NO_VALUE);
    SAoptimum(p);
    fval = func (p);
#endif /* EQUIL_ONLY */
#endif /* MELT_ONLY  */
    SAfree();

/*==================================================================*/
/*  VI)       Print on stdout and write the final configuration    */

    PrintSep(stdout);
    fprintf(stdout," Final Energy         : %15.8f  KJ/mol\n", fval);
    if (TempEnd2 != TempEnd) {
      fprintf(stdout,"\n Warning:\n");
      fprintf(stdout," Final temperature %10.4f K doesn't match the estimated final temperature %10.4f K\n", TempEnd2, TempEnd);
    }
    PrintSepEnd(stdout);
    fclose(stdout);




    for (i = 1; i <= SA_func_ndof; i++) {
      for (j = 1; j <= SA_func_ndof; j++) {
        Rij[i][j] = 0.0;
      }
    }
    for (i = 1; i <= SA_func_ndof - 1; i++) {
      Rij[i+1][i] = p[i];
    }
    Rij[1][SA_func_ndof] = p[SA_func_ndof];

    /* Further check */
    for (i = 1; i <= SA_func_ndof - 1; i++) {
      Rij[i][i+1] = Rij[i+1][i] * SAfunc_parm1[i] / SAfunc_parm1[i+1];
    }
    Rij[SA_func_ndof][1] = Rij[1][SA_func_ndof] * SAfunc_parm1[SA_func_ndof] / SAfunc_parm1[1];

    for (i = 1; i <= SA_func_ndof; i++) {
      temp = 0.0;
      for (j = 1; j <= SA_func_ndof; j++) {
        if (j != i)
          temp = temp - Rij[j][i];
      }
      Rij[i][i] = temp;
    }



    FreeVecR(p);
/*
    WriteXVF (Input->nstep, Time, FileName, x, v, f, Gen);
*/
    return 0;
}
/*----------------- END ROUTINE: SimulatedAnnealing  ---------------*/



real costf(real* p)   /* the cost function */
{
  real lambda = 0.0;
  return func (p) + lambda;
}


real func(real p[])
{
    int  ndof = SA_func_ndof;
    real fval;

    real **Rij, **Rij_sym;
    real *Lambda, **V, **VT, **V2, **Diag, **Diag2;
    real **M_t0, **M_t0_sym;
    real temp;
    int  n, d, i, j;

    Rij      = AllocMatR(ndof,ndof);
    Rij_sym  = AllocMatR(ndof,ndof);
    M_t0     = AllocMatR(ndof,ndof);
    M_t0_sym = AllocMatR(ndof,ndof);
    V2       = AllocMatR(ndof,ndof);
    V        = AllocMatR(ndof,ndof);
    VT       = AllocMatR(ndof,ndof);
    Diag     = AllocMatR(ndof,ndof);
    Diag2    = AllocMatR(ndof,ndof);
    Lambda   = AllocVecR(ndof);

    for (i = 1; i <= ndof; i++) {
      for (j = 1; j <= ndof; j++) {
        Diag[i][j]  = 0.0;
        Diag2[i][j] = 0.0;
        Rij[i][j]   = 0.0;
      }
    }
    for (i = 1; i <= ndof - 1; i++) {
      Rij[i+1][i] = p[i];
/*
      fprintf(stdout,"P %5d  %15.8lf\n", i, p[i]);
*/
    }
    Rij[1][ndof] = p[ndof];
/*
    fprintf(stdout,"P %5d  %15.8lf\n", ndof, p[ndof]);
*/
    /* Further check */
    for (i = 1; i <= ndof - 1; i++) {
      Rij[i][i+1] = Rij[i+1][i] * SAfunc_parm1[i] / SAfunc_parm1[i+1];
    }
    Rij[ndof][1] = Rij[1][ndof] * SAfunc_parm1[ndof] / SAfunc_parm1[1];

    for (i = 1; i <= ndof; i++) {
      temp = 0.0;
      for (j = 1; j <= ndof; j++) {
        if (j != i)
          temp = temp - Rij[j][i];
      }
      Rij[i][i] = temp;
    }

/*
    fprintf(fout,"\n\n\n");
    fprintf(fout," R theo\n");
    for (i = 1; i <= nbin; i++) {
      for (j = 1; j <= nbin; j++) {
        fprintf(fout,"%12.8lf", Rij_theo[i][j]);
      }
      fprintf(fout,"\n");
    }
*/
    for (i = 1; i <= ndof; i++) {
      for (j = 1; j <= ndof; j++) {
        Rij_sym[i][j] = (1.0/sqrt(SAfunc_parm1[i])) * Rij[i][j] * sqrt(SAfunc_parm1[j]);
      }
    }

/*
    fprintf(stdout,"\n\n\n");
    fprintf(stdout," R sym\n");
    for (i = 1; i <= ndof; i++) {
      for (j = 1; j <= ndof; j++) {
        fprintf(stdout,"%12.8lf", Rij_sym[i][j]);
      }
      fprintf(stdout,"\n");
    }
*/
    DiagonalMatrix(ndof, Rij_sym, VT, Lambda);
    TransposeMatrix(ndof, VT, V);

    for (i = 1; i <= ndof; i++) {
        Diag[i][i] = exp(SAfunc_parm2 * Lambda[i]);
    }


    ProductMatrix(ndof, V, Diag, V2);
    ProductMatrix(ndof, V2, VT, M_t0_sym);

    for (i = 1; i <= ndof; i++) {
      Diag[i][i]  = sqrt(SAfunc_parm1[i]);
      Diag2[i][i] = 1.0 / sqrt(SAfunc_parm1[i]);
    }

    ProductMatrix(ndof, Diag, M_t0_sym, V2);
    ProductMatrix(ndof, V2, Diag2, M_t0);

/*
    fprintf(stdout,"\n\n\n");
    fprintf(stdout," Mt0 \n");
    for (i = 1; i <= ndof; i++) {
      for (j = 1; j <= ndof; j++) {
        fprintf(stdout,"%12.8lf", M_t0[i][j]);
      }
      fprintf(stdout,"\n");
    }
*/
    fval = 0.0;
    for (i = 1; i <= ndof; i++) {
      for (j = 1; j <= ndof; j++) {
        fval = fval + 1000.0*(M_t0[i][j] - SAfunc_parm3[i][j]) * (M_t0[i][j] - SAfunc_parm3[i][j]);
      }
    }

    FreeMatR(Rij);      
    FreeMatR(Rij_sym);  
    FreeMatR(M_t0);     
    FreeMatR(M_t0_sym); 
    FreeMatR(V2);       
    FreeMatR(V);        
    FreeMatR(VT);       
    FreeMatR(Diag);     
    FreeMatR(Diag2);    
    FreeVecR(Lambda);   

    return fval;
}

void dfunc(real p[], real xi[])
{
  int  ndof = SA_func_ndof;
  real *pt, delta, fval1, fval2;
  int  n, i;

  delta = 0.001;
  pt = AllocVecR(ndof);

  for (n = 1; n <= ndof; n++) {
    for (i = 1; i <= ndof; i++)
      pt[i] = p[i];

      pt[n] = p[n] - delta;
      fval1 = func(pt);
      pt[n] = p[n] + delta;
      fval2 = func(pt);
      xi[n] = (fval2 - fval1)/(2.0*delta);
  }

  FreeVecR(pt);
}

void force_x_max_rms(int ndof, real *p, real *fmax, real *frms,
                     real *xmax, real *xrms, int First)
{
  static real *pp;
  int n, d, i;

  *fmax=0.0;
  *frms=0.0;
  *xmax=0.0;
  *xrms=0.0;

  if (First == 1)  {
    pp = AllocVecR(ndof);

    for (n = 1; n <= ndof; n++) {
      pp[n] = p[n];
    }

  } else if (First == 2) {
    i = 0;

    for (n = 1; n <= ndof; n++) {

/*
          if (*fmax < fabs(f_l[n][d]))
            *fmax = fabs(f_l[n][d]);
          *frms += f_l[n][d] * f_l[n][d];
*/
          if (*xmax < fabs(p[i]-pp[i]))
            *xmax = fabs(p[i]-pp[i]);
          *xrms += (p[i]-pp[i]) * (p[i]-pp[i]);

    }

    *frms = sqrt(*frms / ndof);
    *xrms = sqrt(*xrms / ndof);

    FreeVecR(pp);
  }

}




/********************************************************************/
/*                                                                  */
/*                  ROUTINE: ExponentialMatrix                      */
/*                                                                  */
/********************************************************************/

int BuildRtheo2(int nbin, real **Rij_theo, real *bin)
{
/*==================================================================*/
/*            Local variable declarations                           */

    real *D, *F;
    real delta2, Delta2;
    real xc;
    int  i, j, n;

/*==================================================================*/
/*  I)        Compute next step                                     */

    D = AllocVecR(nbin);
    F = AllocVecR(nbin);

    delta2 = 2.0 * 3.141592653589793 / ((real) nbin);
    delta2 = delta2 * delta2;
    Delta2 = 2.0 * delta2;

    for (n = 1; n <= nbin; n++) {
      xc = (bin[n+1] + bin[n]) / 2.0;
      F[n] = (-cos(2*xc)+2.0);
      D[n] = 0.1 * (2.0 + sin(xc));
    }

    for (i = 1; i <= nbin - 1; i++) {

      Rij_theo[i+1][i] = -1.0*(F[i+1] - F[i]) / Delta2 + (D[i] + D[i+1]) / Delta2;
      Rij_theo[i][i+1] =  1.0*(F[i+1] - F[i]) / Delta2 + (D[i] + D[i+1]) / Delta2;


      Rij_theo[i][i] = D[i+1]*F[i+1]/Delta2 + 2.0*D[i]*F[i]/Delta2 + D[i-1]*F[i-1]/Delta2 
                       -1.0*F[i]*(D[i+1]+D[i-1])/Delta2 -1.0*D[i]*(F[i+1]+F[i-1])/Delta2 +
                       F[i+1]*D[i]/delta2 - 2.0*D[i]*F[i]/delta2 + F[i-1]*D[i]/delta2 +
                       2.0*F[i]/Delta2 - F[i+1]/Delta2 - F[i-1]/Delta2 - 2.0*D[i]/Delta2 - D[i+1]/Delta2 - D[i-1]/Delta2;
    }

    Rij_theo[1][nbin] = -1.0*(F[1] - F[nbin]) / Delta2 + (D[nbin] + D[1]) / Delta2;
    Rij_theo[nbin][1] =  1.0*(F[1] - F[nbin]) / Delta2 + (D[nbin] + D[1]) / Delta2;

    i = nbin;
    Rij_theo[i][i] = D[1]*F[1]/Delta2 + 2.0*D[i]*F[i]/Delta2 + D[i-1]*F[i-1]/Delta2
                       -1.0*F[i]*(D[1]+D[i-1])/Delta2 -1.0*D[i]*(F[1]+F[i-1])/Delta2 +
                       F[1]*D[i]/delta2 - 2.0*D[i]*F[i]/delta2 + F[i-1]*D[i]/delta2 +
                       2.0*F[i]/Delta2 - F[1]/Delta2 - F[i-1]/Delta2 - 2.0*D[i]/Delta2 - D[1]/Delta2 - D[i-1]/Delta2;
    i = 1;
    Rij_theo[i][i] = D[i+1]*F[i+1]/Delta2 + 2.0*D[i]*F[i]/Delta2 + D[nbin]*F[nbin]/Delta2
                       -1.0*F[i]*(D[i+1]+D[nbin])/Delta2 -1.0*D[i]*(F[i+1]+F[nbin])/Delta2 +
                       F[i+1]*D[i]/delta2 - 2.0*D[i]*F[i]/delta2 + F[nbin]*D[i]/delta2 +
                       2.0*F[i]/Delta2 - F[i+1]/Delta2 - F[nbin]/Delta2 - 2.0*D[i]/Delta2 - D[i+1]/Delta2 - D[nbin]/Delta2;


    FreeVecR(F);
    FreeVecR(D);

    return 0;
}



/********************************************************************/
/*                                                                  */
/*                  ROUTINE: ExponentialMatrix                      */
/*                                                                  */
/********************************************************************/

int BuildRtheo3(int nbin, real **Rij_theo, real *bin)
{
/*==================================================================*/
/*            Local variable declarations                           */

    real *D, *P, Ptot;
    real delta, delta2;
    real xc;
    int  i, j, n;
    int  m1, m2, p1, p2;

/*==================================================================*/
/*  I)        Compute next step                                     */

    D = AllocVecR(nbin);
    P = AllocVecR(nbin);

    delta  = 2.0 * 3.141592653589793 / ((real) nbin);
    delta2 = 12.0 * delta * delta;

    /* Theoretical values */
    for (i = 1; i <= nbin; i++) {
      xc = (bin[i+1] + bin[i]) / 2.0;
      P[i] = delta * exp(-1.0*(-cos(2*xc)+2.0));
      D[i] = 0.1 * (2.0 + sin(xc));
    }

    Ptot = 0.0;
    for (i = 1; i <= nbin; i++)
      Ptot = Ptot + P[i];

    for (i = 1; i <= nbin; i++)
      P[i] = P[i] / Ptot;
    /* End theoretical values */


    for (i = 1; i <= nbin; i++) {
      p1 = i + 1;
      p2 = i + 2;
      m1 = i - 1;
      m2 = i - 2;
      if (p1 > nbin) p1 = p1 - nbin;
      if (p2 > nbin) p2 = p2 - nbin;
      if (m1 < 1) m1 = m1 + nbin;
      if (m2 < 1) m2 = m2 + nbin;


      Rij_theo[i][p2] = -1.0  * D[p2] / delta2;
      Rij_theo[i][p1] =  16.0 * D[p1] / delta2;
      Rij_theo[i][m2] = -1.0  * D[m2] / delta2;
      Rij_theo[i][m1] =  16.0 * D[m1] / delta2;


      Rij_theo[i][i] = D[p2] * (P[p2] / P[i]) / delta2  - 16.0 * D[p1] * (P[p1] / P[i]) / delta2 
                       - 16.0 * D[m1] * (P[m1] / P[i]) / delta2 + D[m2] * (P[m2] / P[i]) / delta2;

    }

    FreeVecR(P);
    FreeVecR(D);

    return 0;
}


/********************************************************************/
/*                                                                  */
/*                  ROUTINE: ExponentialMatrix                      */
/*                                                                  */
/********************************************************************/

int BuildRtheo(int nbin, real **Rij_theo, real *bin)
{
/*==================================================================*/
/*            Local variable declarations                           */

    real *D, *P, *D1, *P1, Ptot;
    real delta, delta2;
    real xc;
    int  i, j, n;
    int  m1, m2, p1, p2;

/*==================================================================*/
/*  I)        Compute next step                                     */

    D  = AllocVecR(nbin);
    P  = AllocVecR(nbin);
    D1 = AllocVecR(nbin);
    P1 = AllocVecR(nbin);

    delta  = 2.0 * 3.141592653589793 / ((real) nbin);
    delta2 = 12.0 * delta * delta;

    /* Theoretical values */
    for (i = 1; i <= nbin; i++) {
      xc = (bin[i+1] + bin[i]) / 2.0;
      P[i]  = -1.0*cos(2*xc) + 2.0;
      P1[i] = 2.0*sin(2*xc);
      D[i]  = 0.1 * (2.0 + sin(xc));
      D1[i] = 0.1 * cos(xc);
    }


    for (i = 1; i <= nbin; i++) {
      p1 = i + 1;
      p2 = i + 2;
      m1 = i - 1;
      m2 = i - 2;
      if (p1 > nbin) p1 = p1 - nbin;
      if (p2 > nbin) p2 = p2 - nbin;
      if (m1 < 1) m1 = m1 + nbin;
      if (m2 < 1) m2 = m2 + nbin;


      Rij_theo[i][p2] = -1.0 * (D[i] * P[i] + D1[i]) / delta2 - 1.0 * D[i] / delta2;
      Rij_theo[i][p1] =  8.0 * (D[i] * P[i] + D1[i]) / delta2 + 16.0 * D[i] / delta2;
      Rij_theo[i][m1] = -8.0 * (D[i] * P[i] + D1[i]) / delta2 + 16.0 * D[i] / delta2;
      Rij_theo[i][m2] =  1.0 * (D[i] * P[i] + D1[i]) / delta2 - 1.0 * D[i] / delta2;



      Rij_theo[i][i] = D1[i] * P[i] + D[i] * P1[i] - 30.0 * D[i] / delta2;

    }

    FreeVecR(P);
    FreeVecR(D);
    FreeVecR(P1);
    FreeVecR(D1);

    return 0;
}



