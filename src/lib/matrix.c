/*****************************************************************************************\
*                                                                                         *
*  Matrix inversion, determinants, and linear equations via LU-decomposition              *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  April 2007                                                                    *
*  Mod   :  June 2008 -- Added TDG's and Cubic Spline to enable snakes and curves         *
*           Dec 2008 -- Refined TDG's and cubic splines to Decompose/Solve paradigm       *
*                                                                                         *
\*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "gene_core.h"
#include "matrix.h"

#define TINY 1.0e-20

/****************************************************************************************
 *                                                                                      *
 *  LU-FACTORIZATION SYSTEM SOLVER                                                      *
 *                                                                                      *
 ****************************************************************************************/


//  M is a square double matrix where the row index moves the fastest.
//    LU_Decompose takes M and produces an LU factorization of M that
//    can then be used to rapidly solve the system for given right hand sides
//    and to compute M's determinant.  The return value is NULL if the matrix
//    is nonsingular.  If the matrix appears unstable (had to use a very nearly
//    zero pivot) then the integer pointed at by stable will be zero, and
//    non-zero otherwise.  M is subsumed and effectively destroyed by the routine.

LU_Factor *LU_Decompose(Double_Matrix *M, int *stable)
{ LU_Factor      *F;
  int             n, i, j;
  int            *p, sign;
  double         *v;
  double         *avec[1001], **a;

  n = M->n;

  if (n > 1000)
    a = Malloc(sizeof(double)*n,"Allocating LU Factor work space");
  else
    a = avec;
  F = Malloc(sizeof(LU_Factor),"Allocating LU Factor");
  p = Malloc((sizeof(int) + sizeof(double))*n,"Allocating LU Factor");
  if (a == NULL || F == NULL || p == NULL)
    exit (1);

  v = (double *) (p+n);

  p[0] = 0;
  a[0] = M->m;
  for (i = 1; i < n; i++)
    { a[i] = a[i-1] + n;
      p[i] = i;
    }

  *stable = 1;
  sign    = 1;
  for (i = 0; i < n; i++)  // Find the scale factors for each row in v.
    { double b, f, *r;

      r = a[i];
      b = 0.;
      for (j = 0; j < n; j++)
        { f = fabs(r[j]);
          if (f > b)
            b = f;
        }
      if (b == 0.0)
        { free(p);
          free(F);
          if (n > 1000)
            free(a);
          return (NULL);
        }
      v[i] = 1./b;
    }

  for (j = 0; j < n; j++)      //  For each column
    { double b, s, *r;
      int    k, w;

      for (i = 0; i < j; i++)    // Determine U
        { r = a[i];
          s = r[j];
          for (k = 0; k < i; k++)
            s -= r[k]*a[k][j];
          r[j] = s;
        }

      b = -1.;
      w = j;
      for (i = j; i < n; i++)      // Determine L without dividing by pivot, in order to
        { r = a[i];                //   determine who the pivot should be.
          s = r[j];
          for (k = 0; k < j; k++)
            s -= r[k]*a[k][j];
          r[j] = s;

          s = v[i]*fabs(s);        // Update best pivot seen thus far
          if (s > b)
            { b = s;
              w = i;
            }
	}

      if (w != j)                  // Pivot if necessary
        { r    = a[w];
          a[w] = a[j];
          a[j] = r;
          k    = p[w];
          p[w] = p[j];
          p[j] = k;
          sign = -sign;
          v[w] = v[j];
        }

      if (fabs(a[j][j]) < TINY)    // Complete column of L by dividing by selected pivot
        { if (a[j][j] < 0.)
            a[j][j] = -TINY;
          else
            a[j][j] = TINY;
          *stable = 0;
        }
      b = 1./a[j][j];
      for (i = j+1; i < n; i++)
        a[i][j] *= b;
    }

#ifdef DEBUG_LU
  { int i, j;

    printf("\nLU Decomposition\n");
    for (i = 0; i < n; i++)
      { printf("  %2d: ",p[i]);
        for (j = 0; j < n; j++)
          printf(" %8g",a[i][j]);
        printf("\n");
      }
  }
#endif

  if (n > 1000)
    free(a);

  F->sign   = sign;
  F->perm   = p;
  F->lu_mat = M;
  return (F);
}


//  Display LU factorization F to specified file

void Show_LU_Product(FILE *file, LU_Factor *F)
{ int    n, i, j, k;
  int   *p;
  double u, **a, *d;

  n = F->lu_mat->n;
  d = F->lu_mat->m;
  p = F->perm;
  a = (double **) (p+n);

  for (i = 0; i < n; i++)
    a[i] = d + p[i]*n;

  fprintf(file,"\nLU Product:\n");
  for (i = 0; i < n; i++)
    { for (j = 0; j < i; j++)
        { u = 0.;
          for (k = 0; k <= j; k++)
            u += a[i][k] * a[k][j];
          fprintf(file," %g",u);
        }
      for (j = i; j < n; j++)
        { u = a[i][j];
          for (k = 0; k < i; k++)
            u += a[i][k] * a[k][j];
          fprintf(file," %g",u);
        }
      fprintf(file,"\n");
    }
}


//  Given rhs vector B and LU-factorization F, solve the system of equations
//    and return the result in B.
//  To invert M = L*U given the LU-decomposition, simply call LU_Solve with
//    b = [ 0^k-1 1 0^n-k] to get the k'th column of the inverse matrix.

Double_Vector *LU_Solve(Double_Vector *B, LU_Factor *F)
{ double   *x;
  int       n, i, j;
  int      *p;
  double   *a, *b, s, *r;

  n = F->lu_mat->n;
  a = F->lu_mat->m;
  p = F->perm;
  b = B->m;
  x = (double *) (p+n);

  for (i = 0; i < n; i++)
    { r = a + p[i]*n;
      s = b[p[i]];
      for (j = 0; j < i; j++)
        s -= r[j] * x[j];
      x[i] = s;
    }

  for (i = n; i-- > 0; )
    { r = a + p[i]*n;
      s = x[i]; 
      for (j = i+1; j < n; j++)
        s -= r[j] * b[j];
      b[i] = s/r[i];
    }

  return (B);
}
  

//  Transpose a matrix M in-place and as a convenience return a pointer to it

Double_Matrix *Transpose_Matrix(Double_Matrix *M)
{ int     n;
  double *a;
  int     p, q;
  int     i, j;

  n = M->n;
  a = M->m;

  p = 0;
  for (j = 0; j < n; j++)                 //  Transpose the result
    { q = j;
      for (i = 0; i < j; i++)
        { double x = a[p];
          a[p++] = a[q];
          a[q] = x;
          q += n;
        }
      p += (n-j);
    }

  return (M);
}


//  Generate the right inverse of the matrix that gave rise to the LU factorization f.
//    That is for matrix A, return matrix A^-1 s.t. A * A^-1 = I.  If transpose is non-zero
//    then the transpose of the right inverse is returned.

Double_Matrix *LU_Invert(LU_Factor *F, int transpose)
{ int            n, i, j;
  Double_Matrix *M, G;
  double        *m, *g;

  n = F->lu_mat->n;

  M = Malloc(sizeof(Double_Matrix),"Allocating matrix");
  m = Malloc(sizeof(double)*n*n,"Allocating matrix");
  if (M == NULL || m == NULL)
    exit (1);

  M->n = n;
  M->m = m;
  G.n = n;

  g = m;
  for (i = 0; i < n; i++)                 //  Find the inverse of each column in the
    { G.m = g;
      for (j = 0; j < n; j++)
        g[j] = 0.;
      g[i] = 1.;
      LU_Solve(&G,F);
      g += n;
    }

  if (!transpose)
    Transpose_Matrix(M);

  return (M);
}


//  Given an LU-factorization F, return the value of the determinant of the
//    original matrix.

double LU_Determinant(LU_Factor *F)
{ int     i, n;
  int    *p;
  double *a, det;

  n = F->lu_mat->n;
  a = F->lu_mat->m;
  p = F->perm;

  det = F->sign;
  for (i = 0; i < n; i++)
    det *= a[p[i]*n+i];
  return (det);
}
