/*****************************************************************************************\
*                                                                                         *
*  Matrix inversion, determinants, and linear equations via LU-decomposition              *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  April 2007                                                                    *
*  Mod   :  June 2008 -- Added TDG's and Cubic Spline to enable snakes and curves         *
*                                                                                         *
\*****************************************************************************************/

#ifndef _MATRIX_LIB

#define _MATRIX_LIB

typedef struct
  { int     n;
    double *m;
  } Double_Matrix;

typedef Double_Matrix Double_Vector;

typedef struct
  { Double_Matrix *lu_mat;  //  LU decomposion: L is below the diagonal and U is on and above it
    int           *perm;    //  Permutation of the original rows of m due to pivoting
    int            sign;    //  Sign to apply to the determinant due to pivoting (+1 or -1)
  } LU_Factor;

LU_Factor     *LU_Decompose(Double_Matrix *M, int *stable);
void           Show_LU_Factor(FILE *file, LU_Factor *F);
Double_Vector *LU_Solve(Double_Vector *B, LU_Factor *F);
Double_Matrix *Transpose_Matrix(Double_Matrix *M);
Double_Matrix *LU_Invert(LU_Factor *F, int transpose);
double         LU_Determinant(LU_Factor *F);

#endif
