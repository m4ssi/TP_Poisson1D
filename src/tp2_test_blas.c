#include "lib_poisson1D.h"
#include <stdio.h>

#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))

void set_GE_operator_poisson1D(double * A, int * la)	{
	int ii, jj, kk;
	for ( ii = 0; ii < (*la); ii++)	{
		for (jj = 0; jj < (*la); jj++)	{
			kk = ii * (*la) +jj;
			A[kk] = 0.0;
			if ( (ii-1) == jj)	A[kk] = -1.0;
			if ( ii == jj)		A[kk] = 2.0;
			if ( (ii+1) == jj)	A[kk] = -1.0;
		}
	}
}

void GE2GB_rowMajor ( double * A, int * la, double * AB, int * lab, int * ku, int * kl)	{
	int ii, jj, kk;
	for (ii = 0; ii < (*la); ii++)	{
		kk = (*kl) - ii;
		for ( jj = max (0, ii-(*kl)); jj < min ((*la), ii+(*ku)+1); jj++)	{
			AB[(kk+jj)+ii*(*lab)] = A[jj + ii*(*la)];
		}
	}
}

void GE2GB_colMajor ( double * A, int * la, double * AB, int * lab, int * ku, int * kl)	{
	int ii, jj, kk;
	for (jj = 0; jj < (*la); jj++)	{
		kk = (*ku) - jj;
		for ( ii = max (0, jj-(*ku)); ii < min ((*la), jj+(*kl)+1); ii++)	{
			AB[(kk+ii)+jj*(*lab)] = A[ii + jj*(*la)];

		}
	}
}

int main ( int argc, char *argv[])	{
	int 	ierr,
			jj,
			nbpoints,
			la, lab,
			ku, kl, kv,
			ldb, info,
			incx, incy,
			NHRS,
			*ipiv;

	double	T0, T1,
			alpha, beta,
			temp, relres,
			*RHS, *EX_SOL, *X,
			*A, *AB, *AB_comp, *XX, *YY;
			 
	incx = 1;
	incy = 1;
	alpha = 1.0;
	beta = 0.0;
	NHRS = 1;
	nbpoints = 102;
	la = nbpoints-2;
	T0 = -5.0;
	T1 = 5.0;
	kv = 1;
	ku = 1;
	kl = 1;
	lab = ku+kl+kv+1;
	
	printf ( "--------- Poisson 1D ---------\n");
	
	RHS = (double*) malloc ( la * sizeof (double));
	EX_SOL = (double*) malloc ( la * sizeof (double));
	X = (double*) malloc ( la * sizeof (double));
	
	set_grid_points_1D ( X, &la);
	set_dense_RHS_DBC_1D ( RHS, &la, &T0, &T1);
	set_analytical_solution_DBC_1D ( EX_SOL, X, &la, &T0, &T1);
	
	write_vec ( RHS, &la, "TEST_BLAS_ROW_MAJOR/RHS.dat");
	write_vec ( EX_SOL, &la, "TEST_BLAS_ROW_MAJOR/EX_SOL.dat");
	write_vec ( X, &la, "TEST_BLAS_ROW_MAJOR/X_grid.dat");
	
	
	A = (double*) malloc ( (la*la) * sizeof(double));
	AB = (double*) malloc ( (la*lab) * sizeof(double));
	AB_comp = (double*) malloc ( (la*lab) * sizeof(double));
	
	XX = (double*) calloc ( la, sizeof(double));
	YY = (double*) calloc ( la, sizeof(double));
	
	ipiv = (int *) calloc ( la, sizeof ( int));
	
	int la2 = la*la,
		labla = lab*la,
		row = 1;
		
	
	if ( row == 1)	{
		set_GE_operator_poisson1D ( A, &la);
		write_vec ( A, &la2, "TEST_BLAS_ROW_MAJOR/A_row_vect.dat");
	
		GE2GB_rowMajor ( A, &la, AB, &lab, &ku, &kl);
		write_GB_operator_rowMajor_poisson1D ( AB, &lab, &la, "TEST_BLAS_ROW_MAJOR/AB_row_vect.dat");
		
		cblas_dgbmv ( CblasRowMajor, CblasNoTrans, la, la, kl, ku, alpha, AB, lab, X, incx, beta, RHS, incy);
		write_vec ( X, &la, "TEST_BLAS_ROW_MAJOR/XX.dat");		
		write_vec ( RHS, &la, "TEST_BLAS_ROW_MAJOR/YY.dat");
	}
	else	{
		set_GE_operator_poisson1D ( A, &la);
		write_vec ( A, &la2, "TEST_BLAS_ROW_MAJOR/A_col_vect.dat");
	
		GE2GB_colMajor ( A, &la, AB, &lab, &ku, &kl);
		write_GB_operator_colMajor_poisson1D ( AB, &lab, &la, "TEST_BLAS_ROW_MAJOR/AB_col_vect.dat");
		
	
		cblas_dgbmv ( CblasColMajor, CblasNoTrans, la, la, kl, ku, alpha, AB, lab, X, incx, beta, RHS, incy);
		write_vec ( X, &la, "TEST_BLAS_ROW_MAJOR/XX.dat");		
		write_vec ( RHS, &la, "TEST_BLAS_ROW_MAJOR/YY.dat");		
		//~ write_vec ( AB, &labla, "TEST_BLAS_ROW_MAJOR/AB_col_sol.dat");
	}
	write_xy( RHS, X, &la, "TEST_BLAS_ROW_MAJOR/SOL.dat");
	free ( RHS);
	free ( EX_SOL);
	free ( X);
	free ( A);
	free ( AB);
	free ( AB_comp);
	free ( XX);
	free ( YY);
	free ( ipiv);
	
	printf("\n\n--------- End -----------\n");	
	
	return 0;
}
