/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"

void jacobi_col_major ( double * AB, double * RHS, int *lab, int * la, int *kv, int it_max, double seuil)	{

	double * x0 = (double * ) calloc ( (*la), sizeof(double)),
		   * x1 = (double * ) calloc ( (*la), sizeof(double));/*,
		   * relres = (double *) malloc ( it_max * sizeof(double));*/
		   
	x0 = (double * ) calloc ( (*la), sizeof(double));
	x1 = (double * ) calloc ( (*la), sizeof(double));

	/* relres = (double *) malloc ( it_max * sizeof(double));*/	   
	for (int i = 0; i < it_max; i++)	{
		x1[0] = ((AB[*kv+2]*x0[1])+RHS[0])/AB[*kv+1];
		for (int jj=1;jj<(*la)-1;jj++){
			int kk = jj*(*lab);
			x1[jj] = (-((-AB[kk+ *kv])*x0[jj-1])-(-(AB[kk+ *kv+2])*x0[jj+1])+RHS[jj])/AB[kk+ *kv+1];
		}
		x1[(*la)-1] = ((AB[(*lab)*((*la)-1)+ *kv]*x0[(*la)-2])+RHS[(*lab)-1])/AB[(*lab)*((*la)-1)+ *kv+1];
		
		for (int jj = 0; jj < (*la); jj++)	{
			x0[jj] = x1[jj];
		}
	}
	write_vec ( x1, la, "/tmp/output.txt");
	
	for (int i = 0; i < *la; i++)
		printf ("%lf %lf \n", x0[i], x1[i]);

	free (x0);
	free (x1);
//	free (relres);
}

int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double *RHS, *EX_SOL, *X, *RHS_BLAS;
  double *AB;

  double temp, relres;

  NRHS=1;
  nbpoints=102;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  RHS_BLAS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);

  info=0;

  /* working array for pivot used by LU Factorization */
  ipiv = (int *) calloc(la, sizeof(int));

  int row = 0; //

  if (row == 1){ // LAPACK_ROW_MAJOR
    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la);
    write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");
    
  
  } 
  else { // LAPACK_COL_MAJOR
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");

	jacobi_col_major ( AB, RHS, &lab, &la, &kv, 1000, 0);
  }    
  
  
  printf("\n INFO DGBSV = %d\n",info);

  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative residual */
  temp = cblas_ddot(la, RHS, 1, RHS,1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1);
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  relres = sqrt(relres);
  relres = relres / temp;
  
  printf("\nThe relative residual error is relres = %e\n",relres);

  free(RHS);
  free(RHS_BLAS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(ipiv);

  printf("\n\n--------- End -----------\n");
}
