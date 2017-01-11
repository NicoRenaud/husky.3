#include "../defMacro.h"
#include "./algebra.h"
#include "../manag_output/print_screen.h"

///////////////////////////////////////////////
// Spec pencil total routine
///////////////////////////////////////////////
void spec_pencil_zinger(double *VAL_PRP, double *VECT_PRP,double *H, double *S, int nb_orb)
{
	// temporary eigenvalues
 	double *val_prp;
	val_prp = calloc(nb_orb,sizeof(double));


	//diagonalize with cblas_dggev
	printf("\n   === start diagonalization of H using cblas_dggev\n");
	spec_pencil(VECT_PRP, val_prp, H, S, nb_orb);
	
	// copy the eigenvalues for sorting
	cblas_dcopy(nb_orb,val_prp,1,VAL_PRP,1);
	
	// sort the eigenvalues/eigenvectors
	qsort(VAL_PRP, nb_orb, sizeof(double), compare_doubles);
	reorder_vect_prp(VECT_PRP,VAL_PRP,val_prp,nb_orb);
	
	// free the memory
	free(val_prp);
	
}


///////////////////////////////////////////////
// Spec SCHMIDT routine
///////////////////////////////////////////////
void spec_shmidt_zinger(double *VAL_PRP, double *VECT_PRP, double *H, double *Vs, double *Us, double *sqrtVs, double *H0, double *S, int nb_orb)
{

	int i;
	double *val_prp;
	double *vect_prp;
	int debug = 0;
	
	// timer
	time_t t1;

	if(debug)
		printf("\n   === start diagonalization of H using cblas_dsyev\n");
			
	// memory allocation
	//H = calloc(nb_orb*nb_orb,sizeof(double));
	val_prp = calloc(nb_orb,sizeof(double));
	vect_prp = calloc(nb_orb*nb_orb,sizeof(double));
	
	// copy the initial Hamiltonian
	cblas_dcopy(nb_orb*nb_orb, H0, 1, H, 1);


	// diagonalize the overlap
	if(debug)
		printf("\t -- diagonalize the overlap matrix\n");
	spec(Us, Vs, S, nb_orb);
	
	// compute 1/S^{1/2}
	if(debug)
		printf("\t -- form S^(-1/2)\n");
	for(i=0;i<nb_orb;i++)
		sqrtVs[i*nb_orb+i] = 1./sqrtf(Vs[i]);				

	
	// form the rotation matrix
	if(debug)
		printf("\t -- form the orthogonalization matrix\n");
	real_eigen2local(sqrtVs, Us, nb_orb);

	// rotate the Hamiltonian
	if(debug)
		printf("\t -- orthogonalize the Hamiltonian\n");
	real_rotateSchmidt(H, sqrtVs, nb_orb);

	// diagonalize the hamiltonian
	if(debug)
		printf("\t -- diagonalize the (orhtogonalized) Hamiltonian\n");
	t1 = clock();
	spec(VECT_PRP,val_prp,H,nb_orb);
	t1 = (int) clock()-t1;
	
	
	// copy the val prp and reorder the spectrum
	cblas_dcopy(nb_orb,val_prp,1,VAL_PRP,1);
	qsort(VAL_PRP, nb_orb, sizeof(double), compare_doubles);
	reorder_vect_prp(VECT_PRP,VAL_PRP,val_prp,nb_orb);
	

	// free local memory
	free(val_prp);
	free(vect_prp);
	
}
