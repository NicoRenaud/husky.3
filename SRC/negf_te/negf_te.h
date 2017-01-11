////////////////////////
// header for negf.c
////////////////////////

#ifndef negf_h
#define negf_h

void dosElec(complex double *Gamma, complex double *Sigma, int size);

void selfElec(complex double *Sigma, double E, complex double *S, complex double *V, complex double *g, int nC, int nX);

void NEGF(complex double *Ga,  complex double *Gr, 
	  double E, double *Smol, double *Hmol, complex double *Sigma_l, complex double *Sigma_r, int size);

void transmission(double *t, double *d, double E, double *Hmol, double *Smol, 
		    complex double *gamal_cplx, complex double *gamar_cplx, 
		    complex double *Sigma_l, complex double *Sigma_r,int size);

void compute_surf_green_func(complex double *GS, double LDOS, int size_elec);

void NEGF_TE(double *TE, double *DOS, double *energies, double LDOS,
	     double *Hmol, double *Smol, 
	     double *Vl, double *Sl,
	     double *Vr, double *Sr, int size_mol, int size_elec_1, int size_elec_2, int nEnergies);

double computeDOS(complex double *Ga, complex double *Gr, int size);

#endif
