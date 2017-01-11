
#include "../defMacro.h"
#include "../algebra_cplx/algebra_cplx.h"
#include "../manag_output/print_screen.h"



///////////////////////////////////////////
// function to comoute the transmission
// using NEGF formalism
///////////////////////////////////////////



//////////////////////////////////////////////////////////////
//	density of states of the electrodes
//////////////////////////////////////////////////////////////
void dosElec(complex double *Gamma, complex double *Sigma, int size)
{
 
 int i;
 for(i=0;i<size*size;i++)
   Gamma[i] = -2*cimagl(Sigma[i]);

}

//////////////////////////////////////////////////////////////
//	self energies of the electrodes
//////////////////////////////////////////////////////////////
void selfElec(complex double *Sigma, double E, complex double *S, complex double *V, complex double *g, int nC, int nX)
{
  int i,j,k=0;
  
  complex double e = (complex double) E;
  
 
  complex double *M, *Mtrans, *TEMP1,*TEMP2;
  
 
  M = calloc(nX*nC,sizeof(complex double));
  Mtrans = calloc(nC*nX,sizeof(complex double));
  TEMP1 = calloc(nC*nX,sizeof(complex double));
  TEMP2 = calloc(nC*nC,sizeof(complex double));
  
  double  a[2] = {1.0,0.0};
  double  b[2] = {0.0,0.0};
  

  // form ES-V
  for(i=0;i<nX;i++)
  {
    for(j=0;j<nC;j++)
    {
      M[k] = e*S[k]-V[k];
      k++;
    }
  }
  
  // transpose the M matrix
  transMatrect(Mtrans, M, nX, nC);
  
  // compute the first temp matrix
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nC,nX,nX,a,Mtrans,nX,g,nX,b,TEMP1,nX);

  // compute the second temp matrix
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,nC,nC,nX,a,TEMP1,nX,M,nC,b,TEMP2,nC);
  
  if(0)
  {
    printf("nc = %d\nnX = %d\n",nC,nX);
    print_mat_C(g,nX,nX,"g");
    print_mat_C(M,nX,nC,"M");
    print_mat_C(Mtrans,nC,nX,"Mtrans");
    print_mat_C(TEMP1,nC,nX,"TMP1");
    print_mat_C(TEMP2,nC,nC,"TMP2");
    exit(1);
  }
  
  // copy the result
   cblas_zcopy(nC*nC,TEMP2,1,Sigma,1);


  // free mem
  free(M);
  free(Mtrans);
  free(TEMP1);
  free(TEMP2);
  

}


//////////////////////////////////////////////////////////////
//	non equilibrium Green function for a give E
//////////////////////////////////////////////////////////////
void NEGF(complex double *Ga,  complex double *Gr, 
	  double E, double *Smol, double *Hmol, complex double *Sigma_l, complex double *Sigma_r, int size)
{
   
  int i,j;
  complex double *ga_temp,*gr_temp;
  
  // alloc memoire pour temp mat
  ga_temp = calloc(size*size,sizeof(complex double));
  gr_temp = calloc(size*size,sizeof(complex double));
  
  // comute temp mat_temp
  for(i=0;i<size;i++)
  {
      for(j=0;j<size;j++)
      {
	  gr_temp[i*size+j] = (complex double) E*Smol[i*size+j] -1.*( (complex double) Hmol[i*size+j] + Sigma_l[i*size+j] + Sigma_r[i*size+j]   );
	  ga_temp[i*size+j] = (complex double) E*Smol[i*size+j] -1.*( (complex double) Hmol[i*size+j] - Sigma_l[i*size+j] - Sigma_r[i*size+j]   );	
	  
      }
  }
   
   
  // inverse the temp matrices 
  invMatComp(Ga, ga_temp, size);
  invMatComp(Gr, gr_temp, size);
   
     
  // free the memeory
  free(ga_temp);
  free(gr_temp);
 
}


//////////////////////////////////////////////////////////////
//	Compute the projected density of state
//////////////////////////////////////////////////////////////
double computeDOS(complex double *Ga, complex double *Gr, int size)
{
  
  int i;
  double D;
  complex double Dcplx;
  complex double N = -I/2./PI;
  complex double *mtemp;
  mtemp = calloc(size*size,sizeof(complex double));
  
  for(i=0;i<size*size;i++)
    mtemp[i] = (Ga[i]-Gr[i])*N;
 
  trace_cplx(&Dcplx, mtemp, size);
  
  D = creall(Dcplx);
  
  free(mtemp);
  
  return(D);
  
  
}


//////////////////////////////////////////////////////////////
//	transmission for a given energy
//////////////////////////////////////////////////////////////
void transmission(double *t, double *d, double E, double *Hmol, double *Smol, 
		    complex double *gamal_cplx, complex double *gamar_cplx, 
		    complex double *Sigma_l, complex double *Sigma_r,int size)
{
  
 double T,D;
 complex double *Ga,*Gr;
 complex double *mtemp1,*mtemp2,*mtemp3;
 double  a[2] = {1.0,0.0};
 double  b[2] = {0.0,0.0};
 complex double tr;
 
 // alloc memeoire
 Ga = calloc(size*size,sizeof(complex double));
 Gr = calloc(size*size,sizeof(complex double));
 mtemp1 = calloc(size*size,sizeof(complex double));
 mtemp2 = calloc(size*size,sizeof(complex double));
 mtemp3 = calloc(size*size,sizeof(complex double));
 
 
 // compute the negf
  NEGF(Ga,  Gr,  E, Smol, Hmol, Sigma_l, Sigma_r, size);
  

  // mtemp1 = GamaL*Gr 
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,
	      a,gamal_cplx,size,Gr,size,b,mtemp1,size);
  
	
  // mtemp2 = GamaR*Ga
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,
	      a,gamar_cplx,size,Ga,size,b,mtemp2,size);
     
  // mtemp3  = GamaL*Gr*GamaR*Ga
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,size,size,size,
	      a,mtemp1,size,mtemp2,size,b,mtemp3,size);
  
  // compute the trace of  GamaL*Gr*GamaR*Ga
  trace_cplx(&tr, mtemp3, size);
 
  // compute the dos
  D = computeDOS(Ga,Gr,size);
  
  // strore the trace
  T = creall(tr);
  
  *t = T;
  *d = D;
  
  // free memory
  free(Ga);
  free(Gr);
  free(mtemp1);
  free(mtemp2);
  free(mtemp3);

}

/////////////////////////////////////////////
//	compute the surface GF WBA
/////////////////////////////////////////////
void compute_surf_green_func(complex double *GS, double LDOS, int size_elec)
{
  
  int i;
  for(i=0;i<size_elec;i++)
    GS[i*(size_elec+1)] = -I * PI *LDOS;
  
}


//////////////////////////////////////////////////////////////
//	compute the TE
//////////////////////////////////////////////////////////////
void NEGF_TE(double *TE, double *DOS, double *energies, double LDOS,
	     double *Hmol, double *Smol, 
	     double *Vl,  double *Sl,
	     double *Vr, double *Sr, int size_mol, int size_elec_1, int size_elec_2, int nEnergies)
{

 double t,d;
 int i;
 complex double *gama_l, *gama_r;
 complex double *sigma_l, *sigma_r;
 complex double *GS_1, *GS_2;
 complex double *vl_cplx, *sl_cplx;
 complex double *vr_cplx, *sr_cplx;
 double percent=0., i_dble=0, thr=0;
 
 int verbose = 1;
 
 
 // alloc memoire
 vl_cplx = calloc(size_elec_1*size_mol,sizeof(complex double));
 sl_cplx = calloc(size_elec_1*size_mol,sizeof(complex double));
 vr_cplx = calloc(size_elec_2*size_mol,sizeof(complex double));
 sr_cplx = calloc(size_elec_2*size_mol,sizeof(complex double));
 
 gama_l = calloc(size_mol*size_mol,sizeof(complex double));
 gama_r = calloc(size_mol*size_mol,sizeof(complex double));
 
 sigma_l = calloc(size_mol*size_mol,sizeof(complex double));
 sigma_r = calloc(size_mol*size_mol,sizeof(complex double));
 
 // surface green function of the electrode
 GS_1 = calloc(size_elec_1*size_elec_1,sizeof(complex double));
 GS_2 = calloc(size_elec_2*size_elec_2,sizeof(complex double));
 
 compute_surf_green_func(GS_1,LDOS,size_elec_1);
 compute_surf_green_func(GS_2,LDOS,size_elec_2);
  
 //cast the matrix
 dbl2cplx(sl_cplx, Sl, size_elec_1*size_mol);
 dbl2cplx(sr_cplx, Sr, size_elec_2*size_mol);
 dbl2cplx(vl_cplx, Vl, size_elec_1*size_mol);
 dbl2cplx(vr_cplx, Vr, size_elec_2*size_mol);
  
 

  // compute the transmission
  for(i=0;i<nEnergies;i++)
  {
        
    // self energies
    selfElec(sigma_l, energies[i], sl_cplx,vl_cplx, GS_1, size_mol, size_elec_1);
    selfElec(sigma_r, energies[i], sr_cplx,vr_cplx, GS_2, size_mol, size_elec_2);

    
    // compute the dos of the electrode    
    dosElec(gama_l,  sigma_l, size_mol);
    dosElec(gama_r,  sigma_r, size_mol);
    
    // compute the transmission and the DOS
    transmission(&t,&d,energies[i],Hmol,Smol,gama_l,gama_r, sigma_l, sigma_r,size_mol);
    
    // store everyhing
    TE[i] = t;
    DOS[i] = d;
    
    // print the progress
    if(verbose)
    {
	 percent =  i_dble/nEnergies;
	 if(percent>thr)
	 {
	   printf(" %3.0f%%", percent*100);
	   fflush(stdout);
	   thr += 0.1;
	 }
	 i_dble = i;
     }
    
    // check if we ve reached a NaN or Inf
    if(isnan(TE[i]) || isinf(TE[i]))
    {
      printf(" ERROR transmission coefficient is not a number\n");
      exit(1);
    }
  }
 
  // done !
  if(verbose)
    printf(" 100%%\n\n"); 
  
  // free memory
  free(gama_l);
  free(gama_r);
  free(sigma_l);
  free(sigma_r);
  free(GS_1);
  free(GS_2);
  free(sl_cplx);
  free(sr_cplx);
  free(vl_cplx);
  free(vr_cplx);
  
}

