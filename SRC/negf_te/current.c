#include "../defMacro.h"


//////////////////////////////////////////////////////////
// 		compute the fermi dirac distrib 
//////////////////////////////////////////////////////////
double fermi(double EF,double E,double T)
{
  double f;
  f = 1./(expf( (E-EF)/KB/T)+1.);
  return(f);
}


//////////////////////////////////////////////////////////
// 		integrate between two point 
//////////////////////////////////////////////////////////
double integrate(double *F, double *X, int ind1, int ind2)
{
  int i;
  double curr = 0;
  double K;
  int ind_min,ind_max;
  
  if(ind1<ind2)
  {
    ind_min = ind1;
    ind_max = ind2;
    K = -0.5;
  }
  else
    {
    ind_min = ind2;
    ind_max = ind1;
    K = 0.5;
  }
   
  if((ind1+1)!=ind2)
  {
    for(i=ind_min+1;i<ind_max;i++)
      curr += K*(X[i]-X[i-1])*(F[i]+F[i-1]);      
    
  }
  else
    curr = 0.0;

  //printf("curr = %lf\n",curr);
  return(curr);
  
}



//////////////////////////////////////////////////////////
// 		compute the current 
//////////////////////////////////////////////////////////
void computeCurrent(double *CURR, double *Vbias, double *TE, double *E,  int nb_nrj, double TEMP, double EF)
{
  int i;
  
  double dEMin = 10000000;
  int indexEF;
  double dE;
  double curr;
  double e_h = 0.486 * 0.001;
  

  // locate the FErmi energy of the electrode
  for(i=0;i<nb_nrj;i++)
  {
    dE = fabs(E[i]-EF);
    if(dE<dEMin)
    {
      indexEF = i;
      dEMin = dE;
    }
  }
  
  
  //for all the energies
  for(i=0;i<nb_nrj;i++)
  {
    Vbias[i] = E[i]-EF;
    if(i<indexEF-1 || i>indexEF+1)
      curr= integrate(TE, E, i, indexEF);
    else
      curr = 0.;
    CURR[i] = e_h*curr;
  } 
}

//////////////////////////////////////////////////////////
// 		compute the dI/dV
//////////////////////////////////////////////////////////
void computedIdV(double *dIdV, double *IV, double *E, int nb_nrj)
{
  
  int i;
  for(i=1;i<nb_nrj-1;i++)
    dIdV[i] = (IV[i+1]-IV[i-1])/(E[i+1]-E[i-1]);
  
  dIdV[0] = dIdV[1];
  dIdV[nb_nrj-1] = dIdV[nb_nrj-2];
     
}
