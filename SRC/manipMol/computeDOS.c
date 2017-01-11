#include "../defMacro.h"
#include "../huckel/header.h"
#include "../algebra_real/algebra.h"
#include "../huckel/huckel.h"
#include "../manag_output/print_screen.h"



void computeDOS_green(double *VP, int nb_orb,double Ef, char *filename)
{

	double minE, maxE,E,dE;
	int nbE;
	double gamma;
	int iE,iorb;
	complex double G;
	double DOS;
	FILE *f;
	
	minE = -5;
	maxE = 10;
	nbE = 1000;
	dE = (maxE-minE)/nbE;
	gamma = 0.05;
	
	
	f = fopen(filename,"w");
  if(!f) 
  {
      printf("couldn't read %s\n",filename);
      exit(1);
  }
	
	for(iE=0;iE<nbE;iE++)
	{

	    E = Ef + minE + iE*dE;
			G = (complex double) 0.;
			for(iorb=0;iorb<nb_orb;iorb++)
				//G += (complex double) 1.0/(E-VP[iorb] + Ef + I*gamma);
				 G += (complex double) 1.0/(E-VP[iorb] + I*gamma);
			DOS = -cimagf(G);
			fprintf(f,"%lf %lf\n",E-Ef,DOS);
	}

  fclose(f);

}
