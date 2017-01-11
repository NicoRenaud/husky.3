#include "../defMacro.h"

////////////////////////////////////
//      write a Matrix            //
//      to a file									//
////////////////////////////////////

void export_mat_tofile(double *MAT,int n, char *fname)
{

	FILE *f;
	int i,j;
	
	f = fopen(fname,"w");
  if (!f) 
  {
      printf("couldn't open %s\n",fname);
      exit(1);
  }
	
	for(i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			fprintf(f,"%lf ",MAT[i*n+j]);
		}
		fprintf(f,"\n");
	}

	fclose(f);

}