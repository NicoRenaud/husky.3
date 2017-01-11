#include "../defMacro.h"
#include "../manag_output/print_screen.h"

/////////////////////////////////////////////////
//
//	extract a sub matrix from a bigger one
//
/////////////////////////////////////////////////
void extract(double *Hout, double *Hin, int *ind1, int n1, int *ind2, int n2, int nin)
{	
  int II,J,i,j;
  for(II=0;II<n1;II++)
  {
		
    i = ind1[II];
    for(J=0;J<n2;J++)
    {
      j = ind2[J];
      Hout[II*n2+J] = Hin[i*nin+j];
    }
  }
}

////////////////////////////////////////////////////////////////////////////
//  PUT A MATRIX IN A LARGER ONE
////////////////////////////////////////////////////////////////////////////
void insert(double *Hdest, double *Hsource, int *ind1, int n1, int *ind2, int n2, int nin)
{
  int ii,jj,i,j;
  for(ii=0;ii<n1;ii++)
  {
    i = ind1[ii];
    for(jj=0;jj<n2;jj++)
    {
      j = ind2[jj];
      //Hout[I*n2+J] = Hin[i*nin+j];
      Hdest[i*nin+j] = Hsource[ii*n2+jj];
    }
  }
}
