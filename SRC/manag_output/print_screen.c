
#include "../defMacro.h"

////////////////////////////////////
//      Print a Matrix            //
////////////////////////////////////

void print_mat(double * MAT,int n1,int n2, char *name){
      
     int i,j,k=0;
     printf("\n %s= [\n",name);
     for(i=0;i<n1;i++)
     {
//                       printf("\t\t");
                      for(j=0;j<n2;j++)
                      {
                                       printf("%1.3f ",MAT[k]);
					k++;
                      }
                      printf("\n");
      }  
      printf("\t     ];\n"); 
}


////////////////////////////////////
//      Print a Vector            //
///////////////////////////////////
void print_vect(double  *VECT, int n,char *name){
     int i;
     printf("\n %s= [ ",name);
     for(i=0;i<n;i++)
                     printf("%2.5g ",VECT[i]);
     printf(" ]\n");
}


////////////////////////////////////
//      Print a integer Vector    //
///////////////////////////////////
void print_vect_real(int  *VECT, int n,char *name){
     int i;
     printf("\n %s= [ ",name);
     for(i=0;i<n;i++)
                     printf("%d ",VECT[i]);
     printf(" ]\n");
}


////////////////////////////////////
//      Print a ComplexMatrix       //
////////////////////////////////////
void print_mat_C(complex double  *MAT, int n1, int n2, char *name)
{
      
     int i,j,k=0;
     printf("\n %s= [\n",name);
     for(i=0;i<n1;i++)
     {
                      printf("\t");
                      for(j=0;j<n2;j++)
                      {
                            printf("%1.6f+%%i*%1.6f ",crealf(MAT[k]),cimagf(MAT[k]));
														k++;
                      }
                      printf("\n");
      }  
      printf("\t ];\n"); 
}


