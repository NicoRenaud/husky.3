#include "../defMacro.h"
#include "./header.h"
#include "./huckel.h"



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//							read the general information of the system
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_general_info(char *file_name, double *emin, double *emax, int *nb_energy, int *export_mat, double *LDOS, char *HKL_PARAM)
{
	
  int INT;
	double FLT;
  char str_cmt[250];
  char str[250];
	
  FILE *f;
  
  // defalut values
  *export_mat = 0;
  *LDOS = 0.18;

  
   // check if the file exists
   f = fopen(file_name,"r");
   if (!f) {
      printf("couldn't read input file name %s\n", file_name);
      exit(1);
   }
   
   //////////////////////////////////////////
   // for all lines of the input file
   //////////////////////////////////////////
   while (!feof(f))
   {
     
			///////////////////////
			// read the string
			///////////////////////
			fscanf(f,"%s",str);
		 
			////////////////////////
			// if it's a comment
			////////////////////////
			strncpy(str_cmt,str,2);
			str_cmt[2] = '\0';
			if (!strcmp(str_cmt,"//"))
				fgets(str,100,f);
			
			///////////////////////////
			// if it's not a comment
			///////////////////////////
			else
			{

					////////////////////////
					// EHMO parameters
					///////////////////////    
					if( !strcmp(str,"parameters") || !strcmp(str,"PARAMETERS")   ||
					!strcmp(str,"param")  || !strcmp(str,"PARAM")  ) 
					{
						fscanf(f,"%s",HKL_PARAM);
					}

					////////////////////////
					// min energy
					///////////////////////    
					 else if( !strcmp(str,"min_energy") || !strcmp(str,"MIN_ENERGY")   || 
										!strcmp(str,"energy_min")  ||  !strcmp(str,"ENERGY_MIN")  ||
										!strcmp(str,"emin") || !strcmp(str,"EMIN")	)
					{  
							fscanf(f,"%lf",&FLT);
							*emin=FLT;
					}
					
					////////////////////////
					// max energy
					///////////////////////    
					 else if( !strcmp(str,"max_energy") || !strcmp(str,"MAX_ENERGY")   || 
										!strcmp(str,"energy_max")  ||  !strcmp(str,"ENERGY_MAX")  ||
										!strcmp(str,"emax") || !strcmp(str,"EMAX")	)
					{  
								fscanf(f,"%lf",&FLT);
								*emax=FLT;
					}	
					
					////////////////////////
					// nb energy
					///////////////////////    
					 else if( !strcmp(str,"nb_energy") || !strcmp(str,"NB_ENERGY")   || 
										!strcmp(str,"energy_nb")  ||  !strcmp(str,"ENERGY_NB")  ||
										!strcmp(str,"nb_e") || !strcmp(str,"NB_E")	)
					{  
							fscanf(f,"%d",&INT);
							*nb_energy=INT;
					}

					////////////////////////
					// export matrix
					///////////////////////    
					 else if( !strcmp(str,"export_matrices") || !strcmp(str,"EXPORT_MATRICES")   || 
										!strcmp(str,"exp_mat")  ||  !strcmp(str,"EXP_MAT") )
					{  
							*export_mat=1;
					}

					////////////////////////
					// LDOS
					///////////////////////    
					 else if( !strcmp(str,"LDOS") || !strcmp(str,"ldos") )
					{  
							fscanf(f,"%lf",&FLT);
							*LDOS=FLT;
					}
			}
			
		}
		// close the file
		fclose(f);
	
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//							read the number of atoms in the different parts
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_number_atoms(char *file_name, int *nb_atom_tot, int *nb_atom_mol, int *nb_atom_elec_1, int *nb_atom_elec_2)
{
	
  int INT,INT2;
	
  //int current_cell = 0;
  //int read_cell[2];
	
  //int nb_atom_tot_read=0;
  //int nb_atom_mol_read=0;
	
  //int *nb_atom_cell_read;
  //nb_atom_cell_read = calloc(2,sizeof(int));
  //nb_atom_cell_read[0]  = 0;
  //nb_atom_cell_read[1]  = 0;
	
  char str_cmt[5];
  char str[100];
	
  //int index_check;
  
  FILE *f;
  
  
   // check if the file exists
   f = fopen(file_name,"r");
   if (!f) {
      printf("couldn't read input file name %s\n", file_name);
      exit(1);
   }
   
   //////////////////////////////////////////
   // for all lines of the input file
   //////////////////////////////////////////
   while (!feof(f))
   {
     
			///////////////////////
			// read the string
			///////////////////////
			fscanf(f,"%s",str);
		 
			////////////////////////
			// if it's a comment
			////////////////////////
			strncpy(str_cmt,str,2);
			str_cmt[2] = '\0';
			if (!strcmp(str_cmt,"//"))
				fgets(str,100,f);
			
			///////////////////////////
			// if it's not a comment
			///////////////////////////
			else
			{
						 
					////////////////////////
					// total number of atoms
					///////////////////////    
					if( !strcmp(str,"nb_atom_tot") || !strcmp(str,"NB_ATOM_TOT")   || 
					!strcmp(str,"nb_atom")  || !strcmp(str,"NB_ATOM")  ) 
					{
						fscanf(f,"%d",nb_atom_tot);
					}  
				
				
					////////////////////////
					// total number electrode
					///////////////////////    
					else  if( !strcmp(str,"nb_atom_elec") || !strcmp(str,"NB_ATOM_ELECTRODE")   || 
										!strcmp(str,"nb_atom_elec")  || !strcmp(str,"NB_ATOM_ELEC")  )
					{
							fscanf(f,"%d %d",&INT,&INT2);
							if(INT>2)
								printf("Only two electrodes are possible, check file %s",file_name);

							if(INT == 1)
								*nb_atom_elec_1 = INT2;
							if(INT == 2)
								*nb_atom_elec_2 = INT2;
					 }
					 
					 ////////////////////////
					 // nb atom molecule
					 ///////////////////////
					 else if( !strcmp(str,"nb_atom_molecule") || !strcmp(str,"NB_ATOM_MOLECULE")   || 
										!strcmp(str,"nb_atom_mol")  || !strcmp(str,"NB_ATOM_MOL")  )
					 {
							fscanf(f,"%d",nb_atom_mol);
					 }
				
			}
			
		}
		// close the file
		fclose(f);
}






////////////////////////////////////////////////////////////////////////////////////////////////////
//
//			Read the xyz information of the system
//
////////////////////////////////////////////////////////////////////////////////////////////////////

void read_xyz(char *file_name, atom* SYS,int nb_atom_tot,int nb_atom_mol,int *index_atom_mol,
									int nb_atom_elec1,int *index_atom_elec_1,
									int nb_atom_elec2,int *index_atom_elec_2,
									char *HKL_PARAM)
{
  
  int i;
  //int INT;
  //int INT2;
  float X,Y,Z;
	
	int current_elec;
	int read_elec[2];
	int read_mol;
	
  int nb_atom_tot_read=0;
	int nb_atom_mol_read=0;
	
  int *nb_atom_elec_read;
  nb_atom_elec_read = calloc(2,sizeof(int));
	nb_atom_elec_read[0]=0;
	nb_atom_elec_read[1]=0;
	
	//int current_cell = 0;
  //int read_cell[2];
  //int *nb_atom_cell_read;
  //nb_atom_cell_read = calloc(2,sizeof(int));
  //nb_atom_cell_read[0]  = 0;
  //nb_atom_cell_read[1]  = 0;
	
  
  char str_cmt[5];
  char str[100];
	
	
  int index_check;
  
  FILE *f;
  
  
  // check if the file exists
   f = fopen(file_name,"r");
   if (!f) {
      printf("couldn't read input file name %s\n", file_name);
      exit(1);
   }
   
   //////////////////////////////////////////
   // for all lines of the input file
   //////////////////////////////////////////
   while (!feof(f))
   {
     
			///////////////////////
			// read the string
			///////////////////////
			fscanf(f,"%s",str);
		 
			////////////////////////
			// if it's a comment
			////////////////////////
			strncpy(str_cmt,str,2);
			str_cmt[2] = '\0';
			if (!strcmp(str_cmt,"//"))
				fgets(str,100,f);
			
			///////////////////////////
			// if it's not a comment
			///////////////////////////
			else
			{
				
					////////////////////////////////////////////////////////////////////////////////
					//				WHICH POSITIONS ARE WE READING
					///////////////////////////////////////////////////////////////////////////////
		  
					////////////////////////
					// electrode
					///////////////////////    
					 if( !strcmp(str,"electrode") || !strcmp(str,"ELECTRODE")   ||
										!strcmp(str,"elec")  || !strcmp(str,"ELEC")  )
					 {
							// read the current elec
							fscanf(f,"%d",&current_elec);

							// reinit the value of read_X
							for(i=0;i<2;i++)
								read_elec[i] = 0;
							read_mol = 0;
							
							// init value of read_elec(current)
							read_elec[current_elec-1]  = 1;
					 }
				
					////////////////////////
					// molecule
					///////////////////////    
					 else if( !strcmp(str,"molecule") || !strcmp(str,"MOLECULE")   || 
										!strcmp(str,"mol")  || !strcmp(str,"MOL")  )
					{
							for(i=0;i<2;i++)
									read_elec[i] = 0;
							read_mol = 1;
					}
			
				 ////////////////////////////////////////////////////////////////////////////////
				 //			READING POSITIONS
				 ///////////////////////////////////////////////////////////////////////////////
										
					else
					{
						 
						 index_check = findIndex_reduced(str,HKL_PARAM);
						 //printf("str : %s index : %d\n",str,index_check);
						 if (index_check > 0)
						 {
			
								// read the data
								fscanf(f,"%f%f%f",&X,&Y,&Z);
		
			
								// store the data
								strcpy(SYS[nb_atom_tot_read].atomTypeChar,str);
								SYS[nb_atom_tot_read].atomtype = index_check;
								SYS[nb_atom_tot_read].x=X;
								SYS[nb_atom_tot_read].y=Y;
								SYS[nb_atom_tot_read].z=Z;
								 
								// update the nuber of atoms and their indexes of the molecule
								if(read_mol)
								{
										// store the index
										index_atom_mol[nb_atom_mol_read] = nb_atom_tot_read;
									
										// one more atom in the molecule
										nb_atom_mol_read ++;
								}
								
								
								// update the #of atoms and their indexes of the first electrode
								if(read_elec[0])
								{

										// store the index
										index_atom_elec_1[ nb_atom_elec_read[0] ] = nb_atom_tot_read;
									
										// one more atom in the first elec
										nb_atom_elec_read[0] ++;
								}
								
								// update the #of atoms and their indexes of the second electrode
								if(read_elec[1])
								{
								
										// store the index
										index_atom_elec_2[ nb_atom_elec_read[1] ] = nb_atom_tot_read;
									
										// one more atom in the secon elec
										nb_atom_elec_read[1] ++;
								}
									
								// one more atom read!
								nb_atom_tot_read ++;
		
						}
					}
				}
			}



		// check if everything is correct
		if(nb_atom_tot_read != nb_atom_tot)
			printf("\nWarning: input claims there is %d atom in the system but %d have been found \n Please check %s\n",nb_atom_tot,nb_atom_tot_read,file_name);
			 
		if(nb_atom_mol_read != nb_atom_mol)
			printf("\nWarning: input claims there is %d atom in the molecule but %d have been found \n Please check %s\n",nb_atom_mol,nb_atom_mol_read,file_name);

		if( nb_atom_elec_read[0]!= nb_atom_elec1)
			printf("\nWarning: input claims there is %d electrode 1 in the molecule but %d have been found \n Please check %s\n",nb_atom_elec1,nb_atom_elec_read[0],file_name);
		
		 if( nb_atom_elec_read[1]!= nb_atom_elec2)
			printf("\nWarning: input claims there is %d electrode 2 in the molecule but %d have been found \n Please check %s\n",nb_atom_elec2,nb_atom_elec_read[1],file_name);
	
	
		// free memory
		free(nb_atom_elec_read);
		
		// close the file
		fclose(f);

}












