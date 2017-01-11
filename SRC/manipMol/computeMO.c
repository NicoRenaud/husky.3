#include "../defMacro.h"
#include "../huckel/header.h"
#include "../algebra_real/algebra.h"
#include "../algebra_real/spec.h"
#include "../huckel/huckel.h"
#include "../manag_output/print_screen.h"


///////////////////////////////////////////////////////////////////////////////////////////////////
// Get the info (valence electron and slater exponents) for specific
// atoms
///////////////////////////////////////////////////////////////////////////////////////////////////

int get_atom_info(double *slater_coeff, double *hard, int index_atom, char *HUCKEL_PARAM)
{
	FILE *Fparameter;
	int count = 0;
  int j;
	int check_occ;
	int index_lumo;
  double coeff1[4];
	double coeff2[4];

	//open the file containing the parameter
	Fparameter=fopen(HUCKEL_PARAM,"r");
	
	// stop reading
	int read=1;
	char atname[5];
	
	// define the variable to read
  int valence_electron;
	
	while(read)
	{
	
	  count++;
		//printf("%d == %d read=%d\n",count,index_atom,read);
		
		// read the atom type and number of valence electron
		fscanf(Fparameter,"%s",atname);
		fscanf(Fparameter,"%d",&valence_electron);
	
		// read the occupation
		index_lumo = -1;
		for(j=0;j<4;j++)
		{
			fscanf(Fparameter,"%d ",&check_occ);
			if(check_occ==0 && index_lumo == -1)
				index_lumo = j;
		}
		
		// read the slater coefficients
		for(j=0;j<4;j++)
		{
			fscanf(Fparameter,"%*f %lf %lf %*f %*f %*f",&coeff1[j],&coeff2[j]);
	  }

		// the coeff of the 2 and p orbitals are on coeff1
		// the coeff of the d and f orbitals are on coeff2
		slater_coeff[0] = coeff1[0];
		slater_coeff[1] = coeff1[1];
		slater_coeff[2] = coeff1[2];
		slater_coeff[3] = coeff1[3];
		
		// check of we still continue to read
		if(count==index_atom)
		{
		 //printf("%d/%d %d %f %f %f %f\n",count,index_atom,valence_electron,slater_coeff[0],slater_coeff[1],slater_coeff[2],slater_coeff[3]);
	   read = 0;
		}
		 
    if(count > 256)
	   read = 0;
		 
	}
	fclose(Fparameter);
	//printf("hard = %f\n",energy[index_lumo-1]-energy[index_lumo]);
	return valence_electron;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Get the atomic info EA IP of given atoms
///////////////////////////////////////////////////////////////////////////////////////////////////

double get_hardness(int index_atom, char *ATOMIC_PARAM)
{
	FILE *Fparameter;
	int count = 0;
	double EA,IP;
	
	//open the file containing the parameter
	Fparameter=fopen(ATOMIC_PARAM,"r");
	
	// stop reading
	int read=1;
	char atname[5];
	
	while(read)
	{
	
	  count++;
		//printf("%d == %d read=%d\n",count,index_atom,read);
		
		// read the atom type and number of valence electron
		fscanf(Fparameter,"%s",atname);
		fscanf(Fparameter,"%lf %lf",&EA,&IP);
		
		// check of we still continue to read
		if(count==index_atom)
	   read = 0;
		
		// if we exeed the number of possible atoms
    if(count > 256)
	   read = 0;
		 
	}
	fclose(Fparameter);
	return(IP-EA);
}




///////////////////////////////////////////////////////////////////////////////////////////////////
//			WRITE THE MO FILE MOPAC STYLE
///////////////////////////////////////////////////////////////////////////////////////////////////

double computeMO(double *H, double *S, atom *molecule, int nb_atom, int nb_orb, int *INDEX_HOMO, char *file_name, char *PARAM)
{

	// varialbles
  int i,j;
	int iprint;
	int iorb,jcomp;
	int count = 0;
	
	// size of the system
  int nb_elec = 0,nb_elec_tot=0,count_elec;
	
	// have to check that
  double *VP_TR;
	
	// slater coefficients
	double *slater_coeff;
	
	// occupations
	int *OCC;
	
	// print info
  int indexHomo;
	int nb_print_orb = 40;
	int half_nb_print_orb = (int) nb_print_orb/2.;
	int *index_print_orb;
	
	// atomic informatio
	int atom_type;
	double hard;
	
	// gap
	double Egap;
	double Ehomo;

	// final file
  FILE *f;
	
	// diagonalization of H
	double *VECT_PRP, *VAL_PRP;
	double *Us, *Vs, *sqrtVs;
	double *Hortho;
	
  // allocate memory
	VECT_PRP = calloc(nb_orb*nb_orb,sizeof(double));
  VAL_PRP = calloc(nb_orb,sizeof(double));
	Us = calloc(nb_orb*nb_orb,sizeof(double));
  Vs = calloc(nb_orb,sizeof(double));
	sqrtVs = calloc(nb_orb*nb_orb,sizeof(double));
	Hortho = calloc(nb_orb*nb_orb,sizeof(double));
  VP_TR = calloc(nb_orb*nb_orb,sizeof(double));
  slater_coeff = calloc(4,sizeof(double));
	OCC = calloc(nb_orb,sizeof(int));
	index_print_orb = calloc(nb_print_orb,sizeof(int));
	
	
  // diagonalize the hamiltonan
	spec_shmidt_zinger(VAL_PRP, VECT_PRP, Hortho, Vs, Us, sqrtVs, H, S, nb_orb);
	
  // read the atomic parameter
  read_atomic_parameters(PARAM);
	
  // debug print
  if (0)
  {
    print_vect(VAL_PRP,nb_orb,"VP");
    print_mat(S,nb_orb,nb_orb,"S");
  }
  
  
  ///////////////////////////////////////////////////////////////////////////////
  //		WRITE THE FILE MOPAC STYLE for JMOL
  ///////////////////////////////////////////////////////////////////////////////
	//printf("\t -- print the MO and the charge densities\n");
	
	
  // open the file
  f = fopen(file_name,"w");
  if(!f) 
  {
      printf("couldn't read %s\n",file_name);
      exit(1);
  }
	
	
  // detect which format
	fprintf(f,"        %d MOPAC-Graphical data Version 2012.13.084W\n",nb_atom);

	
	// print the atoms
	for(i=0;i<nb_atom;i++)
	{
	 //atom_type = findIndex_reduced(molecule[i].atomTypeChar, PARAM);
	 atom_type = findIndexMopac(molecule[i].atomTypeChar);
	 fprintf(f,"%4d %*f%*f%*f  0.0000\n",atom_type,12,molecule[i].x,12,molecule[i].y,12,molecule[i].z);
	}

	 
	 // print the slater exponents
	 for(i=0;i<nb_atom;i++)
	 {
			nb_elec = get_atom_info(slater_coeff, &hard, molecule[i].atomtype,PARAM);
	    nb_elec_tot += nb_elec;
		  fprintf(f,"  %1.7f  %1.7f  %1.7f\n",slater_coeff[0],slater_coeff[1],slater_coeff[2]);
	 }
	 
	 printf("\t -- nb electron counted: %d\n", nb_elec_tot);
	
	 // determine the occupation of each orbital
	 count_elec = nb_elec_tot;
	 for(i=0;i<nb_orb;i++)
	 {
			if(count_elec >= 2)
			{
				OCC[i] = 2;
				count_elec-=2;
			}
			else if(count_elec == 0)
				OCC[i] = 0;
			else
			{
				OCC[i]=1;
				count_elec=0;
			}
			//printf("orbital %d: occupation %d (remaining electrons %d)\n",i+1,OCC[i],count_elec);
	}

	// determine the index of the orbitals to plot
	indexHomo = (int) nb_elec_tot/2.-1;
	if(nb_print_orb<nb_orb)
	{
		for(i=0;i<nb_print_orb;i++)
			index_print_orb[i] = indexHomo-half_nb_print_orb + i;
	}
	else
	{
			nb_print_orb = nb_orb;
			for(i=0;i<nb_orb;i++)
				index_print_orb[i] = i;
	}
	

	// print the orbitals
  for(iprint=0;iprint<nb_print_orb;iprint++)
	{

	  iorb = index_print_orb[iprint];
		fprintf(f," ORBITAL %3d  A  %2.5g\n",OCC[iorb],VAL_PRP[iorb]);
		
		for(jcomp=0;jcomp<nb_orb;jcomp++)
		{
				fprintf(f,"% 1.8E",VECT_PRP[iorb*nb_orb+jcomp]);
				count += 1;
				if(count==5)
				{
				   fprintf(f,"\n");
					 count=0;
				}
		}
		if(count>0)
		{
			fprintf(f,"\n");
			count=0;
		}
	}

	// export the index of the HOMO
	*INDEX_HOMO = indexHomo;


  // print the inverse of S^(1/2)
	count=0;
	fprintf(f,"INVERSE_MATRIX[%dx%d]=\n",nb_orb,nb_orb);
	for(i=0;i<nb_orb;i++)
	{
		for(j=0;j<=i;j++)
		{
			//i = index_print_orb[iprint];
			//j = index_print_orb[jprint];
			fprintf(f,"% 1.8E",sqrtVs[i*nb_orb+j]);
			count += 1;
			if(count==5)
			{
				 fprintf(f,"\n");
				 count=0;
			}
		}
		if(count>0)
		{
			fprintf(f,"\n");
			count=0;
		}
	}
 fprintf(f," Keywords: SYMMETRY GRAPHF");

  
 // close file
 fclose(f);
 
 // print the value of the gap
 Egap = VAL_PRP[indexHomo+1]-VAL_PRP[indexHomo];
 Ehomo = VAL_PRP[indexHomo];
 printf("\t -- Energy Gap %1.3f eV \n",Egap);
 printf("\t        Energy HOMO %1.3f eV \n",VAL_PRP[indexHomo]);
 printf("\t        Energy LUMO %1.3f eV \n",VAL_PRP[indexHomo+1]);
 
 // free memory
 free(VP_TR);
 free(OCC);
 free(slater_coeff);
 free(index_print_orb);
 free(VECT_PRP);
 free(VAL_PRP);
 free(Us);
 free(Vs);
 free(sqrtVs);
 free(Hortho);
 

 
 // return the energy of the HOMO
 return(Ehomo);
}



////////////////////////////////
// compute the occupation
// the atomic charges
// the atomic hardness
////////////////////////////////

int compute_occupation(int *OCC, int *Za, double *HARD, atom *molecule, int nb_orb,int nb_atom,char *HKL_PARAM, char *ATOMIC_PARAM)
{

	int i;
	int count_elec;
	int nb_elec,nb_elec_tot=0;
	double *slater_coeff;
	double hard;
	slater_coeff = calloc(4,sizeof(double));
	
	 // print the slater exponents
	 for(i=0;i<nb_atom;i++)
	 {
			// get the slater coeff and the atomic charges
			nb_elec = get_atom_info(slater_coeff, &hard, molecule[i].atomtype,HKL_PARAM);
			Za[i] = nb_elec;
			
			// get the atomic hardness
			hard = get_hardness(molecule[i].atomtype,ATOMIC_PARAM);
			HARD[i] = hard;
			
			// total number of electrons
	    nb_elec_tot += nb_elec;
		}
	 
	 printf("\t -- nb electron counted: %d\n", nb_elec_tot);
	
	 // determine the occupation of each molecular orbital
	 count_elec = nb_elec_tot;
	 for(i=0;i<nb_orb;i++)
	 {
			if(count_elec >= 2)
			{
				OCC[i] = 2;
				count_elec-=2;
			}
			else if(count_elec == 0)
				OCC[i] = 0;
			else
			{
				OCC[i]=1;
				count_elec=0;
			}
	}
	
	free(slater_coeff);
	return(nb_elec_tot);
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//							write a xyz file of the
//										system
///////////////////////////////////////////////////////////////////////////////////////////////////
void write_system(atom *SYSTEM,int nb_atom_tot,char *file_name)
{
  int  i;
  FILE *f;

  f = fopen(file_name,"w");
  if (!f) 
  {
      printf("can't create %s \n",file_name);
      exit(1);
  } 

  fprintf(f,"   %d\n\n",nb_atom_tot);
  for(i=0;i<nb_atom_tot;i++)
    fprintf(f,"%s  %lf %lf %lf\n",SYSTEM[i].atomTypeChar,SYSTEM[i].x,SYSTEM[i].y,SYSTEM[i].z);
  
  fclose(f);
}







