
//#include "./manag_input/read_input_file.h"
//#include "./manag_input/create_system.h"
//#include "./manag_output/export_scilab.h"
//#include "./algebra_cplx/algebra_cplx.h"
//#include "./huckel/read_input_hkl.h"
//#include "./huckel/read_input_hkl_general.h"
//#include "./manipMol/rotate_molecule.h"

#include "./defMacro.h"
//#include "./manag_output/print_file.h"
#include "./manag_output/print_screen.h"
#include "./manag_output/export_file.h"
#include "./manag_output/export_gnuplot.h"

#include "./algebra_real/algebra.h"
#include "./algebra_real/extract_hamiltonian.h"



#include "./huckel/header.h"
#include "./huckel/huckel.h"
#include "./huckel/read_input.h"

#include "./negf_te/negf_te.h"
#include "./negf_te/current.h"


#include "./manipMol/computeMO.h"


///////////////////////////////////////////////////////////
////////					///////////
////////	MAIN FUNCTION			///////////
////////					///////////
///////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{

  printf("\n====  ===============================  ====\n");    
  printf("====                                   ====\n");   
  printf("====  	        Husky.3	               ====\n");
	printf("====    Transport Calculation with     ====\n");
	printf("====         Extended Huckel           ====\n");
	printf("====         N. Renaud 2015            ====\n");
  printf("====                                   ====\n");   
  printf("====  ===============================  ====\n\n");



  ////////////////////////
  //	DECLARATION	//
  ////////////////////////

  // molecular Hamiltonian
  double *H,*S,*Hmol,*Smol;
	double *Vl,*Sl,*Vr,*Sr;
	int index_homo;

	// atomic structure information
	atom *MOL;
	atom *SYS;
	
	// total number of atoms and orbitals
  int nb_atom_tot;
  int nb_orb_tot = 0;

  // molecule information
  int nb_atom_mol = 0;
  int nb_orb_mol = 0;
  int *index_orb_mol=NULL;
  int *index_atom_mol=NULL;

  
  // index of the electrode atoms in the total system
  int *index_atom_elec_1=NULL;
	int *index_atom_elec_2=NULL;
  int nb_atom_elec_1;
  int nb_atom_elec_2;
  
  // index of the electrode orbitals in the system
  int *index_orb_elec_1=NULL, *index_orb_elec_2=NULL;
  int nb_orb_elec_1 = 0;
  int nb_orb_elec_2 = 0;
 
	// arrays for TE calculations
	double *TE;
	double *DOS;
	double *ENERGY;
	double emin,emax;
	int nb_energy;

	// do we export the matrices
	int export_mat;
 
	// path for the huckel parameters
	char HKL_PARAM[250];
 
  // out files
  char  out_file_system[200];
  strcpy(out_file_system,argv[2]);
  strcat(out_file_system,"system.xyz");
  
  char om_file[200];
  strcpy(om_file,argv[2]);
  strcat(om_file,"orbitals.dat");
  
  char file_te[200];
  strcpy(file_te,argv[2]);
  strcat(file_te,"TE.dat");
  
  char file_dos[200];
  strcpy(file_dos,argv[2]);
  strcat(file_dos,"DOS.dat");
  
  char file_IV[200];
  strcpy(file_IV,argv[2]);
  strcat(file_IV,"IV.dat");
  
  char file_dIdV[200];
  strcpy(file_dIdV,argv[2]);
  strcat(file_dIdV,"dIdV.dat");
	
  char hname[200];
  strcpy(hname,argv[2]);
  strcat(hname,"hmol.dat");

  char sname[200];
  strcpy(sname,argv[2]);
  strcat(sname,"smol.dat");

  // local density of states of the electrodes
  double LDOS; // = 0.18;

  ////////////////////////////////////////////////////////
  //  read main input
	//  the input is now read in different batches
	//	to avoid memory leaks
	//
	//	there is probably a better way of doing this
  ///////////////////////////////////////////////////////

	printf(" == Read the input file\n");
	
	// read general information
	read_general_info(argv[1], &emin, &emax,&nb_energy,&export_mat,&LDOS,HKL_PARAM);
	
	// read the number of atoms in the system
	read_number_atoms(argv[1],  &nb_atom_tot, &nb_atom_mol, &nb_atom_elec_1, &nb_atom_elec_2);
	
	// allocate the memory for the atomic information
	SYS = calloc(nb_atom_tot,sizeof(atom));
	index_atom_mol = calloc(nb_atom_mol,sizeof(int));
	index_atom_elec_1 = calloc(nb_atom_elec_1,sizeof(int));
	index_atom_elec_2 = calloc(nb_atom_elec_2,sizeof(int));
	
	
	// read/store atom-type and position
	read_xyz(argv[1],SYS,nb_atom_tot,
									nb_atom_mol,index_atom_mol,
									nb_atom_elec_1,index_atom_elec_1,
									nb_atom_elec_2,index_atom_elec_2,HKL_PARAM);
	
	// write a xyx of the junction
	write_system(SYS, nb_atom_tot,out_file_system);
		
  ///////////////////////////////////////////////////////
	//	determine the number of orbitals in
	//	the different section of the junction
  ///////////////////////////////////////////////////////
	
	// determine the number of orbitals
	nb_orb_mol = compute_nb_orb(SYS, nb_atom_mol,index_atom_mol,HKL_PARAM);
	nb_orb_elec_1 = compute_nb_orb(SYS, nb_atom_elec_1,index_atom_elec_1,HKL_PARAM);
	nb_orb_elec_2 = compute_nb_orb(SYS, nb_atom_elec_2,index_atom_elec_2,HKL_PARAM);
	
	// determine the total number of orbitals
	nb_orb_tot = nb_orb_mol+nb_orb_elec_1+nb_orb_elec_2;
	
	printf("\t -- %d atom detected in the system: \n \t\t %d in the first electrode \n \t\t %d in the molecule \n \t\t %d in the second electrode\n",nb_atom_tot,nb_atom_elec_1,nb_atom_mol,nb_atom_elec_2);
	printf("\t -- %d orbitals detected in the system: \n \t\t %d in the first electrode \n \t\t %d in the molecule \n \t\t %d in the second electrode\n",nb_orb_tot,nb_orb_elec_1,nb_orb_mol,nb_orb_elec_2);
	
	// create the index of the orbitals of the different part
	index_orb_elec_1 = calloc(nb_orb_elec_1,sizeof(int));
	index_orb_elec_2 = calloc(nb_orb_elec_2,sizeof(int));
	index_orb_mol = calloc(nb_orb_mol,sizeof(int));
	create_index_orbitals(index_orb_elec_1,index_orb_mol,index_orb_elec_2,nb_orb_elec_1,nb_orb_mol,nb_orb_elec_2); 
	
	
  ///////////////////////////////////////////////////////
	//	compute the Hamiltonian and overlaps
	//	in the junctions
	//  Hmol,Smol : hamiltonian/overlap of the molecule
	//  Vl,Sl : interaction/overlap left contact
	//  Vr,Sr : interaction/overlap right contact
  ///////////////////////////////////////////////////////
	
	// allocate the memory for and cmpute the total hamiltonian and overlap matrix
	H = calloc(nb_orb_tot*nb_orb_tot,sizeof(double));
	S = calloc(nb_orb_tot*nb_orb_tot,sizeof(double));
	printf(" == Compute the Hamiltonian/Overlap matrices\n");
	compute_hamiltonian(H, S, SYS, nb_atom_tot, nb_orb_tot, HKL_PARAM);
	
	printf(" == Extract the sub matrices\n");
	// extract the Hamiltonian/Overlap of the molecule
	Hmol = calloc(nb_orb_mol*nb_orb_mol,sizeof(double));
	Smol = calloc(nb_orb_mol*nb_orb_mol,sizeof(double));
	extract(Hmol, H, index_orb_mol,nb_orb_mol, index_orb_mol,nb_orb_mol, nb_orb_tot);
	extract(Smol, S, index_orb_mol, nb_orb_mol, index_orb_mol, nb_orb_mol, nb_orb_tot);
	
	
	// extract the Hamiltonian/Overlap of the molecule
	Vl = calloc(nb_orb_mol*nb_orb_elec_1,sizeof(double));
	Sl = calloc(nb_orb_mol*nb_orb_elec_1,sizeof(double));
	extract(Vl, H, index_orb_mol,nb_orb_mol, index_orb_elec_1,nb_orb_elec_1, nb_orb_tot);
	extract(Sl, S, index_orb_mol,nb_orb_mol, index_orb_elec_1,nb_orb_elec_1, nb_orb_tot);
	
	// extract the Hamiltonian/Overlap of the molecule
	Vr = calloc(nb_orb_mol*nb_orb_elec_2,sizeof(double));
	Sr = calloc(nb_orb_mol*nb_orb_elec_2,sizeof(double));
	extract(Vr, H, index_orb_elec_2,nb_orb_elec_2, index_orb_mol,nb_orb_mol, nb_orb_tot);
	extract(Sr, S, index_orb_elec_2,nb_orb_elec_2, index_orb_mol,nb_orb_mol, nb_orb_tot);
	
	
	
	
	
  ///////////////////////////////////////////////////////
	//	compute the transmission of the junction
	//	using (NE)GF
  ///////////////////////////////////////////////////////
	
	// allocate memory
	TE = calloc(nb_energy,sizeof(double));
	ENERGY = calloc(nb_energy,sizeof(double));
	DOS = calloc(nb_energy,sizeof(double));
	
	// define energy range
	linspace(ENERGY, emin,  emax, nb_energy);
	
	// compute TE
	printf(" == Compute the electronic transmission\n");
	NEGF_TE(TE, DOS,  ENERGY, LDOS, Hmol, Smol, Vl, Sl, Vr, Sr, nb_orb_mol, nb_orb_elec_1, nb_orb_elec_2, nb_energy);
	

	// print the TE
	export_gnulot_semilog(file_te, TE, ENERGY, nb_energy);
	export_gnulot(file_dos, DOS, ENERGY, nb_energy);
	
	printf(" == Transmission spectrum written in      : %s\n",file_te);
	printf(" == Molecule density of state written in  : %s\n",file_dos);
	 
	
	/////////////////////////////////////////////////
  // 	   	compute the MO												 //
  /////////////////////////////////////////////////   
	printf(" == Compute the molecular orbitals\n");
	
	// extract the molecule information
	MOL = calloc(nb_atom_mol,sizeof(atom));
	extract_atom_mol(MOL,index_atom_mol,nb_atom_mol,SYS);

	// compute/print the MO
	computeMO(Hmol, Smol, MOL, nb_atom_mol,nb_orb_mol, &index_homo, om_file, HKL_PARAM);
	printf(" == Molecular Orbitals written in         : %s\n",om_file);    

  /////////////////////////////////////////////////
  // export the matrices if necessary 	   	
  /////////////////////////////////////////////////   
  if(export_mat==1)
  {
  	export_mat_tofile(Hmol,nb_orb_mol,hname);
  	export_mat_tofile(Smol,nb_orb_mol,sname);
  	printf(" == Hamiltonian  matrix written in        : %s\n",hname);
  	printf(" == Overlap  matrix written in            : %s\n",sname);
  	
  }
  

  ///////////////////////////////////
  //	free memory
  ///////////////////////////////////
	free(SYS);
	free(MOL);
	free(index_atom_mol);
	free(index_atom_elec_1);
	free(index_atom_elec_2);
	free(index_orb_mol);
	free(index_orb_elec_1);
	free(index_orb_elec_2);
	free(H);
	free(S);
	free(Hmol);
	free(Smol);
	free(Vl);
	free(Sl);
	free(Vr);
	free(Sr);
	free(TE);
	free(ENERGY);
	free(DOS);
	
	
   ////////////////////////////////////////////////
  // 	   	Say Good bye			//
  /////////////////////////////////////////////// 
  printf("\n====  ===============================  ====\n");
  printf("====  Thanks for using Husky.1 	       ====\n");
  printf("====                    see you later  ====\n");
  printf("====  ===============================  ====\n");

  return(0);
}
 
 
