///////////////////////////////////////////////
// compute MO	
/////////////////////

#ifndef _compute_mo
#define _compute_mo

double computeMO(double *H, double *S, atom *molecule, int nb_atom, int nb_orb, int *INDEX_HOMO, char *file_name, char *PARAM);
int get_atom_info(double *slater_coeff, double *hard, int index_atom,char *HUCKEL_PARAM);
void get_hardness(int index_atom, char *ATOMIC_PARAM);
int compute_occupation(int *OCC, int *Za, double *HARD, atom *molecule, int nb_orb,int nb_atom,char *HKL_PARAM, char *ATOMIC_PARAM);
void write_system(atom *SYSTEM,int nb_atom_tot,char *file_name);
#endif