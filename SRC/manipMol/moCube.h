///////////////////////////////////////////////
// compute MO	
/////////////////////

#ifndef _moCube_
#define _moCube_

void moCube(double *VAL_PRP, double *VECT_PRP, int index_mo, atom *molecule, int nb_atom, int nb_orb, char *file_name, char *PARAM);
int get_atom_info_all(double *slater_coeff, double *qnumb, int index_atom, char *HUCKEL_PARAM);
double factorial(double n);
double slater_radial(int n, int l, double slater_coeff, double DIST);
double slater_harmonic(int l, int m, double TETA, double PHI);
double compute_value_mo(double *VECT_PRP,int index_mo, double X, double Y, double Z, atom* molecule, int natom, int nb_orb, char *PARAM);

#endif