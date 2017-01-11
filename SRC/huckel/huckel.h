
#include "./header.h"

#ifndef _comp_HKL_H_h
#define _comp_HKL_H_h

// fin the index of a given atom type
int findIndex_reduced(char *name, char *PARAM);

// find the index for Mopac printing
int findIndexMopac(char *name);

// extract a subatrix
void extract_atom_mol(atom *mol_only,  int *ind_atom_mol, int no_atom_mol, atom *MOL);

// function to computes overlaps
void overlap(int atom_row,int atom_col,double delx,double dely,double delz,double S[16][16],double H[16][16],double KEHT);
void mov(double *sigma,double *pi,double *delta,double *phi,int atom_col,int atom_row,double rr,int n1,int n2,int l1,int l2);
void abfns(double *a,double *b,double sk1,double sk2,double rr,int maxcal);
double lovlap(double *a,double *b,double sk1,double sk2,double r,int l1,int l2,int m1,int n1,int n2);

// compute the hamiltonian of a given system
void compute_hamiltonian(double *hmat, double *smat, atom *molecule, int no_atoms, int nb_orb_precomp, char *HKL_PARAM);

// reand the huckel aprameters
void read_atomic_parameters(char * HUCKEL_PARAM);

// compte the number of orbitals in a molecule
int compute_nb_orb(atom *SYS,int nb_atom, int *index_atom, char *HKL_PARAM);
void create_index_orbitals(int *index_orb_elec_1,int *index_orb_mol,int *index_orb_elec_2,int nb_orb_elec_1,int nb_orb_mol,int nb_orb_elec_2);

#endif