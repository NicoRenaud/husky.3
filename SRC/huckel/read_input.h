////////////////////////////////////
//	read_data_h
////////////////////////////////////


#ifndef read_data_h
#define read_data_h

void read_number_atoms(char *file_name, int *nb_atom_tot, int *nb_atom_mol, int *nb_atom_elec1, int *no_atom_elec2);
void read_general_info(char *file_name, double *emin, double *emax, int *nb_energy, int *export_mat, double *LDOS, char *HKL_PARAM);
void read_xyz(char *file_name, atom* SYS,int nb_atom_tot,int nb_atom_mol,int *index_atom_mol,
									int nb_atom_elec1,int *index_atom_elec_1,
									int nb_atom_elec2,int *index_atom_elec_2,
									char *HKL_PARAM);

#endif