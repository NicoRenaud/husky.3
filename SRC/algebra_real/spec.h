///////////////////////////////////////
//
//    Proto des fonctions algebra.c
//
///////////////////////////////////////

#ifndef _spec_h
#define _spec_h
void spec_pencil_zinger(double *VAL_PRP, double *VECT_PRP,  double *H, double *S, int nb_orb);
void spec_shmidt_zinger(double *VAL_PRP, double *VECT_PRP, double *Hortho, double *Vs, double *Us, double *sqrtVs, double *H, double *S, int nb_orb);
#endif