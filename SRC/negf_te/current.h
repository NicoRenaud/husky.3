////////////////
// current _H
///////////////

#ifndef _current_h
#define _current_h

double fermi(double EF,double E,double T);
double integrate(double *F, double *X, int ind1, int ind2);
void computeCurrent(double *CURR, double *Vbias, double *TE, double *E, int nb_nrj, double TEMP, double EF);
void computedIdV(double *dIdV, double *IV, double *E, int nb_nrj);
#endif