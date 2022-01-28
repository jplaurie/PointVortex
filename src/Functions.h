#include <armadillo>

using namespace arma;

void initialize_vortices(mat &, rowvec &,double &,int &);
mat step_vortices(mat , rowvec);
mat RK45(mat &, rowvec ,mat & , double &);
void invoke_BC(mat &);
void output_vortex(mat,rowvec,double, int);
void output_vel(mat,double, int);
void output_remove(double,int, int &, int&, double &, double &,double &, double &);
void hamiltonian(mat, rowvec, double &, double &, double &);
void output_hamiltonian(double, int, double, double, double);
void remove(mat & ,rowvec &,int &, int &, double &, double &,double &, double & );
void output_structure(mat, rowvec, mat,int,irowvec);