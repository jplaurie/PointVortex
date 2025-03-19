#include <armadillo>

using namespace arma;

void initializeVortices(mat &, rowvec &,double &,int &);
mat computeVelocity(mat , rowvec);
mat rungeKutta45(mat , rowvec , double &);
void invokeBoundaryCondtions(mat &);

void printVortex(mat,rowvec,double, int);
void printParameters();
void printDiagnostics(double , double , double, double , double , int );
//void output_vel(mat,double, int);
void printRemoval(double,int, int &, int&, double &, double &,double &);
void computeMomentum(mat, rowvec, double & , double &, double & );
void computeHamiltonian(mat, rowvec, double &);
void removeVortices(mat & ,rowvec &,int &, int &, double &, double &,double & );
void mergeVortices(mat & ,rowvec & );
//void output_structure(mat, rowvec, mat,int,irowvec);
