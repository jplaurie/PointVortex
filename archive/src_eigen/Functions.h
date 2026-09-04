#include <eigen3/Eigen/Dense>

using namespace Eigen;

void initializeVortices(MatrixXd &, VectorXd &,double &,int &);
MatrixXd computeVelocity(const MatrixXd &, const VectorXd & );
MatrixXd rungeKutta45(const MatrixXd & , const VectorXd & , double &);
void invokeBoundaryCondtions(MatrixXd &);

void printVortex(const MatrixXd &, const VectorXd & ,double, int);
void printParameters();
void printDiagnostics(double , double , double , double , int );
//void output_vel(mat,double, int);
void printRemoval(double,int, int &, int&, double &, double &,double &);
void computeMomentum(const MatrixXd &,const VectorXd &, double & , double & );
void computeHamiltonian(const MatrixXd &,const  VectorXd &, double &);
void removeVortices(MatrixXd & ,VectorXd  &,int &, int &, double &, double &,double & );
void mergeVortices(MatrixXd & ,VectorXd  & );
//void output_structure(mat, rowvec, mat,int,irowvec);
