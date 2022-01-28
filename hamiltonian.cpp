#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include <time.h>
#include <fstream>
#include "Const.h"
#include "Functions.h"
#include <iomanip>

using namespace std;
using namespace arma;

double Ham =0.0;
double dHam_m=0.0;
double Ham_m=0.0;
double Momx, Momy,xij, yij,rij;
int m,i,j;

void hamiltonian(mat A_xy, rowvec A_G, double & ham_out,double & momx_out,double & momy_out){

Ham = 0.0;
dHam_m=0.0;
Ham_m = 0.0;
Momx= 0.0;
Momy=0.0;
xij =0.0;
yij = 0.0;
rij= 0.0;

if(FLAG_PERIODIC_BC == true){

	//#pragma omp parallel for schedule(dynamic) private(i,j,xij,yij,m,Ham_m,dHam_m) reduction(+:Ham,Momx,Momy)
	for(i =0 ; i < N; i++){

		Momx += A_G(i)*(A_xy(i,0)-pi);
		Momy += A_G(i)*(A_xy(i,1)-pi);

		for(j =0 ; j < N; j++){

			if(j != i){
				xij = A_xy(i,0) - A_xy(j,0);
				yij = A_xy(i,1) - A_xy(j,1);

				// compute m-sum for hamiltonian
				Ham_m=0.0;
				dHam_m=1.e0;
				m=1;
				// m =0
				Ham_m = log(cosh(xij) - cos(yij));  // m=0
			for(m=1; m < M_max; m++){
			//	while(dHam_m > ham_eps){

					dHam_m = log( 1. + (sinh(xij)*sinh(xij)*(1. - tanh(2.0*pi*m)*tanh(2.0*pi*m)) ) + (cos(yij)*cos(yij)/(cosh(2.0*pi*m)*cosh(2.0*pi*m))) - (2.0*cosh(xij)*cos(yij)/cosh(2.0*pi*m))) ;

   					Ham_m += dHam_m;
//						cout << i << " " << j << " " << dHam_m << " " << Ham_m << endl;
   				/*	m++;
   					if(m > 100){
   						cout << "Hamiltonian not converged" << endl;
   						cout << "Aborting summation..." << endl;
   						dHam_m = 0.0;
   					}
   				*/
   				}

   			}
   			//cout << " Hamiltonian converged in m = " << m << endl;

			Ham_m -= (0.5/pi)*pow(xij,2.0);
			Ham_m *= 0.25*A_G(i)*A_G(j)/pi;
			Ham -= Ham_m;

		}

	}
}
else{

	//#pragma omp parallel for schedule(dynamic) private(i,j,xij,yij,rij) reduction(+:Ham,Momx,Momy)
	for(i =0 ; i < N; i++){

		Momx += A_G(i)*A_xy(i,0);
		Momy += A_G(i)*A_xy(i,1);

		for(j =0 ; j < N; j++){

			if(i != j){
				xij = A_xy(i,0) - A_xy(j,0);
				yij = A_xy(i,1) - A_xy(j,1);
				rij = pow( xij*xij + yij*yij    ,0.5);
				Ham -= (0.25/pi) * A_G(i)*A_G(j)*log(rij);
			}
		}

	}


}

ham_out = Ham;
momx_out = Momx;
momy_out = Momy;

return;
}

void output_hamiltonian(double runtime, int counter, double ham_out, double momx_out, double momy_out){

ostringstream out_Ham;
out_Ham << "./Output/ham." << setw(5) << setfill('0') << counter << ends;                  //creates file name for outputting da    ta at time slice
string filename_Ham = out_Ham.str();
ofstream fout_ham(filename_Ham.c_str());
fout_ham << scientific;
fout_ham.precision(12);

fout_ham << runtime << " " << ham_out << " " << momx_out << " " << momy_out << endl;

fout_ham.close();

//cout << "Hamiltonian = " << Ham << " Momentum x = " << Momx << " Momentum y = " << Momy << endl;

return;

}
