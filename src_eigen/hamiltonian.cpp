#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <complex>
#include <time.h>
#include <fstream>
#include "Const.h"
#include "Functions.h"


using namespace std;
using namespace Eigen;

static double xij, yij,rij,csh,tnh;

void computeMomentum(const MatrixXd & A_xy, const VectorXd & A_G, double & momentumxValue, double & momentumyValue){
	momentumxValue = 0.0;
	momentumyValue = 0.0;
	//compute the linear momentum
	for(int i =0 ; i < N; i++){
		momentumxValue += A_G(i) * (A_xy(i,0) - 0.5*Lx);
		momentumyValue += A_G(i) * (A_xy(i,1) - 0.5*Ly);
	}

	return;
}

void computeHamiltonian(const MatrixXd & A_xy, const  VectorXd & A_G, double & hamiltonianValue){

hamiltonianValue = 0.0;

xij = 0.0;
yij = 0.0;
rij = 0.0;

//compute the hamiltonian based on the boundary conditions
if(FLAG_BOUNDARY == "infinite"){
	
	for(int j =0 ; j < N; j++){
		for(int i = j+1 ; i < N; i++){
			xij = A_xy(i,0) - A_xy(j,0);
			yij = A_xy(i,1) - A_xy(j,1);
			rij = pow( xij*xij + yij*yij , 0.5);
			hamiltonianValue -= (0.5/pi) * A_G(i) * A_G(j) * log(rij);

		}
	}
}
else if(FLAG_BOUNDARY == "periodic-x"){
	for(int j =0 ; j < N; j++){
		for(int i = j+1 ; i < N; i++){
			xij = A_xy(i,0) - A_xy(j,0);
			yij = A_xy(i,1) - A_xy(j,1);
				
			hamiltonianValue -= (0.5/pi) * A_G(i) * A_G(j) * log(abs(sin((pi/Lx)*complex<double>(xij,yij))));	
		}
		
	}
	if(FLAG_BACKGROUND_FLOW == "sinh"){
		for(int j =0 ; j < N; j++){
			hamiltonianValue += A_G(j)*backgroundFlowStrength*cosh(A_xy(j,1)-0.5*Ly);
		}
	}
	else if(FLAG_BACKGROUND_FLOW == "shear"){
		for(int j =0 ; j < N; j++){
			hamiltonianValue += A_G(j)*backgroundFlowStrength*0.5*pow(A_xy(j,1)-0.5*Ly,2.0);
		}
	}


}
else if(FLAG_BOUNDARY == "periodic"){

	//#pragma omp parallel for schedule(dynamic) private(i,j,xij,yij,m,Ham_m,dHam_m) reduction(+:Ham,Momx,Momy)
	

	for(int j = 0 ; j < N; j++){
		for(int i = j+1 ; i < N; i++){
		
			xij = A_xy(i,0) - A_xy(j,0);
			yij = A_xy(i,1) - A_xy(j,1);

			hamiltonianValue -= (0.25/pi) * A_G(i) * A_G(j) * log(cosh(xij) - cos(yij));

			hamiltonianValue += (0.125 / pow(pi,2.0)) * A_G(i) * A_G(j) * pow(xij,2.0);

			for(int nImage = 1; nImage < imageTotal; nImage++){
				csh = cosh(2.0*pi*nImage);
				tnh = tanh(2.0*pi*nImage);
				
				hamiltonianValue -= (0.25/pi) * A_G(i) * A_G(j) * log( 1.0 + pow(sinh(xij),2.0)*(1.0 - pow(tnh,2.0)) + pow(cos(yij)/csh,2.0) - (2.0*cosh(xij)*cos(yij)/csh)) ;   	
   			}
				
		}
	}
	
}
else{
	cout << "error in computeHamiltonian" << endl;
	exit(1);
}

return;
}


