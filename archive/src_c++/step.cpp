
#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include <time.h>
#include <fstream>
#include <omp.h>
#include "Const.h"
#include "Functions.h"
#include <iomanip>

using namespace std;
using namespace arma;

static double xij, yij, r2ij, cir_i, cir_j;
static double errorSolution,dtNew;
static const double safety = 0.9;

static mat dv(N,2,fill::zeros),k1(N,2,fill::zeros),k2(N,2,fill::zeros),k3(N,2,fill::zeros),k4(N,2,fill::zeros),k5(N,2,fill::zeros),k6(N,2,fill::zeros),sol4(N,2,fill::zeros),sol5(N,2,fill::zeros);




mat rungeKutta45(mat xy, rowvec G, double & dt){


k1.zeros();
k2.zeros();
k3.zeros();
k4.zeros();
k5.zeros();
k6.zeros();

sol4.zeros();
sol5.zeros();


//=========== Runge-Kutta Cash-Karp Algorithm ===========
k1 = dt*computeVelocity(xy, G);
k2 = dt * computeVelocity(xy + (1./5.)*k1, G);
k3 = dt * computeVelocity(xy + (3./40.)*k1 + (9./40.)*k2, G);
k4 = dt * computeVelocity(xy + (3./10.)*k1 - (9./10.)*k2 + (6./5.)*k3, G);
k5 = dt * computeVelocity(xy - (11./54.)*k1 + (5./2.)*k2 - (70./27.)*k3 + (35./27.)*k4, G);
k6 = dt * computeVelocity(xy + (1631./55296.)*k1 + (175./512.0)*k2 + (575./13824.)*k3 + (44275./110592.)*k4 + (253./4096.)*k5, G);

sol4 = xy + (37./378.)*k1 + (250./612.)*k3 + (125./594.)*k4 + (512./1771.)*k6;
sol5 = xy + (2825./27648.)*k1 + (18575./48384.)*k3 + (13525./55296.)*k4 + (277./14336.)*k5 + (1./4.)*k6;



//======== Adaptive Time Stepping ===========

errorSolution = abs(sol4 - sol5).max();


if( errorSolution > 1.e12){ //abort if error is extemely large
	cout << "Error in timestepping routine is too large...code aborted" << endl;
	cout << "errorSolution = " << errorSolution << endl;
	exit(1);
}
else if(errorSolution < rkTolerance){ //if error is less than tolerance, increase dt
	dtNew = min( min( safety * dt * pow( abs(rkTolerance / errorSolution ),0.25), 5.0*dt), dtMax);
	dt = dtNew;
	return sol5;
}
else{ // if error is larger than tolerance then, reduce dt
	dtNew = safety * dt * pow(abs( rkTolerance / errorSolution),0.2);

	dt = max(dtNew,0.1*dt);
	return rungeKutta45(xy,G,dt);
}

}


mat computeVelocity(mat A_xy, rowvec A_G){
	dv.zeros();
		//#pragma omp parallel for schedule(dynamic)
	//	for(i = 0; i < N; i++){
	//		dv.row(i) = compute_vel_BC(A_xy, A_G,i);
	//	}

	if(FLAG_BOUNDARY == "infinite"){
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				if(j != i){
					xij = A_xy(i,0) - A_xy(j,0);
					yij = A_xy(i,1) - A_xy(j,1);
					r2ij = pow(xij,2.0) + pow(yij,2.0);
					dv(i,0) -= (0.5/pi) * A_G(j) * yij / r2ij;
					dv(i,1) += (0.5/pi) * A_G(j) * xij / r2ij;
				}
			}
		}
	}
	else if(FLAG_BOUNDARY == "periodic-x"){
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				if(j != i){
					xij = A_xy(i,0) - A_xy(j,0);
					yij = A_xy(i,1) - A_xy(j,1);

					dv(i,0) += (0.5/Lx) * A_G(j) * imag(1.0/tan( (pi/Lx)*complex<double>(xij,yij)));
					dv(i,1) += (0.5/Lx) * A_G(j) * real(1.0/tan( (pi/Lx)*complex<double>(xij,yij)));

				}
			}
		}
	}
	else if(FLAG_BOUNDARY == "periodic"){
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				if(j != i){
					xij = A_xy(i,0) - A_xy(j,0);
					yij = A_xy(i,1) - A_xy(j,1);

					dv(i,0) -= (0.25/pi) * A_G(j) * sin(yij) / (cosh(xij) - cos(yij));
					dv(i,1) += (0.25/pi) * A_G(j) * sin(xij) / (cosh(yij) - cos(xij));

					for(int nImage = 1 ; nImage < imageTotal; nImage++){
						
						dv(i,0) -= (0.25/pi) * A_G(j) * ( (sin(yij) / (cosh(xij - 2.0*pi*nImage) - cos(yij))  )  + (sin(yij) / (cosh(xij + 2.0*pi*nImage) - cos(yij) )  ));

						dv(i,1) +=   (0.25/pi) * A_G(j) * ( (sin(xij) / (cosh(yij - 2.0*pi*nImage) - cos(xij))  )  + (sin(xij) / (cosh(yij + 2.0*pi*nImage) - cos(xij) ) ) );
					}

				}
			}
		}
	}
	else if (FLAG_BOUNDARY == "disc" ){

		for(int i = 0; i < N; i++){
			cir_i = 1.0 - ((pow(A_xy(i,0),2.0)  + pow(A_xy(i,1),2.0) )/pow(R,2.0));
			
			dv(i,0) -= (0.5/pi) * A_G(i) * A_xy(i,1)  / (pow(R,2.0) * cir_i);
			dv(i,1) += (0.5/pi) * A_G(i) * A_xy(i,0)  / (pow(R,2.0) * cir_i);
		}

		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				if(j != i){
					xij = A_xy(i,0) - A_xy(j,0);
					yij = A_xy(i,1) - A_xy(j,1);

					r2ij = pow(xij,2.0) + pow(yij,2.0);
					cir_i = 1.0 - ((pow(A_xy(i,0),2.0)  + pow(A_xy(i,1),2.0) )/pow(R,2.0));
					cir_j = 1.0 - ((pow(A_xy(j,0),2.0)  + pow(A_xy(j,1),2.0) )/pow(R,2.0));

				
					dv(i,0) -= (0.5/pi) * A_G(j) *  cir_j * ( A_xy(i,1) + (pow(R,2.0) * cir_i * yij / r2ij) )  / (r2ij + (pow(R,2.0) * cir_i * cir_j));
				

					dv(i,1) +=  (0.5/pi) * A_G(j) * cir_j * ( A_xy(i,0)  + (pow(R,2.0) * cir_i * xij / r2ij ) ) / (r2ij + (pow(R,2.0) * cir_i * cir_j));
				}	
				
			}
		}



	}

	if(FLAG_BACKGROUND_FLOW == "shear"){
			for(int i = 0; i < N; i++){
				// shear flow v = [A*y,0]
				dv(i,0) +=  backgroundFlowStrength * (A_xy(i,1));
			}
	}
	else if(FLAG_BACKGROUND_FLOW == "periodic-shear"){
		for(int i = 0; i < N; i++){
	    		if(A_xy(i,1) < 0){        
				dv(i,0) +=  backgroundFlowStrength * (A_xy(i,1) + 0.25*Ly);
            }
			else{
				dv(i,0) -=  backgroundFlowStrength * (A_xy(i,1) - 0.25*Ly);
			}
		}
	}
	else if(FLAG_BACKGROUND_FLOW == "sinh"){
			for(int i = 0; i < N; i++){
				// sinh flow v = [A*sinh(y),0]
				dv(i,0) +=  backgroundFlowStrength * sinh(A_xy(i,1));
			}
	}



	return dv;
	}



void invokeBoundaryCondtions(mat & A_xy){

	if(FLAG_BOUNDARY == "periodic-x"){
		for(int i =0 ; i< N; i++){
			if(A_xy(i,0) < -Lx/2.0){
				A_xy(i,0) += Lx;
			}
			else if(A_xy(i,0) >= Lx/2.0){
				A_xy(i,0) -= Lx;
			}
		}
	}
	else if(FLAG_BOUNDARY == "periodic"){
		for(int i =0 ; i< N; i++){
			if(A_xy(i,0) < -Lx/2.0){
				A_xy(i,0) += Lx;
			}
			else if(A_xy(i,0) >= Lx/2.0){
				A_xy(i,0) -= Lx;
			}
		
			if(A_xy(i,1) < -Ly/2.0){
				A_xy(i,1) += Ly;
			}
			else if(A_xy(i,1) >= Ly/2.0){
				A_xy(i,1) -= Ly;
			}
		}
	}

	return;
}
