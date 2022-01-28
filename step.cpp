
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

mat compute_vel_BC(mat,rowvec, int );

static mat dv(N,2,fill::zeros),periodic_box(N,2,fill::zeros),k1(N,2,fill::zeros),k2(N,2,fill::zeros),k3(N,2,fill::zeros),k4(N,2,fill::zeros),k5(N,2,fill::zeros),k6(N,2,fill::zeros),sol4(N,2,fill::zeros),sol5(N,2,fill::zeros),xy_temp(N,2,fill::zeros);

static double error,dt_new,Safety;



mat RK45(mat & xy, rowvec G, mat & vortex_vel, double &dt){


k1.zeros();
k2.zeros();
k3.zeros();
k4.zeros();
k5.zeros();
k6.zeros();

sol4.zeros();
sol5.zeros();

xy_temp.zeros();

Safety=0.9;

xy_temp = xy;

//=========== RK45-Cash-Karp ===================================
k1 = dt*step_vortices(xy_temp, G);

xy_temp = xy + (1./5.)*k1;
k2 = dt*step_vortices(xy_temp, G);
xy_temp = xy + (3./40.)*k1 + (9./40.)*k2;
k3 = dt*step_vortices(xy_temp, G);
xy_temp = xy + (3./10.)*k1 - (9./10.)*k2 + (6./5.)*k3;
k4 = dt*step_vortices(xy_temp, G);
xy_temp = xy - (11./54.)*k1 + (5./2.)*k2 - (70./27.)*k3 + (35./27.)*k4;
k5 = dt*step_vortices(xy_temp, G);
xy_temp = xy + (1631./55296.)*k1 + (175./512.0)*k2 + (575./13824.)*k3 + (44275./110592.)*k4 + (253./4096.)*k5;
k6 = dt*step_vortices(xy_temp, G);

sol4 = xy + (37./378.)*k1 + (250./612.)*k3 + (125./594.)*k4 + (512./1771.)*k6;
sol5 = xy + (2825./27648.)*k1 + (18575./48384.)*k3 + (13525./55296.)*k4 + (277./14336.)*k5 + (1./4.)*k6;
vortex_vel = (1./dt)*((2825./27648.)*k1 + (18575./48384.)*k3 + (13525./55296.)*k4 + (277./14336.)*k5 + (1./4.)*k6);
//=============================================================

error = abs(sol4 - sol5).max();

if( error > 1.e12){
	cout << "Error in timestepping routine is too large....aborting...." << endl;
	cout << "error = " << error << endl;
	exit(1);
}
else if(error < rktol){
	xy = sol5;
	dt_new = min(Safety*dt*pow(abs(rktol/error),0.25), 5.0*dt);
	dt = dt_new;
	return xy;
}
else{

	dt_new = Safety*dt*pow(abs(rktol/error),0.2);
	dt = max(dt_new,0.1*dt);
	xy = RK45(xy,G,vortex_vel,dt);
	return xy;
}

}


mat step_vortices(mat A_xy, rowvec A_G){
		int i;
		dv.zeros();
		#pragma omp parallel for schedule(dynamic)
		for(i = 0; i < N; i++){
			dv.row(i) = compute_vel_BC(A_xy, A_G,i);
		}
	/*	if(FLAG_BACKGROUND_SHEAR == true){
			for(i = 0; i < N; i++){
				dv(i,0) += shearStrength*(A_xy(i,1)-0.5*Ly);
			}
		}
*/
	return dv;
	}




mat compute_vel_BC(mat A_xy,rowvec A_G,int i){

	mat dv_temp(1,2,fill::zeros),dv_m(1,2,fill::zeros);
	int m;
	double xij=0.0;
	double yij=0.0;
	double r2ij = 0.0;
	//double dvmax = 1.0;
	if(FLAG_PERIODIC_BC == true){

		for(int j = 0; j < N; j++){

			if(j != i){
				xij = A_xy(i,0) - A_xy(j,0);
				yij = A_xy(i,1) - A_xy(j,1);

				dv_m(0,0) = -(0.25/pi) * A_G(j) * sin(yij) / (cosh(xij) - cos(yij));
				dv_m(0,1) = (0.25/pi) * A_G(j) * sin(xij) / (cosh(yij) - cos(xij));

				dv_temp += dv_m;

				m=1;
		//		dvmax= 1.0;
				for(m =1 ; m < M_max; m++){
				//while( dvmax > vel_eps){
					dv_m(0,0) = - (sin(yij) / (cosh(xij - 2.0*pi*m) - cos(yij))  )  - (sin(yij) / (cosh(xij + 2.0*pi*m) - cos(yij) )  );
					dv_m(0,1) =   (sin(xij) / (cosh(yij - 2.0*pi*m) - cos(xij))  )  + (sin(xij) / (cosh(yij + 2.0*pi*m) - cos(xij)  ) );

					dv_temp += (0.25/pi) * A_G(j) * dv_m;

		//			dvmax =  abs(dv_m).max();
				/*	m++;
					if(m > 100){
						cout << "Velocity not converged" << endl;
						cout << "Aborting summation..." << endl;
						dvmax = 0.0;
					}
					*/
				}
				//	cout << "m = " << m << " i = " << i << " j = " << j << endl;
			}
		}
		return dv_temp;
	}
	else{

		for(int j = 0; j < N; j++){

			if(j != i){
				xij = A_xy(i,0) - A_xy(j,0);
				yij = A_xy(i,1) - A_xy(j,1);
				r2ij = pow(xij,2.0) + pow(yij,2.0);
				dv_temp(0,0) -= (0.5/pi) * A_G(j) * yij  / r2ij;
				dv_temp(0,1) += (0.5/pi) * A_G(j) * xij / r2ij;
			}

		}


		return dv_temp;
	}
}


void invoke_BC(mat & A_xy){


		for(int i =0 ; i< N; i++){

			if(A_xy(i,0) < 0){
				A_xy(i,0) += Lx;
			}
			else if(A_xy(i,0) >= Lx){
				A_xy(i,0) -= Lx;
			}

			if(A_xy(i,1) < 0){
				A_xy(i,1) += Ly;
			}
			else if(A_xy(i,1) >= Ly){
				A_xy(i,1) -= Ly;
			}
		}

	return;
}
