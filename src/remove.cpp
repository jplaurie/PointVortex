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

static double dist,dx,dy, phase, dist_closest,prob,ham_old,ham_new,temp, kappa_temp_i,kappa_temp_j;
static int j_closest;


void remove(mat & A_xy,rowvec & A_G, int& remove_low_counter, int& remove_high_counter, double & energy_low_diss_counter, double & energy_high_diss_counter,double & energy_low_inject_counter, double & energy_high_inject_counter){

phase = 0.0;
j_closest=0;
dist_closest=0.0;
dist = 0.0;
dx=0.0;
dy=0.0;

for(int i = 0; i < N; i++){

	dist_closest = 2.0*fmax(Lx,Ly);

	for(int j = 0; j < N; j++){

		if( signbit(A_G(i)/A_G(j)) == true){
			dx = A_xy(i,0) - A_xy(j,0);
			dy = A_xy(i,1) - A_xy(j,1);
			if(dx  >= 0.5*Lx){
				dx -= Lx;
			}
			else if(dx < -0.5*Lx){
				dx += Lx;
			}
			if(dy  >= 0.5*Ly){
				dy -= Ly;
			}
			else if(dy < -0.5*Ly){
				dy += Ly;
			}
			

			dist = pow(dx*dx + dy*dy,0.5);

			if(dist < dist_closest){
				dist_closest = dist;
				j_closest = j;
			}
		}
	}
	//cout << "dist_closest = " << dist_closest << endl;
	if(dist_closest < remove_dist_low){
		remove_low_counter++;
		if(FLAG_ADD == false){
			A_G(i) = 0.0;
			A_G(j_closest)=0.0;
		}
		else{


		ham_old=0.0;
		ham_new=0.0;
		hamiltonian(A_xy, A_G, ham_old, temp,temp);
		kappa_temp_i = A_G(i);
		kappa_temp_j = A_G(j_closest);
		A_G(i)=0.0;
		A_G(j_closest)=0.0;
		hamiltonian(A_xy, A_G, ham_new, temp,temp);
		energy_low_diss_counter += ham_new-ham_old;
		A_G(i)=kappa_temp_i;
		A_G(j_closest)=kappa_temp_j;

		//====reposition point i =========================
			A_xy(i,0) = Lx*double(rand())/(RAND_MAX + 1.0);
			A_xy(i,1) = Ly*double(rand())/(RAND_MAX + 1.0);

			phase = 2.0*pi*double(rand())/(RAND_MAX + 1.0);
			A_xy(j_closest,0) = A_xy(i,0) + L*cos(phase);
			A_xy(j_closest,1) =  A_xy(i,1) + L*sin(phase);

			if(A_xy(j_closest,0)  >= Lx) A_xy(j_closest,0) -= Lx;
			else if(A_xy(j_closest,0)  < 0) A_xy(j_closest,0) += Lx;

			if(A_xy(j_closest,1)  >= Ly) A_xy(j_closest,1) -= Ly;
			else if(A_xy(j_closest,1)  < 0) A_xy(j_closest,1) += Ly;


		hamiltonian(A_xy, A_G, ham_old, temp,temp);
		
		energy_low_inject_counter += ham_old-ham_new;	

		}
	}
	else if(dist_closest > remove_dist_high){
		prob = double(rand())/(RAND_MAX + 1.0);
		
		if(prob < prob_remove){
		
			remove_high_counter++;
			if(FLAG_ADD == false){
				A_G(i) = 0.0;
				A_G(j_closest)= 0.0;
			}
			else{
				ham_old=0.0;
				ham_new=0.0;
				hamiltonian(A_xy, A_G, ham_old, temp,temp);
				kappa_temp_i = A_G(i);
				kappa_temp_j = A_G(j_closest);
				A_G(i)=0.0;
				A_G(j_closest)=0.0;
				hamiltonian(A_xy, A_G, ham_new, temp,temp);
				energy_high_diss_counter += ham_new-ham_old;
				A_G(i)=kappa_temp_i;
				A_G(j_closest)=kappa_temp_j;



				//====reposition point i =========================
				A_xy(i,0) = Lx*double(rand())/(RAND_MAX + 1.0);
				A_xy(i,1) = Ly*double(rand())/(RAND_MAX + 1.0);

				phase = 2.0*pi*double(rand())/(RAND_MAX + 1.0);
				A_xy(j_closest,0) = A_xy(i,0) + L*cos(phase);
				A_xy(j_closest,1) =  A_xy(i,1) + L*sin(phase);

				if(A_xy(j_closest,0)  >= Lx) A_xy(j_closest,0) -= Lx;
				else if(A_xy(j_closest,0)  < 0) A_xy(j_closest,0) += Lx;

				if(A_xy(j_closest,1)  >= Ly) A_xy(j_closest,1) -= Ly;
				else if(A_xy(j_closest,1)  < 0) A_xy(j_closest,1) += Ly;
			
				hamiltonian(A_xy, A_G, ham_old, temp,temp);

				energy_high_inject_counter += ham_old-ham_new;

			}
		}
	}	
}

	return;
}








