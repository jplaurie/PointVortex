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

static double dist,dx,dy, phase, dist_closest,hamiltonValueBefore,hamiltonValueAfter, Gi,Gj,r;
static int jClosest;
static irowvec neighbour(N,fill::value(-1)), removeCandidate(N,fill::value(0));
static rowvec nearestDistance(N,fill::zeros);
static bool FLAG_RECURSION, FLAG_VORTEX_INSIDE;

double distance(mat, int ,int );
void nearestNeighbour(mat , irowvec & );





void removeVortices(mat & A_xy,rowvec & A_G, int& removeLowerCounter, int& removeUpperCounter, double & energyLowerRemoval, double & energyUpperRemoval,double & energyInjection){




phase = 0.0;
jClosest=0;
dist_closest=0.0;
dist = 0.0;
dx=0.0;
dy=0.0;


removeCandidate.zeros();

nearestNeighbour(A_xy, removeCandidate);

//record Hamiltonian before deletion of vortices
if(FLAG_RECORD_ENERGY == true){ 
	computeHamiltonian(A_xy, A_G, hamiltonValueBefore);
}

//delete vortices that are too close
for(int i = 0; i < N; i++){
	if( removeCandidate(i) == -1 ){
		//sets circulation of removed vortices to zero
		A_G(i) = 0.0;
		removeLowerCounter++;
	}
}
//record energy difference						
if(FLAG_RECORD_ENERGY == true){ 
	computeHamiltonian(A_xy, A_G, hamiltonValueAfter);
	energyLowerRemoval += hamiltonValueAfter-hamiltonValueBefore;
}	
//delete vortices that are too far
for(int i = 0; i < N; i++){
	if( removeCandidate(i) == 1 ){
		//sets circulation of removed vortices to zero
		A_G(i) = 0.0;
		removeUpperCounter++;
	}
}
//record energy difference						
if(FLAG_RECORD_ENERGY == true){ 
	computeHamiltonian(A_xy, A_G, hamiltonValueBefore);
	energyUpperRemoval += hamiltonValueBefore-hamiltonValueAfter;
}	

if(FLAG_VORTEX_ADD == true){
	//restore values if circulations for reinjection
	for(int i = 0; i < Np; i++){
		if( abs(removeCandidate(i)) == 1){
			A_G(i) = 1.0/double(N);
		}
	}
	for(int i = Np; i < N; i++){
		if( abs(removeCandidate(i)) == 1){
			A_G(i) = -1.0/double(N);
		}
	}

	for(int i = 0; i< N; i++){
		if( abs(removeCandidate(i)) == 1){

			if(FLAG_BOUNDARY == "periodic"){

				//Re-position vortex i 
				A_xy(i,0) = Lx * ( double(rand())/(RAND_MAX + 1.0) - 0.5);
				A_xy(i,1) = Ly * ( double(rand())/(RAND_MAX + 1.0) - 0.5 );
		
				phase = 2.0*pi*double(rand())/(RAND_MAX + 1.0);
				A_xy(neighbour(i),0) = A_xy(i,0) + L*cos(phase);
				A_xy(neighbour(i),1) =  A_xy(i,1) + L*sin(phase);

				if(A_xy(neighbour(i),0)  >= Lx/2.0) A_xy(neighbour(i),0) -= Lx;
				else if(A_xy(neighbour(i),0)  < -Lx/2.0) A_xy(neighbour(i),0) += Lx;

				if(A_xy(neighbour(i),1)  >= Ly/2.0) A_xy(neighbour(i),1) -= Ly;
				else if(A_xy(neighbour(i),1)  < -Ly/2.0) A_xy(neighbour(i),1) += Ly;
			}
			else if(FLAG_BOUNDARY == "disc"){

				//Re-position vortex i 
				r = R * double(rand())/(RAND_MAX + 1.0);
				phase = 2.0*pi*double(rand())/(RAND_MAX + 1.0);
				A_xy(i,0) = r*cos(phase);
				A_xy(i,1) = r*sin(phase);
				FLAG_VORTEX_INSIDE = false;
				
				while(FLAG_VORTEX_INSIDE == false){
					phase = 2.0*pi*double(rand())/(RAND_MAX + 1.0);
					A_xy(neighbour(i),0) = A_xy(i,0) + Lr*cos(phase);
					A_xy(neighbour(i),1) =  A_xy(i,1) + Lr*sin(phase);

					r = pow( pow(A_xy(neighbour(i),0),2.0) + pow(A_xy(neighbour(i),1),2.0) , 0.5);
					
					if( r < R ){
						FLAG_VORTEX_INSIDE = true;
					}
				}



			} 
			removeCandidate(i) = 0;
			removeCandidate(neighbour(i)) = 0;
		}
	}

	if(FLAG_RECORD_ENERGY == true){
		computeHamiltonian(A_xy, A_G, hamiltonValueAfter);	
		energyInjection += hamiltonValueAfter-hamiltonValueBefore;
	}
}	

return;
}



//Computes the distance between two vortices
double distance(mat A_xy, int ii, int jj){
	
	dx = A_xy(ii,0) - A_xy(jj,0);
	dy = A_xy(ii,1) - A_xy(jj,1);
	
	if(FLAG_BOUNDARY == "periodic"){
		//adjust distance due to periodic boundaries
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
	}
	else if(FLAG_BOUNDARY == "periodic-x"){
		//adjust distance due to periodic boundaries
		if(dx  >= 0.5*Lx){
			dx -= Lx;
		}
		else if(dx < -0.5*Lx){
			dx += Lx;
		}

	}
			
	//return distance between vortices
	return pow(dx*dx + dy*dy,0.5);
}



void nearestNeighbour(mat A_xy, irowvec & removeCandidate){

	//iterate over positive vortices checking closest negative vortex
	for(int i = 0; i < Np; i++){
		if(removeCandidate(i) == 0 ){
			//set arbitrary large distance (this will always be written over) 
			if(FLAG_BOUNDARY == "periodic"){
				dist_closest = 2*(Lx+Ly); 
			}
			else if(FLAG_BOUNDARY == "disc"){
				dist_closest = 4 * R; 
			}
			
			//find the closest vortex (j) of opposite sign if vortex not already candidate for removal
			for(int j = Np; j < N; j++){
				
				if(removeCandidate(j) == 0){
					//calculate distance
					dist = distance(A_xy,i,j);

					if(dist < dist_closest){
						dist_closest = dist;
						jClosest = j;
					}
				}
			}
			//record distance and index of nearest neighbour
			neighbour(i)  = jClosest;
			nearestDistance(i) = dist_closest;
		}
	}
	//iterate over negative vortices checking closest positive vortex
	for(int i = Np; i < N; i++){
		if(removeCandidate(i) == 0){
			//set arbitrary large distance (this will always be written over)
			dist_closest = 2*(Lx+Ly); 
			//find the closest vortex (j) of opposite sign if vortex not already candidate for removal
			for(int j = 0; j < Np; j++){
				
				if(removeCandidate(j) == 0){
					dist = distance(A_xy,i,j);

					if(dist < dist_closest){
						dist_closest = dist;
						jClosest = j;
					}
				}
			}
			//record distance and index of nearest neighbour
			neighbour(i)  = jClosest;
			nearestDistance(i) = dist_closest;
		}
	}
	
	//set flag for recursion of function
	FLAG_RECURSION = false;

	//if nearest neighbours select each other (pairs) set for removal, if not recall function taking these pairs out
	for(int i = 0; i < N; i++){

		if(removeCandidate(i) == 0){
			if( nearestDistance(i) < removalDistanceLower ){
				if( i == neighbour(neighbour(i)) ){
					removeCandidate(i) = -1;
				}
				else{
					FLAG_RECURSION = true;
				} 
			}
			else if( nearestDistance(i) > removalDistanceUpper ){
				if( i == neighbour(neighbour(i)) ){
					removeCandidate(i) = 1;
				}
				else{
					FLAG_RECURSION = true;
				}
			}
		}
	}
	if( FLAG_RECURSION== true){ //call function again if there are triples of vortices close to each other
		nearestNeighbour(A_xy, removeCandidate);
	}
	else{
		return;
	}
		
}




