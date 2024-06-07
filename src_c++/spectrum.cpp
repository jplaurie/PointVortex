#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include <time.h>
#include <fstream>
#include "Const.h"
#include "Functions.h"

using namespace std;
using namespace arma;

static rowvec Ek(2.0*Nk), Tk(2.0*Nk);//, Ek2(2.0*Nk),Tk2(2.0*Nk);
static blitz::Array<int,1> count_k(2.0*Nk);
static double lij,xij,yij,k,self_energy, energy_tot, energyk_tot,rms;

void output_spectrum(blitz::Array<double,2> vortex_xy,blitz::Array<double,1> vortex_G,blitz::Array<double,2> vortex_vel, blitz::Array<double,1> Ek_avg,blitz::Array<double,1> Tk_avg,  double runtime, int counter, int &avg_counter){


Ek = 0.0;
//Ek2=0.0;
Tk=0.0;
// Tk2=0.0;
//Ek += (double) N;
lij=0.0;
xij=0.0;
yij=0.0;
k=0.0;
energy_tot = 0.0;
energyk_tot = 0.0;


for(int k=0; k<2*Nk ; k++){
count_k(k) = 0;
}
self_energy = 0.0;

for(int i =0; i < N; i++){
     self_energy += pow(vortex_G(i),2.0);
   }

for(int kx =-Nk; kx < Nk; kx++){
  for(int ky =-Nk; ky < Nk; ky++){
      
      k = (int) pow( (double) kx*kx + ky*ky  ,0.5);
   
      for(int i =0; i < N-1; i++){
          for(int j =i+1; j < N; j++){

           xij = vortex_xy(i,0) - vortex_xy(j,0);
           yij = vortex_xy(i,1) - vortex_xy(j,1);
           Ek(k) += 2.0*vortex_G(i)*vortex_G(j)*cos(kx*xij + ky*yij) ;
           Tk(k) -= 2.0*vortex_G(i)*vortex_G(j)*sin(kx*xij + ky*yij)*(kx*(vortex_vel(i,0) - vortex_vel(j,0)) + ky*(vortex_vel(i,1) - vortex_vel(j,1)));
        }
      }
      count_k(k) = count_k(k) + 1; 
      Ek(k) += self_energy;
      

}
}

for(int k =1; k < 2*Nk; k++){

    if(count_k(k) != 0){
      Ek(k) /= 4.0 * pi * (double) count_k(k) * k;
      Tk(k) /= 4.0 * pi * (double) count_k(k) * k;
    }

    energy_tot += Ek(k);
     energyk_tot += Ek(k) * k;

}



//cout << "k_E = " << energyk_tot/ energy_tot << endl;//" " <<  energyk_tot << " " <<  energy_tot << endl;
//========================================  Analytical 1D spectra ============================================
//Ek=0.0;

/*for(int k =0; k < Nk; k++){
    for(int i =0; i < N-1; i++){
      for(int j =i+1; j < N; j++){

        lij = pow( pow( vortex_xy(i,0) - vortex_xy(j,0), 2.0)  + pow(vortex_xy(i,1) - vortex_xy(j,1),2.0) ,0.5);
       // cout << lij << endl;

          Ek2(k) += 2.0*vortex_G(i)*vortex_G(j)*j0( (double) k*lij);
          Tk2(k) -= 2.0*vortex_G(i)*vortex_G(j)*j1( (double) k*lij)*k*pow( pow( vortex_vel(i,0) - vortex_vel(j,0), 2.0)  + pow(vortex_vel(i,1) - vortex_vel(j,1),2.0) ,0.5);
        }

      }

    }

 for(int k =0; k < Nk; k++){
 for(int i =0; i < N; i++){
 Ek2(k) += pow(vortex_G(i),2.0);
}
}
 for(int k =1; k < Nk; k++){

  Ek2(k) /= 4.0*pi*k;
  Tk2(k) /= 4.0*pi*k;
}
*/


		Ek_avg += Ek;
		Tk_avg += Tk;
		avg_counter++;



		ostringstream out_spec;
    	out_spec << "./Output/Energy_spec." << setw(5) << setfill('0') << counter << ends;                  //creates file name for outputting da    ta at time slice
    	string filename_spec = out_spec.str();
    	ofstream fout_spec(filename_spec.c_str());
    	fout_spec << scientific;
      fout_spec.precision(12);

    	for(int i =0; i < Nk; i++){
   		fout_spec << i << " " << Ek(i) << " " << (Ek_avg(i)/(double) avg_counter) << " " << Tk(i) << " " << (Tk_avg(i)/(double) avg_counter) << endl;
   		}
     
     	fout_spec.close();


     	ostringstream out_scale;
    	out_scale << "./Output/Energy_scale." << setw(5) << setfill('0') << counter << ends;                  //creates file name for outputting da    ta at time slice
    	string filename_scale = out_scale.str();
    	ofstream fout_scale(filename_scale.c_str());
    	fout_scale << scientific;
      fout_scale.precision(12);

    	
   		fout_scale << runtime << " " << (energyk_tot/ energy_tot) << " " <<  energyk_tot << " " <<  energy_tot << endl;
   		
     
     	fout_scale.close();
	



	return;
}

void output_vel(mat vortex_vel, double runtime, int counter){


//===============  compute rms of the velocity ===========================
rms = 0.0;

for(int i=0; i < N; i++){

    rms += pow( pow(vortex_vel(i,0),2.0) + pow(vortex_vel(i,1),2.0) , 0.5);

}
rms /= (double) N;


      ostringstream out_vel;
      out_vel << "./Output/rms_velocity." << setw(5) << setfill('0') << counter << ends;                  //creates file name for outputting da    ta at time slice
      string filename_vel = out_vel.str();
      ofstream fout_vel(filename_vel.c_str());
      fout_vel << scientific;
      fout_vel.precision(12);

      
      fout_vel << runtime << " " << rms << endl;
      
     
      fout_vel.close();


  return;
}

