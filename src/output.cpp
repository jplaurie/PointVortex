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

static double rms,dxs,dys,drs,di, st_norm;
static mat vortex_vel(N,2,fill::zeros);

void output_vortex(mat vortex_xy, rowvec vortex_G, double runtime, int counter){


		ostringstream out_vortex;
    	out_vortex << "./data/vortex_xy." << setw(5) << setfill('0') << counter << ends;                  //creates file name for outputting da    ta at time slice
    	string filename_vortex = out_vortex.str();
    	ofstream fout_vortex(filename_vortex.c_str());
    	fout_vortex.precision(12);
 		fout_vortex << scientific;
    	for(int i =0; i < N; i++){
   		   //fout_vortex << runtime << " " << vortex_xy(i,0) << " " << vortex_xy(i,1) << " " << vortex_G(i) << endl;
		   fout_vortex << vortex_xy(i,0) << " " << vortex_xy(i,1) << " " << vortex_G(i) << endl;
   		}
     
     	fout_vortex.close();

     	ofstream fout_count("./data/curframe.dat");
        fout_count << scientific;
        fout_count.precision(12);
        fout_count << runtime << " " << counter << endl;

	return;
}


void output_vel(mat vortex_vel, double runtime, int counter){


//===============  compute rms of the velocity ===========================
rms = 0.0;

for(int i=0; i < N; i++){

    rms += pow( pow(vortex_vel(i,0),2.0) + pow(vortex_vel(i,1),2.0) , 0.5);

}
rms /= double(N);


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


void output_structure(mat vortex_xy, rowvec vortex_G, mat structure_function, int counter, irowvec structure_counter){


vortex_vel = step_vortices(vortex_xy, vortex_G);
st_norm=0.0;
for(int i =0; i < N; i++){
	for(int j =i+1; j < N; j++){

		dxs = vortex_xy(j,0) - vortex_xy(i,0);
		dys = vortex_xy(j,1) - vortex_xy(i,1);

    if(dxs >= Lx/2) dxs -= Lx;
    else if(dxs < -Lx/2) dxs += Lx;
    if(dys >= Lx/2) dys -= Ly;
    else if(dys < -Ly/2) dys += Ly;

		drs = pow( dxs*dxs + dys*dys  ,0.5);


		di = int(drs*double(NBIN)/Lx);
		structure_function(di,0) = pow( (vortex_vel(j,0) - vortex_vel(i,0))*dxs + (vortex_vel(j,1) - vortex_vel(i,1))*dys,3.0);
		structure_function(di,1) = ((vortex_vel(j,0) - vortex_vel(i,0))*dxs + (vortex_vel(j,1) - vortex_vel(i,1))*dys)*pow(vortex_G(j)-vortex_G(i),2.0);
		structure_function(di,2) += pow( (vortex_vel(j,0) - vortex_vel(i,0))*dxs + (vortex_vel(j,1) - vortex_vel(i,1))*dys,3.0);
		structure_function(di,3) += ((vortex_vel(j,0) - vortex_vel(i,0))*dxs + (vortex_vel(j,1) - vortex_vel(i,1))*dys)*pow(vortex_G(j)-vortex_G(i),2.0);
		structure_counter(di) += 1;
	}
}


ostringstream out_struct;
out_struct << "./Output/structure." << setw(5) << setfill('0') << counter << ends;                  //creates file name for outputting da    ta at time slice
string filename_struct = out_struct.str();
ofstream fout_struct(filename_struct.c_str());
fout_struct << scientific;
fout_struct.precision(12);

for(int i =0; i < NBIN; i++){
	if(structure_counter(i) < 1){
		st_norm=1.0;
	}
	else{
	 st_norm = double(structure_counter(i));
	}
	fout_struct << i*Lx/double(NBIN) << " " << structure_function(i,0) << " " << structure_function(i,1) << " " << structure_function(i,2)/st_norm << " " << structure_function(i,3)/st_norm << endl;
}
fout_struct.close();


  return;
}


void output_remove(double runtime, int counter,int& remove_low_counter, int& remove_high_counter, double& energy_low_diss_counter, double & energy_high_diss_counter,double& energy_low_inject_counter, double & energy_high_inject_counter ){

      ostringstream out_remove;
      out_remove << "./Output/remove." << setw(5) << setfill('0') << counter << ends;                  //creates file name for outputting da    ta at time slice
      string filename_remove = out_remove.str();
      ofstream fout_remove(filename_remove.c_str());
      fout_remove << scientific;
      fout_remove.precision(12);

      fout_remove << runtime << " " << remove_low_counter << " " << remove_high_counter << " " << energy_low_diss_counter << " " << energy_high_diss_counter << " " << energy_low_inject_counter << " " << energy_high_inject_counter << endl;
       fout_remove.close();
      
      remove_low_counter=0;
      remove_high_counter=0;
      energy_low_diss_counter=0.0;
      energy_high_diss_counter=0.0;
      energy_low_inject_counter=0.0;
      energy_high_inject_counter=0.0;


return;


}

