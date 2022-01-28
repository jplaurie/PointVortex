
#include <iostream>
#include <cmath>
#include <armadillo>
#include <time.h>
#include <fstream>
#include "Const.h"
#include "Functions.h"
#include <iomanip>

using namespace std;
using namespace arma;

static double dipole_dist, dipole_shift;

void random_4_vortex(mat &, rowvec &);



void initialize_vortices(mat & A_xy, rowvec & A_G, double &runtime, int &out_count){


if(FLAG_INITIAL_CONDITION == 0){

fstream filein("./data/curframe.dat");
filein.precision(12);
filein >> runtime >> out_count;

ostringstream in_data;
in_data << "./data/vortex_xy." << setw(5) << setfill('0') << out_count << ends;
string filename = in_data.str();
ifstream filein2(filename.c_str());
filein2.precision(12);

cout << "loading data in file vortex_xy." << setw(5) << setfill('0') << out_count << endl;

for(int i = 0; i < N; i++){
	filein2 >> A_xy(i,0) >> A_xy(i,1) >> A_G(i);
}

}
else if(FLAG_INITIAL_CONDITION == 1){
	
	out_count = 0;
	runtime = 0.0;

	unsigned seed2 = std::chrono::system_clock::now().time_since_epoch().count();
		mt19937 generator2(seed2);
		uniform_real_distribution<double> distribution2(0.0,1.0);
	//============RANDOM =============================================
		for(int i = 0; i < N; i++){
			A_xy(i,0) = Lx*distribution2(generator2);
			A_xy(i,1) = Ly*distribution2(generator2);	
		}

		for(int i = 0; i < Np; i++){
			A_G(i)= 1.0 /  double(N);
		}
		for(int i = 0; i < Nm; i++){

			A_G(i+ Np) = -1.0 / double(N);
		}

		ofstream fout_write("./data/vortex_xy.00000");
		fout_write << scientific;
		fout_write.precision(12);

		for(int i = 0; i < N; i++){
			fout_write << A_xy(i,0) << " " << A_xy(i,1) << " " << A_G(i) << endl;
		}

		fout_write.close();
}
else if(FLAG_INITIAL_CONDITION == 2){

cout << "Reading initial condition from file..." << endl;

out_count = 0;
runtime = 0.0;
dipole_dist=0.0;
dipole_shift=0.0;

fstream filein("./data/dipole_data.dat");
filein.precision(12);
filein >> dipole_dist >> dipole_shift;

A_xy(0,0) = 0.0;
A_xy(0,1) = 0.0;
A_G(0)= 1.0 /  double(N);

for(int i = 1; i < 3; i++){
		A_xy(i,0) = -pi;
		A_xy(i,1) = 0.0 + 0.5*dipole_dist*pow(-1.0,double(i)) + dipole_shift;	
		A_G(i)= pow(-1.0,double(i)) /  double(N);
}

	ofstream fout_read("./data/vortex_xy.00000");
	fout_read << scientific;
	fout_read.precision(12);

	for(int i = 0; i < N; i++){
		fout_read << A_xy(i,0) << " " << A_xy(i,1) << " " << A_G(i) << endl;
	}

	fout_read.close();



}
else if(FLAG_INITIAL_CONDITION == 3){

random_4_vortex(A_xy, A_G);

out_count = 0;
runtime = 0.0;

ofstream fout_write("./data/vortex_xy.00000");
	fout_write << scientific;
	fout_write.precision(12);

	for(int i = 0; i < N; i++){
		fout_write << A_xy(i,0) << " " << A_xy(i,1) << " " << A_G(i) << endl;
	}

	fout_write.close();

	out_count = 0;
	runtime = 0.0;

}
else{


	cout << "no defined initialization routine definded, exiting code" << endl;
		
		exit(1);

}
}


void random_4_vortex(mat & A_xy,rowvec & A_G){

unsigned seed2 = std::chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator2(seed2);
	uniform_real_distribution<double> distribution2(0.0,1.0);

gotostart:

double l =0.0;

for(int i = 0; i < N; i++){
		A_xy(i,0) = Lx*distribution2(generator2);
		A_xy(i,1) = Ly*distribution2(generator2);	
	}

for(int i = 0; i < Np; i++){
		A_G(i)= 1.0 /  double(N);
	}
for(int i = 0; i < Nm; i++){
	A_G(i+ Np) = -1.0 / double(N);
}

//compute lengths

for(int i = 0; i < N; i++){
	for(int j=i+1; j< N; j++){

	l= pow(pow(A_xy(i,0) - A_xy(j,0),2.0) + pow(A_xy(i,1) - A_xy(j,1),2.0),0.5);

	if(l < 0.8*Lx/2.0) goto gotostart;
	
	}
}
cout << "done!" << endl; 

return;
}



