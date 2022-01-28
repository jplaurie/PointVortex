
#include <iostream>
#include <cmath>
#include <armadillo>
#include <complex>
#include <time.h>
#include <fstream>
#include "Const.h"
#include "Functions.h"
#include <random>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace arma;

int main(){

double start = 0;
double finish = 0;
start = omp_get_wtime();

double lmin=0;
//================set up openmp ===============
int nthreads = 0;
if(num_threads == 0){
    nthreads = omp_get_max_threads();
}
else{
    nthreads = num_threads;
}
omp_set_num_threads(nthreads);

//====================================================
mat vortex_xy(N,2,fill::zeros), vortex_vel(N,2,fill::zeros), eta(N,2,fill::zeros),structure_function(NBIN,4,fill::zeros);
rowvec vortex_G(N,fill::zeros);
irowvec structure_counter(NBIN,fill::zeros);

double dt = maxdt;
int counter = 0;
//int avg_counter = 0;
double runtime =0.0;
int remove_low_counter = 0;
int remove_high_counter = 0;
double energy_low_diss_counter = 0.0;
double energy_high_diss_counter = 0.0;
double energy_low_inject_counter = 0.0;
double energy_high_inject_counter = 0.0;
double ham_out =0.0;
double momx_out = 0.0;
double momy_out = 0.0;

ofstream fout_data("./parameters.txt");
fout_data.precision(12);
fout_data << "Number of positive vortices = " << Np << endl;
fout_data << "Number of negative vortices = " << Nm << endl;
fout_data << "Total number of vortices =  " << N << endl;
fout_data << "Length in x = " << Lx << endl;
fout_data << "Length in y = " << Ly << endl;
fout_data << "Vortex density = " << N / (Lx*Ly) << endl;
fout_data << "mean intervortex spacing = " << sqrt(Lx*Ly/(double) N) << endl;
fout_data << "===========time stepping=========" << endl;
fout_data << "Runge-Kutta tolerance = " << rktol << endl;
fout_data << "Total number of steps =  " << nsteps << endl;
fout_data << "FLAG_OUTPUT = " << FLAG_OUTPUT << endl;
fout_data << "Output step = " << outstep << endl;
fout_data << "FLAG_PERIODIC_BC = " << FLAG_PERIODIC_BC << endl;
fout_data << "FLAG_STRUCTURE = " << FLAG_STRUCTURE << endl;
fout_data << "============forcing==============" << endl;
fout_data << "FLAG_NOISE = " << FLAG_NOISE << endl;
fout_data << "FLAG_REMOVE = " << FLAG_REMOVE << endl;
//fout_data << "FLAG_REMOVE_OPP_SIGN_ONLY = " << FLAG_REMOVE_OPP_SIGN_ONLY << endl;
fout_data << "FLAG_ADD = " << FLAG_ADD << endl;
fout_data << "Noise amplitude =  " << D << endl;
fout_data << "remove_dist_low = " << remove_dist_low << endl;
fout_data << "remove_dist_high = " << remove_dist_high << endl;
fout_data << "Reintroduce vortices at scale L = " << L << endl;
fout_data << "============parallelism============" << endl;
fout_data << "number of threads = " << num_threads << endl;

//output data in the terminal
cout << "Number of positive vortices = " << Np << endl;
cout << "Number of negative vortices = " << Nm << endl;
cout << "Total number of vortices =  " << N << endl;
cout << "Length in x = " << Lx << endl;
cout << "Length in y = " << Ly << endl;
cout << "Vortex density = " << N / (Lx*Ly) << endl;
cout << "mean intervortex spacing = " << sqrt(Lx*Ly/(double) N) << endl;
cout << "===========time stepping=========" << endl;
cout << "Runge-Kutta tolerance = " << rktol << endl;
cout << "Total number of steps =  " << nsteps << endl;
cout << "FLAG_OUTPUT = " << FLAG_OUTPUT << endl;
cout << "Output step = " << outstep << endl;
cout << "FLAG_PERIODIC_BC = " << FLAG_PERIODIC_BC << endl;
cout << "FLAG_STRUCTURE = " << FLAG_STRUCTURE << endl;
cout << "============forcing==============" << endl;
cout << "FLAG_NOISE = " << FLAG_NOISE << endl;
cout << "FLAG_REMOVE = " << FLAG_REMOVE << endl;
//cout << "FLAG_REMOVE_OPP_SIGN_ONLY = " << FLAG_REMOVE_OPP_SIGN_ONLY << endl;
cout << "FLAG_ADD = " << FLAG_ADD << endl;
cout << "Noise amplitude =  " << D << endl;
cout << "remove_dist_low = " << remove_dist_low << endl;
cout << "remove_dist_high = " << remove_dist_high << endl;
cout << "Reintroduce vortices at scale L = " << L << endl;
cout << "============parallelism============" << endl;
cout << "number of threads = " << num_threads << endl;


system("mkdir ./Output/");


//================initial random variables=======================
unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
mt19937 generator(seed1);
normal_distribution<double> distribution(0.0,1.0);

//=============set up initial condition======================
initialize_vortices(vortex_xy, vortex_G, runtime, counter);

//==================begin timestepping ==========================================
for(int timestep = 1; timestep <= nsteps; timestep++){

    if(FLAG_NOISE == true){

        for(int i = 0; i < N; i++){
            eta(i,0) = distribution(generator);
            eta(i,1) = distribution(generator);
        }
        vortex_xy += pow(2.0*D*dt,0.5)*eta;
    }

    vortex_xy = RK45(vortex_xy, vortex_G, vortex_vel, dt);


    runtime += dt;

    if(FLAG_PERIODIC_BC == true){
        invoke_BC(vortex_xy); // invoke the periodic boundary conditions
    }

    if(FLAG_REMOVE == true){
        remove(vortex_xy, vortex_G,remove_low_counter,remove_high_counter,energy_low_diss_counter,energy_high_diss_counter,energy_low_inject_counter,energy_high_inject_counter);		// removes oppositely signed vortices and place then randomly back into the box
    }

    if(timestep % outstep == 0){

        counter++;
        output_vortex(vortex_xy,vortex_G,runtime, counter);

        if(FLAG_STRUCTURE == true){
        	output_structure(vortex_xy,vortex_G,structure_function, counter,structure_counter);
        }

       //output_vel(vortex_vel,runtime,counter);


        if(FLAG_REMOVE == true){
        	output_remove(runtime,counter,remove_low_counter,remove_high_counter, energy_low_diss_counter, energy_high_diss_counter,energy_low_inject_counter, energy_high_inject_counter);
        }
        hamiltonian(vortex_xy, vortex_G, ham_out, momx_out, momy_out);
        output_hamiltonian(runtime, counter, ham_out, momx_out, momy_out);

        cout << "file = " << counter << " time = " << runtime << " Ham = " << ham_out << " dt = " << dt << endl;
	/*
	
	lmin = pow( pow(vortex_xy(3,0) - vortex_xy(0,0),2.0) + pow(vortex_xy(3,1)-vortex_xy(0,1),2.0) ,0.5);
	for(int ii=1; ii<3; ii++){
	if(lmin>pow( pow(vortex_xy(3,0) - vortex_xy(ii,0),2.0) + pow(vortex_xy(3,1)-vortex_xy(ii,1),2.0) ,0.5) ){
		lmin=pow( pow(vortex_xy(3,0) - vortex_xy(ii,0),2.0) + pow(vortex_xy(3,1)-vortex_xy(ii,1),2.0) ,0.5);
	}
	
	}

	if(lmin < 0.15*pi) return 0;


    }
*/
}
}

finish = omp_get_wtime();
cout << "time taken for code is " << finish-start << endl;

return 0;

}
