
#include <iostream>
#include <cmath>
#include <string>
#include <eigen3/Eigen/Dense>
#include <complex>
#include <time.h>
#include <fstream>
#include "Const.h"
#include "Functions.h"
#include <random>
#include <chrono>
#include <omp.h>


using namespace std;
using namespace Eigen;



int main(){

double startTime = 0;
double finishTime = 0;
startTime = omp_get_wtime();


//============== Setup OpenMP ==============
int nthreads = 0;
if(num_threads == 0){
    nthreads = omp_get_max_threads();
}
else{
    nthreads = num_threads;
}
omp_set_num_threads(nthreads);

//============== Array Declarations ==============
MatrixXd vortex_xy(N,2), vortex_vel(N,2), eta(N,2),structure_function(NBIN,4);
VectorXd vortex_G(N);
VectorXi structure_fileNumber(NBIN);

double dt = dtMax;
int fileNumber = 0;
//int avg_fileNumber = 0;
double runTime = 0.0;
int removeLowerCounter = 0;
int removeUpperCounter = 0;
double energyLowerRemoval = 0.0;
double energyUpperRemoval = 0.0;
double energyInjection = 0.0;


double hamiltonianValue = 0.0;
double momentumxValue = 0.0;
double momentumyValue = 0.0;

system("mkdir ./output/"); //Create folder for diagonistics
 

//============== Initialize Random Seed==============
unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
mt19937 generator(seed1);
normal_distribution<double> distribution(0.0,1.0);

//============== Setup Initial State ==============
printParameters();
initializeVortices(vortex_xy, vortex_G, runTime, fileNumber);

//============== Begin Timestepping Loop ==============
for(int timeStep = 1; timeStep <= totalSteps; timeStep++){
    
    if(FLAG_NOISE == "gaussian"){
        for(int i = 0; i < N; i++){
            eta(i,0) = distribution(generator);
            eta(i,1) = distribution(generator);
        }
        vortex_xy += pow(2.0* noiseStrength * dt,0.5) * eta;
    }

    // main time stepping
    vortex_xy = rungeKutta45(vortex_xy, vortex_G, dt);
  
    runTime += dt;  //Add dt to runTime

    //Apply boundary conditions
    if(FLAG_BOUNDARY != "infinite"){
        invokeBoundaryCondtions(vortex_xy); // invoke the periodic boundary conditions
    }

    // removes oppositely signed vortices and place then randomly back into the box
    if(FLAG_VORTEX_REMOVE == true){
        removeVortices(vortex_xy, vortex_G,removeLowerCounter,removeUpperCounter,energyLowerRemoval,energyUpperRemoval,energyInjection);		
    }
    


    if( int(runTime / outputTime) == fileNumber + 1){

        fileNumber++;
        printVortex(vortex_xy,vortex_G,runTime, fileNumber);

       // if(FLAG_STRUCTURE == true){
       // 	output_structure(vortex_xy,vortex_G,structure_function, fileNumber,structure_fileNumber);
       // }

       //output_vel(vortex_vel,runTime,fileNumber);

       
        if(FLAG_DIAGNOSTICS == true){ 
            computeMomentum(vortex_xy, vortex_G, momentumxValue, momentumyValue);

            computeHamiltonian(vortex_xy, vortex_G, hamiltonianValue);

            printDiagnostics(hamiltonianValue, momentumxValue, momentumyValue,runTime,fileNumber);
        }

        if(FLAG_RECORD_ENERGY == true){
             printRemoval(runTime,fileNumber,removeLowerCounter,removeUpperCounter, energyLowerRemoval, energyUpperRemoval,energyInjection);
        }
        
        cout << "file = " << fileNumber << " time = " << runTime << " dt = " << dt << " Ham = " << hamiltonianValue << " MomX = " << momentumxValue << " MomY = " << momentumyValue << endl;
    }
}

finishTime = omp_get_wtime();
cout << "time taken for code is " << finishTime-startTime << endl;

return 0;

}

