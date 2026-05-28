#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <numbers>

#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include <random>
#include <chrono>
#include "params.h"
#include "vortex.h"
#include "read.h"
#include "print.h"

int main() {

    //  Initialize Random Seed
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed1);
    std::normal_distribution<double> distribution(0.0,1.0);

    //  Load Parameters 
    auto params = loadParams("params.txt");

    double dt = params.timeStep;
    double runTime = 0.0;
    int fileNumber = 0;
    
    
    std::cout << "Loaded " << params.N << " vortices.\n";
    // Now use these parameters in your simulation


    //  Setup OpenMP
    int nthreads = 0;
    if(params.numThreads == 0){
        nthreads = omp_get_max_threads();
    }
    else{
        nthreads = params.numThreads;
    }
    omp_set_num_threads(nthreads);


    std::cout << params.N << std::endl;
    VortexSystem vortices(params.N);

    // Initialize vortices in a circle with alternating signs
    initializeVortices(vortices, params);


    //  Timestep Loop
    for (int step = 0; step < params.numSteps; ++step) {
        

         // main time stepping
    //   vortices = rungeKutta45(vortices, dt);
  
        runTime += dt;  //Add dt to runTime

        
        
       

       if( int(runTime / params.OutputTime) == fileNumber + 1){

            fileNumber++;
         //   printVortex(vortices, runTime, fileNumber);

         //   computeMomentum(vortices, momentumxValue, momentumyValue, momentumangularValue);

          //  computeHamiltonian(vortices, hamiltonianValue);

          //  printDiagnostics(hamiltonianValue, momentumxValue, momentumyValue, momentumangularValue, runTime,fileNumber);

         //   if(FLAG_RECORD_ENERGY == true){
          //      printRemoval(runTime,fileNumber,removeLowerCounter,removeUpperCounter, energyLowerRemoval, energyUpperRemoval,energyInjection);
          //  }
        
          //  std::cout << "file = " << fileNumber << " time = " << runTime << " dt = " << dt << " Ham = " << hamiltonianValue << " MomX = " << momentumxValue << " MomY = " << momentumyValue << " MomAngular = " << momentumangularValue << endl;

          //  std::cout << "Step " << step << " completed\n";
        }
    }

    // Output final vortex positions
 /*   for (const auto& v : vortices) {
        std::cout << v.x << " " << v.y << "\n";
    }
*/
    return 0;
}