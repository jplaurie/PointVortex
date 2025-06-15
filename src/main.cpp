#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include <random>
#include <chrono>
#include "params.h"
#include "vortex.h"
#include "read.h"

int main() {



    //  Initialize Random Seed
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed1);
    std::normal_distribution<double> distribution(0.0,1.0);


    //  Load Parameters 
    auto params = loadParams("params.txt");


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
    std::vector<Vortex> vortices(params.N);

    // Initialize vortices in a circle with alternating signs
    initializeVortices(vortices);
    /*for (int i = 0; i < N; ++i) {
        double angle = 2 * PI * i / N;
        vortices[i].x = cos(angle);
        vortices[i].y = sin(angle);
        vortices[i].gamma = (i % 2 == 0) ? 1.0 : -1.0;
    }
    */



    //  Timestep Loop
    for (int step = 0; step < params.numSteps; ++step) {
        
        
        
        // Compute velocities using OpenMP
        #pragma omp parallel for
        for (int i = 0; i < params.N; ++i) {
            for (int j = 0; j < params.N; ++j) {
                if (i == j) continue;

                double dx = vortices[i].x - vortices[j].x;
                double dy = vortices[i].y - vortices[j].y;
                double r2 = dx * dx + dy * dy + params.coreSize;

                double coeff = vortices[j].gamma / (2 * params.PI * r2);

                vortices[i].u+= coeff * (-dy);
                vortices[i].v += coeff * dx;
            }
            
        }

        // Update positions
        for (int i = 0; i < params.N; ++i) {
            vortices[i].x += dt * vortices[i].u;
            vortices[i].y += dt * vortices[i].v;
        }

        // Optional: print status
        if (step % 100 == 0) {
            std::cout << "Step " << step << " completed\n";
        }
    }

    // Output final vortex positions
    for (const auto& v : vortices) {
        std::cout << v.x << " " << v.y << "\n";
    }

    return 0;
}