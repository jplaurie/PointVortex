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

const double PI = 3.14159265358979323846;

struct Vortex {
    double x, y;
    double gamma;
    double u,v;
};


std::map<std::string, std::string> loadParams(const std::string& filename) {
    std::ifstream file(filename);
    std::map<std::string, std::string> params;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key, eq, value;
        if (iss >> key >> eq >> value && eq == "=") {
            params[key] = value;
        }
    }

    return params;
}


int main() {



    //  Initialize Random Seed
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed1);
    std::normal_distribution<double> distribution(0.0,1.0);


    //  Load Parameters 
    auto params = loadParams("params.txt");

    int N = std::stoi(params["num_vortices"]);
    double dt = std::stod(params["dt"]);
    int numSteps = std::stoi(params["num_steps"]);
    double coreSize = std::stod(params["core_size"]);
    int num_threads = std::stod(params["num_threads"]);

    std::cout << "Loaded " << N << " vortices.\n";
    // Now use these parameters in your simulation


    //  Setup OpenMP
    int nthreads = 0;
    if(num_threads == 0){
        nthreads = omp_get_max_threads();
    }
    else{
        nthreads = num_threads;
    }
    omp_set_num_threads(nthreads);


    std::cout << N << std::endl;
    std::vector<Vortex> vortices(N);

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
    for (int step = 0; step < numSteps; ++step) {
        
        
        
        // Compute velocities using OpenMP
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i == j) continue;

                double dx = vortices[i].x - vortices[j].x;
                double dy = vortices[i].y - vortices[j].y;
                double r2 = dx * dx + dy * dy + coreSize;

                double coeff = vortices[j].gamma / (2 * PI * r2);

                vortices[i].u+= coeff * (-dy);
                vortices[i].v += coeff * dx;
            }
            
        }

        // Update positions
        for (int i = 0; i < N; ++i) {
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