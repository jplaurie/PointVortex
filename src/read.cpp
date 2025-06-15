#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <cmath>
#include <vector>
#include "params.h"
#include "vortex.h"
#include "read.h"
#include <unordered_map>
#include <functional>

/*
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
*/

void initializeVortices(std::vector<Vortex>& vortices, SimParams params){

for (int i = 0; i < params.N; ++i) {
    double angle = 2 * params.PI * i / params.N;
    vortices[i].x = cos(angle);
    vortices[i].y = sin(angle);
    vortices[i].circ = (i % 2 == 0) ? 1.0 : -1.0;
}

return;

}



SimParams loadParams(const std::string& filename) {

    SimParams params;
    std::ifstream file(filename);
    std::string line;

    // Create a map of handlers
    std::unordered_map<std::string, std::function<void(const std::string&)>> setters;

    // Define how to handle each key
    setters["N"] = [&](const std::string& val) {
        params.N = std::stod(val);
    };
    setters["timeStep"] = [&](const std::string& val) {
        params.timeStep = std::stod(val);
    };
    setters["numSteps"] = [&](const std::string& val) {
        params.numSteps = std::stoi(val);
    };


    // Read the file
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key, value;
        if (!(iss >> key >> value)) continue;

        auto it = setters.find(key);
        if (it != setters.end()) {
            it->second(value); // call the setter function
        } else {
            std::cerr << "Unknown parameter: " << key << "\n";
        }
    }


    return params;
}