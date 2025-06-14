#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

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
