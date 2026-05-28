// params.h
#ifndef PARAMS_H
#define PARAMS_H

struct SimParams {
    int N;
    double PI = 3.14159265358979323846; 
    double timeStep;
    double OutputTime;
    int numSteps;
    double coreSize;
    int numThreads;
    std::string boundaryCondition;
};

#endif