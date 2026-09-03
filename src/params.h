#ifndef POINT_VORTEX_PARAMS_H
#define POINT_VORTEX_PARAMS_H
#include <cstddef>
#include <cstdint>
#include <string>
enum class IntegratorKind { rk4, dopri5 };
enum class ReinjectionMode { none, independent, paired };
struct SimParams {
    // Simulation and integrator controls.
    std::size_t vortexCount = 100;
    double timeStep = 1.0e-3;
    double endTime = 1.0;
    double outputTime = 0.1;
    double coreRadius = 0.0;
    double absoluteTolerance = 1.0e-10;
    double relativeTolerance = 1.0e-8;
    double minimumTimeStep = 1.0e-12;
    double maximumTimeStep = 0.1;
    int numThreads = 0;
    IntegratorKind integrator = IntegratorKind::dopri5;

    // Geometry controls. Only the fields selected by boundaryCondition are active.
    std::string boundaryCondition = "infinite";
    double boxLengthX = 2.0;
    double boxLengthY = 2.0;
    double diskRadius = 1.0;
    int periodicImageLayers = 8;
    std::uint64_t randomSeed = 1234567;

    // Optional small-dipole removal and geometry-aware reinjection.
    bool dipoleRemoval = false;
    double dipoleRemovalDistance = 0.01;
    ReinjectionMode dipoleReinjection = ReinjectionMode::none;

    // Input, output, and restart paths.
    std::string initialConditionFile;
    std::string restartFile;
    std::string checkpointDirectory = "checkpoints";
    std::string outputFile = "vortices.csv";
    std::string diagnosticsFile = "diagnostics.csv";
    bool overwriteOutput = false;
    bool overwriteCheckpoints = false;
    void validate() const;
};
#endif
