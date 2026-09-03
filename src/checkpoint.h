#ifndef POINT_VORTEX_CHECKPOINT_H
#define POINT_VORTEX_CHECKPOINT_H
#include "compute.h"
#include "dipole.h"
#include "params.h"
#include "vortex.h"
#include <cstddef>
#include <filesystem>
struct Checkpoint {
    // Restart state includes output/integrator progress as well as vortex data.
    VortexSystem vortices;
    Invariants initialInvariants;
    Invariants segmentInvariants;
    bool hasSegmentInvariants = false;
    double time = 0.0;
    double suggestedTimeStep = 0.0;
    double nextOutputTime = 0.0;
    std::size_t acceptedSteps = 0;
    std::size_t outputIndex = 0;
    double coreRadius = 0.0;
    IntegratorKind integrator = IntegratorKind::dopri5;
    std::string boundaryCondition = "infinite";
    double geometryLengthX = 0.0;
    double geometryLengthY = 0.0;
    int periodicImageLayers = 0;
    bool dipoleRemoval = false;
    double dipoleRemovalDistance = 0.0;
    ReinjectionMode dipoleReinjection = ReinjectionMode::none;
    DipoleEventState dipoleState;
};
Checkpoint loadCheckpoint(const std::filesystem::path &filename);
void writeCheckpoint(const std::filesystem::path &directory, const VortexSystem &vortices,
                     const Invariants &initialInvariants, double time, double suggestedTimeStep,
                     double nextOutputTime, std::size_t acceptedSteps, std::size_t outputIndex,
                     double coreRadius, IntegratorKind integrator,
                     const std::string &boundaryCondition, double geometryLengthX,
                     double geometryLengthY, int periodicImageLayers, bool dipoleRemoval,
                     double dipoleRemovalDistance, ReinjectionMode dipoleReinjection,
                     const DipoleEventState &dipoleState, const Invariants &segmentInvariants,
                     bool overwrite = false);
std::filesystem::path checkpointPath(const std::filesystem::path &directory,
                                     std::size_t outputIndex);
#endif
