#ifndef POINT_VORTEX_READ_H
#define POINT_VORTEX_READ_H
#include "params.h"
#include "vortex.h"
#include <cstdint>
#include <optional>
#include <string>
struct InitialConditionMetadata {
    std::optional<std::string> geometry;
    std::optional<double> boxLength;
    std::optional<double> diskRadius;
};
SimParams loadParams(const std::string &filename);
VortexSystem loadVortices(const std::string &filename);
InitialConditionMetadata readInitialConditionMetadata(const std::string &filename);
void initializeVortices(VortexSystem &vortices, double radius = 1.0);
void initializePeriodicVortices(VortexSystem &vortices, double lengthX, double lengthY,
                                std::uint64_t seed);
#endif
