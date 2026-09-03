#ifndef POINT_VORTEX_READ_H
#define POINT_VORTEX_READ_H
#include "params.h"
#include "vortex.h"
#include <cstdint>
#include <string>
SimParams loadParams(const std::string &filename);
VortexSystem loadVortices(const std::string &filename);
void initializeVortices(VortexSystem &vortices, double radius = 1.0);
void initializePeriodicVortices(VortexSystem &vortices, double lengthX, double lengthY,
                                std::uint64_t seed);
#endif
