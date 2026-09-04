#ifndef POINT_VORTEX_BACKEND_H
#define POINT_VORTEX_BACKEND_H

#include "compute.h"
#include "params.h"
#include <memory>

void backendInitialize(int &argc, char **&argv);
void backendFinalize();
void backendAbort(int exitCode);
[[nodiscard]] bool backendIsRoot();
[[nodiscard]] const char *backendName();
std::unique_ptr<VelocityKernel> makeReferenceKernel(const SimParams &params);
std::unique_ptr<VelocityKernel> makeBackendKernel(const SimParams &params);

#endif
