#include "backend.h"

void backendInitialize(int &, char **&) {}
void backendFinalize() {}
void backendAbort(int) {}
bool backendIsRoot() { return true; }
const char *backendName() { return "CPU/OpenMP"; }
std::unique_ptr<VelocityKernel> makeBackendKernel(const SimParams &params) {
    return makeReferenceKernel(params);
}
