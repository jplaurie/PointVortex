#include "backend.h"

void backendInitialize(int &, char **&) {}
void backendFinalize() {}
void backendAbort(int) {}
bool backendIsRoot() { return true; }
const char *backendName() {
#ifdef _OPENMP
    return "CPU/OpenMP";
#else
    return "CPU";
#endif
}
std::unique_ptr<VelocityKernel> makeBackendKernel(const SimParams &params) {
    return makeReferenceKernel(params);
}
