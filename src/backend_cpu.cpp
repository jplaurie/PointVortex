#include "backend.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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
std::string backendRuntimeDetails() {
#ifdef _OPENMP
    return "openmp_enabled true\nopenmp_threads " + std::to_string(omp_get_max_threads());
#else
    return "openmp_enabled false\nopenmp_threads 1";
#endif
}
std::unique_ptr<VelocityKernel> makeBackendKernel(const SimParams &params) {
    return makeReferenceKernel(params);
}
