#include "backend.h"
#include <mpi.h>
#include <stdexcept>

namespace {
int rank = 0;
int ranks = 1;

class MpiKernel final : public VelocityKernel {
  public:
    explicit MpiKernel(std::unique_ptr<VelocityKernel> local) : local_(std::move(local)) {}

    void evaluateRange(const std::vector<double> &x, const std::vector<double> &y,
                       const std::vector<double> &circulation, VelocityField &velocity,
                       std::size_t begin, std::size_t end) const override {
        if (begin > end || end > x.size())
            throw std::out_of_range("invalid target-vortex range");
        const std::size_t targets = end - begin;
        const std::size_t localBegin =
            begin + targets * static_cast<std::size_t>(rank) / static_cast<std::size_t>(ranks);
        const std::size_t localEnd =
            begin + targets * static_cast<std::size_t>(rank + 1) / static_cast<std::size_t>(ranks);
        VelocityField local(x.size());
        int localFailure = 0;
        try {
            local_->evaluateRange(x, y, circulation, local, localBegin, localEnd);
        } catch (const std::exception &) {
            localFailure = 1;
        }
        int anyFailure = 0;
        MPI_Allreduce(&localFailure, &anyFailure, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (anyFailure)
            throw std::runtime_error("velocity evaluation failed on one or more MPI ranks");
        velocity.resize(x.size());
        MPI_Allreduce(local.x.data(), velocity.x.data(), static_cast<int>(x.size()), MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(local.y.data(), velocity.y.data(), static_cast<int>(x.size()), MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
    }

    double hamiltonian(const VortexSystem &state) const override {
        return local_->hamiltonian(state);
    }

  private:
    std::unique_ptr<VelocityKernel> local_;
};
} // namespace

void backendInitialize(int &argc, char **&argv) {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
        MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
}
void backendFinalize() {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();
}
void backendAbort(int exitCode) { MPI_Abort(MPI_COMM_WORLD, exitCode); }
bool backendIsRoot() { return rank == 0; }
const char *backendName() { return "MPI"; }
std::unique_ptr<VelocityKernel> makeBackendKernel(const SimParams &params) {
    return std::make_unique<MpiKernel>(makeReferenceKernel(params));
}
