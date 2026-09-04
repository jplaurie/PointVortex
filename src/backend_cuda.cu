#include "backend.h"
#include <algorithm>
#include <cmath>
#include <cuda_runtime.h>
#include <math_constants.h>
#include <stdexcept>
#include <string>

namespace {

void cudaCheck(cudaError_t status, const char *operation) {
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
}

enum class Geometry : int { infinite, periodic, disk };

__global__ void velocityKernel(const double *x, const double *y, const double *gamma, double *u,
                               double *v, std::size_t count, std::size_t begin, std::size_t end,
                               Geometry geometry, double first, double second, int imageLayers,
                               int *singular) {
    const std::size_t target = begin + blockIdx.x * blockDim.x + threadIdx.x;
    if (target >= end)
        return;
    constexpr double inverseTwoPi = 0.15915494309189533576888376337251;
    double velocityX = 0.0, velocityY = 0.0;

    for (std::size_t source = 0; source < count; ++source) {
        if (geometry == Geometry::infinite) {
            if (source == target)
                continue;
            const double dx = x[target] - x[source];
            const double dy = y[target] - y[source];
            const double denominator = dx * dx + dy * dy + first;
            if (denominator == 0.0) {
                atomicExch(singular, 1);
                continue;
            }
            const double coefficient = inverseTwoPi * gamma[source] / denominator;
            velocityX -= coefficient * dy;
            velocityY += coefficient * dx;
        } else if (geometry == Geometry::periodic) {
            constexpr double twoPi = 6.283185307179586476925286766559;
            const double waveNumber = twoPi / first;
            const double scale = 0.5 / first;
            const double dx = waveNumber * remainder(x[target] - x[source], first);
            const double dy = waveNumber * remainder(y[target] - y[source], second);
            const double sineX = sin(dx), sineY = sin(dy);
            const double cosineX = cos(dx), cosineY = cos(dy);
            for (int image = -imageLayers; image <= imageLayers; ++image) {
                if (source == target && image == 0)
                    continue;
                const double shiftedX = dx - twoPi * image;
                const double shiftedY = dy - twoPi * image;
                const double denominatorU =
                    fabs(shiftedX) > 40.0 ? CUDART_INF : cosh(shiftedX) - cosineY;
                const double denominatorV =
                    fabs(shiftedY) > 40.0 ? CUDART_INF : cosh(shiftedY) - cosineX;
                if (denominatorU == 0.0 || denominatorV == 0.0) {
                    atomicExch(singular, 1);
                    continue;
                }
                velocityX -= scale * gamma[source] * sineY / denominatorU;
                velocityY += scale * gamma[source] * sineX / denominatorV;
            }
        } else {
            if (source != target) {
                const double dx = x[target] - x[source];
                const double dy = y[target] - y[source];
                const double denominator = dx * dx + dy * dy;
                if (denominator == 0.0) {
                    atomicExch(singular, 1);
                } else {
                    const double coefficient = inverseTwoPi * gamma[source] / denominator;
                    velocityX -= coefficient * dy;
                    velocityY += coefficient * dx;
                }
            }
            const double sourceRadiusSquared = x[source] * x[source] + y[source] * y[source];
            if (sourceRadiusSquared != 0.0) {
                const double imageScale = first / sourceRadiusSquared;
                const double dx = x[target] - imageScale * x[source];
                const double dy = y[target] - imageScale * y[source];
                const double denominator = dx * dx + dy * dy;
                const double coefficient = -inverseTwoPi * gamma[source] / denominator;
                velocityX -= coefficient * dy;
                velocityY += coefficient * dx;
            }
        }
    }
    u[target] = velocityX;
    v[target] = velocityY;
}

class CudaKernel final : public VelocityKernel {
  public:
    explicit CudaKernel(const SimParams &params)
        : params_(params), cpu_(makeReferenceKernel(params)) {
        if (params.boundaryCondition == "infinite")
            geometry_ = Geometry::infinite;
        else if (params.boundaryCondition == "periodic")
            geometry_ = Geometry::periodic;
        else
            geometry_ = Geometry::disk;
    }
    ~CudaKernel() override { release(); }

    void evaluateRange(const std::vector<double> &x, const std::vector<double> &y,
                       const std::vector<double> &gamma, VelocityField &velocity, std::size_t begin,
                       std::size_t end) const override {
        const std::size_t count = x.size();
        if (y.size() != count || gamma.size() != count || begin > end || end > count)
            throw std::invalid_argument("invalid CUDA vortex arrays or target range");
        if (geometry_ == Geometry::periodic) {
            double total = 0.0, magnitude = 0.0;
            for (double value : gamma) {
                total += value;
                magnitude += std::abs(value);
            }
            if (std::abs(total) > 1e-12 * std::max(1.0, magnitude))
                throw std::invalid_argument("periodic box requires zero total circulation");
        }
        if (geometry_ == Geometry::disk)
            for (std::size_t i = 0; i < count; ++i)
                if (x[i] * x[i] + y[i] * y[i] >= params_.diskRadius * params_.diskRadius)
                    throw std::invalid_argument("vortex lies on or outside the disk");

        velocity.resize(count);
        if (begin == end)
            return;
        const std::size_t bytes = count * sizeof(double);
        ensureCapacity(count);
        cudaCheck(cudaMemcpy(deviceX_, x.data(), bytes, cudaMemcpyHostToDevice), "copy x to GPU");
        cudaCheck(cudaMemcpy(deviceY_, y.data(), bytes, cudaMemcpyHostToDevice), "copy y to GPU");
        cudaCheck(cudaMemcpy(deviceGamma_, gamma.data(), bytes, cudaMemcpyHostToDevice),
                  "copy circulation to GPU");
        cudaCheck(cudaMemset(deviceSingular_, 0, sizeof(int)), "clear GPU error flag");
        const int threads = 256;
        const int blocks = static_cast<int>((end - begin + threads - 1) / threads);
        const double first =
            geometry_ == Geometry::infinite
                ? params_.coreRadius * params_.coreRadius
                : (geometry_ == Geometry::periodic ? params_.boxLengthX
                                                   : params_.diskRadius * params_.diskRadius);
        velocityKernel<<<blocks, threads>>>(deviceX_, deviceY_, deviceGamma_, deviceU_, deviceV_,
                                            count, begin, end, geometry_, first, params_.boxLengthY,
                                            params_.periodicImageLayers, deviceSingular_);
        cudaCheck(cudaGetLastError(), "launch velocity kernel");
        cudaCheck(cudaMemcpy(velocity.x.data() + begin, deviceU_ + begin,
                             (end - begin) * sizeof(double), cudaMemcpyDeviceToHost),
                  "copy u from GPU");
        cudaCheck(cudaMemcpy(velocity.y.data() + begin, deviceV_ + begin,
                             (end - begin) * sizeof(double), cudaMemcpyDeviceToHost),
                  "copy v from GPU");
        int singular = 0;
        cudaCheck(cudaMemcpy(&singular, deviceSingular_, sizeof(int), cudaMemcpyDeviceToHost),
                  "copy GPU error flag");
        if (singular)
            throw std::runtime_error("coincident vortices in CUDA velocity kernel");
    }

    double hamiltonian(const VortexSystem &state) const override {
        return cpu_->hamiltonian(state);
    }

  private:
    void ensureCapacity(std::size_t count) const {
        if (count <= capacity_)
            return;
        release();
        const std::size_t bytes = count * sizeof(double);
        try {
            cudaCheck(cudaMalloc(&deviceX_, bytes), "cudaMalloc(x)");
            cudaCheck(cudaMalloc(&deviceY_, bytes), "cudaMalloc(y)");
            cudaCheck(cudaMalloc(&deviceGamma_, bytes), "cudaMalloc(circulation)");
            cudaCheck(cudaMalloc(&deviceU_, bytes), "cudaMalloc(u)");
            cudaCheck(cudaMalloc(&deviceV_, bytes), "cudaMalloc(v)");
            cudaCheck(cudaMalloc(&deviceSingular_, sizeof(int)), "cudaMalloc(error flag)");
            capacity_ = count;
        } catch (...) {
            release();
            throw;
        }
    }
    void release() const noexcept {
        cudaFree(deviceX_);
        cudaFree(deviceY_);
        cudaFree(deviceGamma_);
        cudaFree(deviceU_);
        cudaFree(deviceV_);
        cudaFree(deviceSingular_);
        deviceX_ = deviceY_ = deviceGamma_ = deviceU_ = deviceV_ = nullptr;
        deviceSingular_ = nullptr;
        capacity_ = 0;
    }
    SimParams params_;
    Geometry geometry_ = Geometry::infinite;
    std::unique_ptr<VelocityKernel> cpu_;
    mutable double *deviceX_ = nullptr, *deviceY_ = nullptr, *deviceGamma_ = nullptr;
    mutable double *deviceU_ = nullptr, *deviceV_ = nullptr;
    mutable int *deviceSingular_ = nullptr;
    mutable std::size_t capacity_ = 0;
};
} // namespace

void backendInitialize(int &, char **&) { cudaCheck(cudaFree(nullptr), "initialize CUDA"); }
void backendFinalize() {}
void backendAbort(int) {}
bool backendIsRoot() { return true; }
const char *backendName() { return "CUDA"; }
std::unique_ptr<VelocityKernel> makeBackendKernel(const SimParams &params) {
    return std::make_unique<CudaKernel>(params);
}
