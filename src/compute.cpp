#include "compute.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace {

Invariants computeMoments(const VortexSystem &vortices) {
    vortices.validate();
    Invariants result;
    for (std::size_t i = 0; i < vortices.size(); ++i) {
        const double gamma = vortices.circulation[i];
        result.circulation += gamma;
        result.linearImpulseX += gamma * vortices.x[i];
        result.linearImpulseY += gamma * vortices.y[i];
        result.angularImpulse +=
            gamma * (vortices.x[i] * vortices.x[i] + vortices.y[i] * vortices.y[i]);
    }
    for (double value :
         {result.circulation, result.linearImpulseX, result.linearImpulseY, result.angularImpulse})
        if (!std::isfinite(value))
            throw std::runtime_error("non-finite vortex moment; check input scales");
    return result;
}

} // namespace

InfinitePlaneKernel::InfinitePlaneKernel(double coreRadius)
    : coreRadiusSquared_(coreRadius * coreRadius) {
    if (!(coreRadius >= 0.0) || !std::isfinite(coreRadiusSquared_) ||
        (coreRadius > 0.0 && coreRadiusSquared_ == 0.0))
        throw std::invalid_argument(
            "core radius is negative or outside the supported numeric range");
}
void VelocityKernel::evaluate(const VortexSystem &vortices, VelocityField &velocity) const {
    evaluate(vortices.x, vortices.y, vortices.circulation, velocity);
}
void VelocityKernel::evaluate(const std::vector<double> &x, const std::vector<double> &y,
                              const std::vector<double> &circulation,
                              VelocityField &velocity) const {
    evaluateRange(x, y, circulation, velocity, 0, x.size());
    if (velocity.x.size() != x.size() || velocity.y.size() != x.size())
        throw std::runtime_error("kernel returned invalid velocity arrays");
    for (std::size_t i = 0; i < x.size(); ++i)
        if (!std::isfinite(velocity.x[i]) || !std::isfinite(velocity.y[i]))
            throw std::runtime_error("non-finite velocity; check scales and close encounters");
}
void VelocityKernel::uploadDeviceState(const VortexSystem &) const {
    throw std::logic_error("selected backend does not support device-resident integration");
}
void VelocityKernel::downloadDeviceState(VortexSystem &) const {
    throw std::logic_error("selected backend does not support device-resident integration");
}
void VelocityKernel::evaluateDeviceState(VelocityField &) const {
    throw std::logic_error("selected backend does not support device-resident integration");
}
void VelocityKernel::deviceRk4Step(double) const {
    throw std::logic_error("selected backend does not support device-resident integration");
}
DeviceStepResult VelocityKernel::deviceDopri5Step(double, double, double, double, double) const {
    throw std::logic_error("selected backend does not support device-resident integration");
}
double InfinitePlaneKernel::hamiltonian(const VortexSystem &vortices) const {
    return computeInvariants(vortices, std::sqrt(coreRadiusSquared_)).hamiltonian;
}
void InfinitePlaneKernel::evaluateRange(const std::vector<double> &x, const std::vector<double> &y,
                                        const std::vector<double> &circulation,
                                        VelocityField &velocity, std::size_t begin,
                                        std::size_t end) const {
    const std::size_t count = x.size();
    validateVortexArrays(x, y, circulation);
    if (begin > end || end > count)
        throw std::out_of_range("invalid target-vortex range");
    velocity.resize(count);
    constexpr double inverseTwoPi = 0.5 / std::numbers::pi;
    int singularPair = 0;
    const double *const px = x.data();
    const double *const py = y.data();
    const double *const pg = circulation.data();
    double *const pu = velocity.x.data();
    double *const pv = velocity.y.data();
    // Each target owns its accumulators, avoiding atomics and reduction arrays.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (count >= 256) reduction(| : singularPair)
#endif
    for (std::ptrdiff_t target = static_cast<std::ptrdiff_t>(begin);
         target < static_cast<std::ptrdiff_t>(end); ++target) {
        double u = 0.0, v = 0.0;
#ifdef _OPENMP
#pragma omp simd reduction(+ : u, v) reduction(| : singularPair)
#endif
        for (std::ptrdiff_t source = 0; source < static_cast<std::ptrdiff_t>(count); ++source) {
            if (source == target)
                continue;
            const double dx = px[target] - px[source];
            const double dy = py[target] - py[source];
            const double denominator = dx * dx + dy * dy + coreRadiusSquared_;
            if (denominator == 0.0) {
                singularPair = 1;
                continue;
            }
            const double coefficient = inverseTwoPi * pg[source] / denominator;
            u -= coefficient * dy;
            v += coefficient * dx;
        }
        pu[target] = u;
        pv[target] = v;
    }
    if (singularPair != 0)
        throw std::runtime_error("coincident vortices in singular point-vortex kernel");
}
Invariants computeInvariants(const VortexSystem &vortices, double coreRadius) {
    if (!(coreRadius >= 0.0) || !std::isfinite(coreRadius) ||
        !std::isfinite(coreRadius * coreRadius) ||
        (coreRadius > 0.0 && coreRadius * coreRadius == 0.0))
        throw std::invalid_argument(
            "core radius is negative or outside the supported numeric range");
    Invariants result = computeMoments(vortices);
    const double epsilonSquared = coreRadius * coreRadius;
    for (std::size_t i = 0; i < vortices.size(); ++i) {
        for (std::size_t j = i + 1; j < vortices.size(); ++j) {
            const double dx = vortices.x[i] - vortices.x[j];
            const double dy = vortices.y[i] - vortices.y[j];
            const double r2 = dx * dx + dy * dy + epsilonSquared;
            if (r2 == 0.0)
                throw std::runtime_error("coincident vortices");
            result.hamiltonian -= vortices.circulation[i] * vortices.circulation[j] * std::log(r2) /
                                  (4.0 * std::numbers::pi);
        }
    }
    if (!std::isfinite(result.hamiltonian))
        throw std::runtime_error("non-finite Hamiltonian; check input scales");
    return result;
}
Invariants computeInvariants(const VortexSystem &vortices, const VelocityKernel &kernel) {
    Invariants result = computeMoments(vortices);
    for (double value :
         {result.circulation, result.linearImpulseX, result.linearImpulseY, result.angularImpulse})
        if (!std::isfinite(value))
            throw std::runtime_error("non-finite vortex moment; check input scales");
    result.hamiltonian = kernel.hamiltonian(vortices);
    if (!std::isfinite(result.hamiltonian))
        throw std::runtime_error("non-finite Hamiltonian; check input scales");
    return result;
}

PeriodicBoxKernel::PeriodicBoxKernel(double lengthX, double lengthY, int imageLayers)
    : lengthX_(lengthX), lengthY_(lengthY), imageLayers_(imageLayers) {
    if (!(lengthX > 0.0) || !(lengthY > 0.0) || !std::isfinite(lengthX) ||
        !std::isfinite(lengthY) || !std::isfinite(2.0 * std::numbers::pi / lengthX) ||
        imageLayers < 0 || imageLayers > 64)
        throw std::invalid_argument("invalid periodic-box parameters");
    if (std::abs(lengthX - lengthY) > 1e-13 * std::max(lengthX, lengthY))
        throw std::invalid_argument("Weiss-McWilliams kernel requires a square periodic box");
}
double PeriodicBoxKernel::hamiltonian(const VortexSystem &vortices) const {
    vortices.validate();
    double value = 0.0;
    const double waveNumber = 2.0 * std::numbers::pi / lengthX_;
    const auto logCosh = [](double argument) {
        const double magnitude = std::abs(argument);
        return magnitude + std::log1p(std::exp(-2.0 * magnitude)) - std::log(2.0);
    };
    for (std::size_t i = 0; i < vortices.size(); ++i) {
        for (std::size_t j = i + 1; j < vortices.size(); ++j) {
            const double dx = waveNumber * std::remainder(vortices.x[i] - vortices.x[j], lengthX_);
            const double dy = waveNumber * std::remainder(vortices.y[i] - vortices.y[j], lengthY_);
            // The quadratic term fixes the conditionally convergent lattice-sum limit.
            double pairEnergy = -dx * dx / (2.0 * std::numbers::pi);
            for (int image = -imageLayers_; image <= imageLayers_; ++image) {
                const double shifted = dx - 2.0 * std::numbers::pi * image;
                const double sinhHalf = std::sinh(0.5 * shifted);
                const double sinHalf = std::sin(0.5 * dy);
                const double denominator = 2.0 * (sinhHalf * sinhHalf + sinHalf * sinHalf);
                if (denominator == 0.0)
                    throw std::runtime_error("coincident periodic vortices");
                pairEnergy += std::log(denominator) - logCosh(2.0 * std::numbers::pi * image);
            }
            value -= vortices.circulation[i] * vortices.circulation[j] * pairEnergy /
                     (4.0 * std::numbers::pi);
        }
    }
    return value;
}
void PeriodicBoxKernel::evaluateRange(const std::vector<double> &x, const std::vector<double> &y,
                                      const std::vector<double> &circulation,
                                      VelocityField &velocity, std::size_t begin,
                                      std::size_t end) const {
    const std::size_t count = x.size();
    validateVortexArrays(x, y, circulation);
    if (begin > end || end > count)
        throw std::out_of_range("invalid target-vortex range");
    double total = 0.0, absoluteTotal = 0.0;
    for (double gamma : circulation) {
        total += gamma;
        absoluteTotal += std::abs(gamma);
    }
    if (std::abs(total) > 1e-12 * std::max(1.0, absoluteTotal))
        throw std::invalid_argument("periodic box requires zero total circulation");
    velocity.resize(count);
    const double waveNumber = 2.0 * std::numbers::pi / lengthX_;
    const double scale = 0.5 / lengthX_;
    int singularPair = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (count >= 256) reduction(| : singularPair)
#endif
    for (std::ptrdiff_t target = static_cast<std::ptrdiff_t>(begin);
         target < static_cast<std::ptrdiff_t>(end); ++target) {
        double u = 0.0, v = 0.0;
        for (std::ptrdiff_t source = 0; source < static_cast<std::ptrdiff_t>(count); ++source) {
            const double scaledDx = waveNumber * std::remainder(x[target] - x[source], lengthX_);
            const double scaledDy = waveNumber * std::remainder(y[target] - y[source], lengthY_);
            const double sineX = std::sin(scaledDx), sineY = std::sin(scaledDy);
            const double sinHalfX = std::sin(0.5 * scaledDx);
            const double sinHalfY = std::sin(0.5 * scaledDy);
            // Complementary sum orientations minimize finite-truncation drift.
            for (int image = -imageLayers_; image <= imageLayers_; ++image) {
                if (source == target && image == 0)
                    continue;
                const double shiftedX = scaledDx - 2.0 * std::numbers::pi * image;
                const double shiftedY = scaledDy - 2.0 * std::numbers::pi * image;
                // cosh(a)-cos(b) loses all significant digits for close pairs.
                const auto denominator = [](double a, double sinHalfB) {
                    if (std::abs(a) > 40.0)
                        return std::numeric_limits<double>::infinity();
                    const double sinhHalfA = std::sinh(0.5 * a);
                    return 2.0 * (sinhHalfA * sinhHalfA + sinHalfB * sinHalfB);
                };
                const double denominatorU = denominator(shiftedX, sinHalfY);
                const double denominatorV = denominator(shiftedY, sinHalfX);
                if (denominatorU == 0.0 || denominatorV == 0.0) {
                    singularPair = 1;
                    continue;
                }
                u -= scale * circulation[source] * sineY / denominatorU;
                v += scale * circulation[source] * sineX / denominatorV;
            }
        }
        velocity.x[target] = u;
        velocity.y[target] = v;
    }
    if (singularPair != 0)
        throw std::runtime_error("coincident periodic vortices");
}

DiskKernel::DiskKernel(double radius) : radius_(radius), radiusSquared_(radius * radius) {
    if (!(radius > 0.0) || !std::isfinite(radiusSquared_) || radiusSquared_ == 0.0)
        throw std::invalid_argument("invalid disk parameters");
}
double DiskKernel::hamiltonian(const VortexSystem &vortices) const {
    vortices.validate();
    double value = 0.0;
    for (std::size_t i = 0; i < vortices.size(); ++i) {
        const double ri2 = vortices.x[i] * vortices.x[i] + vortices.y[i] * vortices.y[i];
        if (ri2 >= radiusSquared_)
            throw std::invalid_argument("vortex outside disk");
        value += vortices.circulation[i] * vortices.circulation[i] *
                 std::log((radiusSquared_ - ri2) / radius_) / (4.0 * std::numbers::pi);
        for (std::size_t j = i + 1; j < vortices.size(); ++j) {
            const double dx = vortices.x[i] - vortices.x[j], dy = vortices.y[i] - vortices.y[j];
            const double numerator = (dx * dx + dy * dy) * radiusSquared_;
            const double dot = vortices.x[i] * vortices.x[j] + vortices.y[i] * vortices.y[j];
            const double cross = vortices.y[i] * vortices.x[j] - vortices.x[i] * vortices.y[j];
            const double denominator =
                (radiusSquared_ - dot) * (radiusSquared_ - dot) + cross * cross;
            if (numerator == 0.0)
                throw std::runtime_error("coincident vortices in disk");
            value -= vortices.circulation[i] * vortices.circulation[j] *
                     std::log(numerator / denominator) / (4.0 * std::numbers::pi);
        }
    }
    return value;
}
void DiskKernel::evaluateRange(const std::vector<double> &x, const std::vector<double> &y,
                               const std::vector<double> &circulation, VelocityField &velocity,
                               std::size_t begin, std::size_t end) const {
    const std::size_t count = x.size();
    validateVortexArrays(x, y, circulation);
    if (begin > end || end > count)
        throw std::out_of_range("invalid target-vortex range");
    for (std::size_t i = 0; i < count; ++i)
        if (x[i] * x[i] + y[i] * y[i] >= radiusSquared_)
            throw std::invalid_argument("vortex lies on or outside the disk");
    velocity.resize(count);
    constexpr double inverseTwoPi = 0.5 / std::numbers::pi;
    int singularPair = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (count >= 256) reduction(| : singularPair)
#endif
    for (std::ptrdiff_t target = static_cast<std::ptrdiff_t>(begin);
         target < static_cast<std::ptrdiff_t>(end); ++target) {
        double u = 0.0, v = 0.0;
        for (std::ptrdiff_t source = 0; source < static_cast<std::ptrdiff_t>(count); ++source) {
            if (source != target) {
                const double dx = x[target] - x[source];
                const double dy = y[target] - y[source];
                const double denominator = dx * dx + dy * dy;
                if (denominator == 0.0) {
                    singularPair = 1;
                    continue;
                }
                const double coefficient = inverseTwoPi * circulation[source] / denominator;
                u -= coefficient * dy;
                v += coefficient * dx;
            }
            // Evaluate the image contribution without constructing an inverse point.
            // This stays finite when a source is at or extremely close to the center.
            const double tx = x[target] / radius_, ty = y[target] / radius_;
            const double sx = x[source] / radius_, sy = y[source] / radius_;
            const double a = 1.0 - (tx * sx + ty * sy);
            const double b = ty * sx - tx * sy;
            const double denominator = a * a + b * b;
            const double imageX = -a * sx - b * sy;
            const double imageY = -a * sy + b * sx;
            const double coefficient = -inverseTwoPi * circulation[source] / radius_ / denominator;
            u -= coefficient * imageY;
            v += coefficient * imageX;
        }
        velocity.x[target] = u;
        velocity.y[target] = v;
    }
    if (singularPair != 0)
        throw std::runtime_error("coincident vortices in disk");
}
