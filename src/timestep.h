#ifndef POINT_VORTEX_TIMESTEP_H
#define POINT_VORTEX_TIMESTEP_H
#include "compute.h"
#include "params.h"
#include <array>

namespace integrator_detail {
inline constexpr std::array<std::array<double, 7>, 7> dopriCoefficients = {
    {{},
     {1.0 / 5.0},
     {3.0 / 40.0, 9.0 / 40.0},
     {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0},
     {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0},
     {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0},
     {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0}}};
} // namespace integrator_detail

class RungeKuttaIntegrator {
  public:
    explicit RungeKuttaIntegrator(std::size_t vortexCount);
    void rk4Step(VortexSystem &state, double dt, const VelocityKernel &kernel);
    StepResult dopri5Step(VortexSystem &state, double dt, const VelocityKernel &kernel,
                          const SimParams &params);
    // Call after a discrete event changes positions or circulations between timesteps.
    void invalidateCachedDerivative() noexcept { fsalValid_ = false; }

  private:
    PositionField initial_;
    PositionField temporary_;
    std::array<VelocityField, 7> stages_;
    bool fsalValid_ = false; // stages_[0] is the derivative at the current state.
    void ensureSize(const VortexSystem &state);
    void makeStage(std::size_t stage, const std::array<double, 7> &coefficients, double dt);
};
#endif
