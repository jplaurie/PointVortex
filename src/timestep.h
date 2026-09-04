#ifndef POINT_VORTEX_TIMESTEP_H
#define POINT_VORTEX_TIMESTEP_H
#include "compute.h"
#include "params.h"
#include <array>
struct StepResult {
    // The caller advances by acceptedTimeStep and proposes suggestedTimeStep next.
    double acceptedTimeStep;
    double suggestedTimeStep;
    double normalizedError;
    unsigned rejectedSteps;
};
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
    void makeStage(double dt, std::size_t stage, const std::array<double, 7> &coefficients);
};
#endif
