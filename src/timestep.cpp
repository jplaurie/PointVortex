#include "timestep.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
RungeKuttaIntegrator::RungeKuttaIntegrator(std::size_t n) : initial_(n), temporary_(n) {
    for (auto &stage : stages_)
        stage.resize(n);
}
void RungeKuttaIntegrator::ensureSize(const VortexSystem &state) {
    state.validate();
    if (temporary_.size() != state.size()) {
        initial_.resize(state.size());
        temporary_.resize(state.size());
        for (auto &stage : stages_)
            stage.resize(state.size());
        fsalValid_ = false;
    }
}
void RungeKuttaIntegrator::makeStage(double dt, std::size_t stage,
                                     const std::array<double, 7> &coefficients) {
    for (std::size_t i = 0; i < initial_.size(); ++i) {
        double dx = 0.0, dy = 0.0;
        for (std::size_t j = 0; j < stage; ++j) {
            dx += coefficients[j] * stages_[j].x[i];
            dy += coefficients[j] * stages_[j].y[i];
        }
        temporary_.x[i] = initial_.x[i] + dt * dx;
        temporary_.y[i] = initial_.y[i] + dt * dy;
    }
}
void RungeKuttaIntegrator::rk4Step(VortexSystem &state, double dt, const VelocityKernel &kernel) {
    ensureSize(state);
    initial_.x = state.x;
    initial_.y = state.y;
    fsalValid_ = false;
    kernel.evaluate(initial_.x, initial_.y, state.circulation, stages_[0]);
    const std::array<double, 7> a2 = {0.5};
    makeStage(dt, 1, a2);
    kernel.evaluate(temporary_.x, temporary_.y, state.circulation, stages_[1]);
    const std::array<double, 7> a3 = {0.0, 0.5};
    makeStage(dt, 2, a3);
    kernel.evaluate(temporary_.x, temporary_.y, state.circulation, stages_[2]);
    const std::array<double, 7> a4 = {0.0, 0.0, 1.0};
    makeStage(dt, 3, a4);
    kernel.evaluate(temporary_.x, temporary_.y, state.circulation, stages_[3]);
    for (std::size_t i = 0; i < state.size(); ++i) {
        state.x[i] = initial_.x[i] + dt *
                                         (stages_[0].x[i] + 2.0 * stages_[1].x[i] +
                                          2.0 * stages_[2].x[i] + stages_[3].x[i]) /
                                         6.0;
        state.y[i] = initial_.y[i] + dt *
                                         (stages_[0].y[i] + 2.0 * stages_[1].y[i] +
                                          2.0 * stages_[2].y[i] + stages_[3].y[i]) /
                                         6.0;
    }
}
StepResult RungeKuttaIntegrator::dopri5Step(VortexSystem &state, double dt,
                                            const VelocityKernel &kernel, const SimParams &p) {
    ensureSize(state);
    initial_.x = state.x;
    initial_.y = state.y;
    static constexpr std::array<std::array<double, 7>, 7> a = {
        {{},
         {1.0 / 5.0},
         {3.0 / 40.0, 9.0 / 40.0},
         {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0},
         {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0},
         {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0},
         {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0}}};
    static constexpr double b5[7] = {
        35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0};
    static constexpr double b4[7] = {
        5179.0 / 57600.0, 0.0,       7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0,
        187.0 / 2100.0,   1.0 / 40.0};
    unsigned rejected = 0;
    // Dormand--Prince is FSAL: k7 from an accepted step is k1 of the next step.
    if (!fsalValid_)
        kernel.evaluate(initial_.x, initial_.y, state.circulation, stages_[0]);
    for (;;) {
        for (std::size_t s = 1; s < 7; ++s) {
            makeStage(dt, s, a[s]);
            kernel.evaluate(temporary_.x, temporary_.y, state.circulation, stages_[s]);
        }
        double error = 0.0;
        for (std::size_t i = 0; i < state.size(); ++i) {
            double x5 = initial_.x[i], y5 = initial_.y[i];
            double x4 = initial_.x[i], y4 = initial_.y[i];
            for (std::size_t s = 0; s < 7; ++s) {
                x5 += dt * b5[s] * stages_[s].x[i];
                y5 += dt * b5[s] * stages_[s].y[i];
                x4 += dt * b4[s] * stages_[s].x[i];
                y4 += dt * b4[s] * stages_[s].y[i];
            }
            temporary_.x[i] = x5;
            temporary_.y[i] = y5;
            const double xScale =
                p.absoluteTolerance +
                p.relativeTolerance * std::max(std::abs(initial_.x[i]), std::abs(x5));
            const double yScale =
                p.absoluteTolerance +
                p.relativeTolerance * std::max(std::abs(initial_.y[i]), std::abs(y5));
            error = std::max(error, std::abs(x5 - x4) / xScale);
            error = std::max(error, std::abs(y5 - y4) / yScale);
        }
        const double factor =
            error == 0.0 ? 5.0 : std::clamp(0.9 * std::pow(error, -0.2), 0.2, 5.0);
        const double suggested = std::clamp(dt * factor, p.minimumTimeStep, p.maximumTimeStep);
        if (error <= 1.0) {
            state.x.swap(temporary_.x);
            state.y.swap(temporary_.y);
            stages_[0].x.swap(stages_[6].x);
            stages_[0].y.swap(stages_[6].y);
            fsalValid_ = true;
            return {dt, suggested, error, rejected};
        }
        if (dt <= p.minimumTimeStep || ++rejected > 32)
            throw std::runtime_error("adaptive integrator could not satisfy tolerance");
        dt = std::max(p.minimumTimeStep, std::min(suggested, dt * 0.9));
    }
}
