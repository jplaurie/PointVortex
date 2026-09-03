#ifndef POINT_VORTEX_COMPUTE_H
#define POINT_VORTEX_COMPUTE_H
#include "vortex.h"
#include <vector>
// Geometry-independent right-hand side used by both time integrators.
class VelocityKernel {
  public:
    virtual ~VelocityKernel() = default;
    void evaluate(const VortexSystem &, VelocityField &) const;
    virtual void evaluate(const std::vector<double> &x, const std::vector<double> &y,
                          const std::vector<double> &circulation,
                          VelocityField &velocity) const = 0;
    virtual double hamiltonian(const VortexSystem &) const = 0;
};
class InfinitePlaneKernel final : public VelocityKernel {
  public:
    using VelocityKernel::evaluate;
    explicit InfinitePlaneKernel(double coreRadius = 0.0);
    void evaluate(const std::vector<double> &, const std::vector<double> &,
                  const std::vector<double> &, VelocityField &) const override;
    double hamiltonian(const VortexSystem &) const override;

  private:
    double coreRadiusSquared_;
};
// Weiss--McWilliams rapidly convergent image sum for a square torus.
class PeriodicBoxKernel final : public VelocityKernel {
  public:
    using VelocityKernel::evaluate;
    PeriodicBoxKernel(double lengthX, double lengthY, int imageLayers = 8);
    void evaluate(const std::vector<double> &, const std::vector<double> &,
                  const std::vector<double> &, VelocityField &) const override;
    double hamiltonian(const VortexSystem &) const override;

  private:
    double lengthX_;
    double lengthY_;
    int imageLayers_;
};
// Impermeable circular wall represented by opposite-sign inverse-point images.
class DiskKernel final : public VelocityKernel {
  public:
    using VelocityKernel::evaluate;
    explicit DiskKernel(double radius);
    void evaluate(const std::vector<double> &, const std::vector<double> &,
                  const std::vector<double> &, VelocityField &) const override;
    double hamiltonian(const VortexSystem &) const override;

  private:
    double radius_;
    double radiusSquared_;
};
struct Invariants {
    double circulation = 0.0, linearImpulseX = 0.0, linearImpulseY = 0.0;
    double angularImpulse = 0.0, hamiltonian = 0.0;
};
Invariants computeInvariants(const VortexSystem &, double coreRadius = 0.0);
Invariants computeInvariants(const VortexSystem &, const VelocityKernel &kernel);
#endif
