#ifndef POINT_VORTEX_VORTEX_H
#define POINT_VORTEX_VORTEX_H
#include <cstddef>
#include <stdexcept>
#include <vector>
struct VortexSystem {
    // Structure-of-arrays storage keeps each hot numerical stream contiguous.
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> circulation;
    explicit VortexSystem(std::size_t n = 0) : x(n), y(n), circulation(n) {}
    [[nodiscard]] std::size_t size() const noexcept { return x.size(); }
    void resize(std::size_t n) {
        x.resize(n);
        y.resize(n);
        circulation.resize(n);
    }
    void validate() const {
        if (x.size() != y.size() || x.size() != circulation.size())
            throw std::invalid_argument("vortex arrays have different lengths");
    }
};
// Runge--Kutta stage positions do not need a copy of the constant circulations.
struct PositionField {
    std::vector<double> x;
    std::vector<double> y;
    explicit PositionField(std::size_t n = 0) : x(n), y(n) {}
    [[nodiscard]] std::size_t size() const noexcept { return x.size(); }
    void resize(std::size_t n) {
        x.resize(n);
        y.resize(n);
    }
};
struct VelocityField {
    std::vector<double> x;
    std::vector<double> y;
    explicit VelocityField(std::size_t n = 0) : x(n), y(n) {}
    void resize(std::size_t n) {
        x.resize(n);
        y.resize(n);
    }
};
#endif
