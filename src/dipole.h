#ifndef POINT_VORTEX_DIPOLE_H
#define POINT_VORTEX_DIPOLE_H

#include "params.h"
#include "vortex.h"

#include <cstddef>
#include <cstdint>
#include <random>
#include <string>

struct DipoleEventState {
    std::size_t removedPairs = 0;
    std::size_t reinjectedPairs = 0;
    std::string randomEngineState;
};

class DipoleManager {
  public:
    explicit DipoleManager(const SimParams &params);
    DipoleManager(const SimParams &params, const DipoleEventState &state);

    // Remove all closest-first, disjoint opposite-sign pairs below the threshold.
    // Returns the number of pairs processed during this call.
    std::size_t process(VortexSystem &vortices);
    [[nodiscard]] DipoleEventState state() const;

  private:
    const SimParams &params_;
    std::mt19937_64 random_;
    std::size_t removedPairs_ = 0;
    std::size_t reinjectedPairs_ = 0;

    void injectPair(VortexSystem &vortices, double firstCirculation, double secondCirculation,
                    std::size_t population);
};

#endif
