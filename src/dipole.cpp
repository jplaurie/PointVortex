#include "dipole.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
struct Candidate {
    std::size_t first;
    std::size_t second;
    double distanceSquared;
};

double displacement(double difference, double length, bool periodic) {
    return periodic ? std::remainder(difference, length) : difference;
}
} // namespace

DipoleManager::DipoleManager(const SimParams &params)
    : params_(params), random_(params.randomSeed ^ 0xd1b54a32d192ed03ULL) {}

DipoleManager::DipoleManager(const SimParams &params, const DipoleEventState &state)
    : params_(params), random_(params.randomSeed ^ 0xd1b54a32d192ed03ULL),
      removedPairs_(state.removedPairs), reinjectedPairs_(state.reinjectedPairs) {
    if (!state.randomEngineState.empty()) {
        std::istringstream input(state.randomEngineState);
        if (!(input >> random_))
            throw std::runtime_error("invalid random-generator state in checkpoint");
    }
}

std::size_t DipoleManager::process(VortexSystem &vortices) {
    if (!params_.dipoleRemoval || vortices.size() < 2)
        return 0;

    const bool periodic = params_.boundaryCondition == "periodic";
    const double thresholdSquared = params_.dipoleRemovalDistance * params_.dipoleRemovalDistance;
    std::vector<Candidate> candidates;
    for (std::size_t i = 0; i < vortices.size(); ++i) {
        if (vortices.circulation[i] == 0.0)
            continue;
        for (std::size_t j = i + 1; j < vortices.size(); ++j) {
            if (vortices.circulation[j] == 0.0 ||
                std::signbit(vortices.circulation[i]) == std::signbit(vortices.circulation[j]))
                continue;
            const double dx =
                displacement(vortices.x[i] - vortices.x[j], params_.boxLengthX, periodic);
            const double dy =
                displacement(vortices.y[i] - vortices.y[j], params_.boxLengthY, periodic);
            const double distanceSquared = dx * dx + dy * dy;
            if (distanceSquared < thresholdSquared)
                candidates.push_back({i, j, distanceSquared});
        }
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate &left, const Candidate &right) {
                  if (left.distanceSquared != right.distanceSquared)
                      return left.distanceSquared < right.distanceSquared;
                  return std::pair{left.first, left.second} < std::pair{right.first, right.second};
              });

    std::vector<bool> selected(vortices.size(), false);
    std::vector<std::pair<double, double>> circulations;
    for (const Candidate &candidate : candidates) {
        if (selected[candidate.first] || selected[candidate.second])
            continue;
        selected[candidate.first] = selected[candidate.second] = true;
        circulations.emplace_back(vortices.circulation[candidate.first],
                                  vortices.circulation[candidate.second]);
    }
    if (circulations.empty())
        return 0;

    const std::size_t originalPopulation = vortices.size();
    VortexSystem survivors;
    survivors.x.reserve(originalPopulation - 2 * circulations.size());
    survivors.y.reserve(originalPopulation - 2 * circulations.size());
    survivors.circulation.reserve(originalPopulation - 2 * circulations.size());
    for (std::size_t i = 0; i < originalPopulation; ++i) {
        if (!selected[i]) {
            survivors.x.push_back(vortices.x[i]);
            survivors.y.push_back(vortices.y[i]);
            survivors.circulation.push_back(vortices.circulation[i]);
        }
    }
    vortices = std::move(survivors);
    removedPairs_ += circulations.size();

    if (params_.dipoleReinjection != ReinjectionMode::none) {
        for (const auto &[first, second] : circulations)
            injectPair(vortices, first, second, originalPopulation);
        reinjectedPairs_ += circulations.size();
    }
    return circulations.size();
}

void DipoleManager::injectPair(VortexSystem &vortices, double firstCirculation,
                               double secondCirculation, std::size_t population) {
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const bool periodic = params_.boundaryCondition == "periodic";
    const double area = periodic ? params_.boxLengthX * params_.boxLengthY
                                 : std::numbers::pi * params_.diskRadius * params_.diskRadius;
    const double spacing = std::sqrt(area / static_cast<double>(population));

    const auto randomPosition = [&]() {
        if (periodic)
            return std::pair{(unit(random_) - 0.5) * params_.boxLengthX,
                             (unit(random_) - 0.5) * params_.boxLengthY};
        const double radius = params_.diskRadius * std::sqrt(unit(random_));
        const double angle = 2.0 * std::numbers::pi * unit(random_);
        return std::pair{radius * std::cos(angle), radius * std::sin(angle)};
    };

    auto first = randomPosition();
    auto second = randomPosition();
    if (params_.dipoleReinjection == ReinjectionMode::paired) {
        constexpr std::size_t maximumAttempts = 100000;
        bool placed = false;
        for (std::size_t attempt = 0; attempt < maximumAttempts; ++attempt) {
            first = randomPosition();
            const double angle = 2.0 * std::numbers::pi * unit(random_);
            second = {first.first + spacing * std::cos(angle),
                      first.second + spacing * std::sin(angle)};
            if (periodic) {
                second.first = std::remainder(second.first, params_.boxLengthX);
                second.second = std::remainder(second.second, params_.boxLengthY);
                placed = true;
                break;
            }
            if (second.first * second.first + second.second * second.second <
                params_.diskRadius * params_.diskRadius) {
                placed = true;
                break;
            }
        }
        if (!placed)
            throw std::runtime_error("could not find a valid paired reinjection position");
    }
    vortices.x.push_back(first.first);
    vortices.y.push_back(first.second);
    vortices.circulation.push_back(firstCirculation);
    vortices.x.push_back(second.first);
    vortices.y.push_back(second.second);
    vortices.circulation.push_back(secondCirculation);
}

DipoleEventState DipoleManager::state() const {
    std::ostringstream output;
    output << random_;
    return {removedPairs_, reinjectedPairs_, output.str()};
}
