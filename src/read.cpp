#include "read.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numbers>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>

namespace {
double parseDouble(const std::string &value) {
    std::size_t parsed = 0;
    const double result = std::stod(value, &parsed);
    if (parsed != value.size() || !std::isfinite(result))
        throw std::invalid_argument("expected a finite floating-point value");
    return result;
}
unsigned long long parseUnsigned(const std::string &value) {
    if (!value.empty() && value.front() == '-')
        throw std::invalid_argument("expected a non-negative integer");
    std::size_t parsed = 0;
    const auto result = std::stoull(value, &parsed);
    if (parsed != value.size())
        throw std::invalid_argument("expected an integer");
    return result;
}
int parseInt(const std::string &value) {
    std::size_t parsed = 0;
    const int result = std::stoi(value, &parsed);
    if (parsed != value.size())
        throw std::invalid_argument("expected an integer");
    return result;
}
} // namespace

void SimParams::validate() const {
    if (vortexCount == 0 && initialConditionFile.empty() && restartFile.empty())
        throw std::invalid_argument("N must be positive");
    if (!(timeStep > 0.0))
        throw std::invalid_argument("timeStep must be positive");
    if (!(endTime >= 0.0))
        throw std::invalid_argument("endTime must be non-negative");
    if (!(outputTime > 0.0))
        throw std::invalid_argument("outputTime must be positive");
    if (!(coreRadius >= 0.0))
        throw std::invalid_argument("coreRadius must be non-negative");
    if (!(absoluteTolerance > 0.0) || !(relativeTolerance >= 0.0))
        throw std::invalid_argument("invalid integration tolerances");
    if (!(minimumTimeStep > 0.0) || !(maximumTimeStep >= minimumTimeStep))
        throw std::invalid_argument("invalid timestep bounds");
    if (numThreads < 0)
        throw std::invalid_argument("numThreads must be non-negative");
    if (boundaryCondition != "infinite" && boundaryCondition != "periodic" &&
        boundaryCondition != "disk")
        throw std::invalid_argument("invalid boundaryCondition");
    if (!(boxLengthX > 0.0) || !(boxLengthY > 0.0) || !(diskRadius > 0.0) ||
        periodicImageLayers < 0 || periodicImageLayers > 64)
        throw std::invalid_argument("invalid geometry parameters or periodicImageLayers > 64");
    if (boundaryCondition != "infinite" && coreRadius != 0.0)
        throw std::invalid_argument("periodic and disk boundaries currently require coreRadius 0");
    if (!(dipoleRemovalDistance > 0.0))
        throw std::invalid_argument("dipoleRemovalDistance must be positive");
    if (!dipoleRemoval && dipoleReinjection != ReinjectionMode::none)
        throw std::invalid_argument("dipoleReinjection requires dipoleRemoval true");
    if (boundaryCondition == "infinite" && dipoleReinjection != ReinjectionMode::none)
        throw std::invalid_argument(
            "dipole reinjection is available only for periodic and disk geometries");
    if (boundaryCondition == "periodic" &&
        std::abs(boxLengthX - boxLengthY) > 1e-13 * std::max(boxLengthX, boxLengthY))
        throw std::invalid_argument("Weiss-McWilliams periodic geometry requires a square box");
    if (outputFile.empty() || diagnosticsFile.empty() || checkpointDirectory.empty())
        throw std::invalid_argument("output and checkpoint paths must not be empty");
    if (outputFile == diagnosticsFile)
        throw std::invalid_argument("outputFile and diagnosticsFile must be different");
}
SimParams loadParams(const std::string &filename) {
    std::ifstream input(filename);
    if (!input)
        throw std::runtime_error("cannot open parameter file: " + filename);
    SimParams p;
    std::optional<std::size_t> legacyNumSteps;
    bool explicitEndTime = false;
    std::string line;
    std::size_t lineNumber = 0;
    // Strict parsing prevents a misspelled option from silently using a default.
    while (std::getline(input, line)) {
        ++lineNumber;
        const auto comment = line.find('#');
        if (comment != std::string::npos)
            line.erase(comment);
        std::istringstream fields(line);
        std::string key, value;
        if (!(fields >> key))
            continue;
        if (!(fields >> value))
            throw std::runtime_error("missing value on parameter line " +
                                     std::to_string(lineNumber));
        try {
            if (key == "N") {
                const auto count = parseUnsigned(value);
                if (count > std::numeric_limits<std::size_t>::max())
                    throw std::out_of_range("N is too large");
                p.vortexCount = static_cast<std::size_t>(count);
            } else if (key == "timeStep")
                p.timeStep = parseDouble(value);
            else if (key == "endTime") {
                p.endTime = parseDouble(value);
                explicitEndTime = true;
            } else if (key == "numSteps")
                legacyNumSteps = parseUnsigned(value);
            else if (key == "outputTime" || key == "OutputTime")
                p.outputTime = parseDouble(value);
            else if (key == "coreRadius")
                p.coreRadius = parseDouble(value);
            else if (key == "coreSize")
                p.coreRadius = std::sqrt(parseDouble(value));
            else if (key == "absoluteTolerance")
                p.absoluteTolerance = parseDouble(value);
            else if (key == "relativeTolerance")
                p.relativeTolerance = parseDouble(value);
            else if (key == "minimumTimeStep")
                p.minimumTimeStep = parseDouble(value);
            else if (key == "maximumTimeStep")
                p.maximumTimeStep = parseDouble(value);
            else if (key == "numThreads")
                p.numThreads = parseInt(value);
            else if (key == "boxLengthX")
                p.boxLengthX = parseDouble(value);
            else if (key == "boxLengthY")
                p.boxLengthY = parseDouble(value);
            else if (key == "diskRadius")
                p.diskRadius = parseDouble(value);
            else if (key == "periodicImageLayers")
                p.periodicImageLayers = parseInt(value);
            else if (key == "randomSeed")
                p.randomSeed = parseUnsigned(value);
            else if (key == "dipoleRemovalDistance")
                p.dipoleRemovalDistance = parseDouble(value);
            else if (key == "initialConditionFile")
                p.initialConditionFile = value;
            else if (key == "restartFile")
                p.restartFile = value;
            else if (key == "checkpointDirectory")
                p.checkpointDirectory = value;
            else if (key == "outputFile")
                p.outputFile = value;
            else if (key == "diagnosticsFile")
                p.diagnosticsFile = value;
            else if (key == "dipoleRemoval" || key == "overwriteOutput" ||
                     key == "overwriteCheckpoints") {
                bool enabled;
                if (value == "true" || value == "1")
                    enabled = true;
                else if (value == "false" || value == "0")
                    enabled = false;
                else
                    throw std::invalid_argument(key + " must be true or false");
                if (key == "dipoleRemoval")
                    p.dipoleRemoval = enabled;
                else if (key == "overwriteOutput")
                    p.overwriteOutput = enabled;
                else
                    p.overwriteCheckpoints = enabled;
            } else if (key == "dipoleReinjection") {
                if (value == "none")
                    p.dipoleReinjection = ReinjectionMode::none;
                else if (value == "independent")
                    p.dipoleReinjection = ReinjectionMode::independent;
                else if (value == "paired")
                    p.dipoleReinjection = ReinjectionMode::paired;
                else
                    throw std::invalid_argument(
                        "dipoleReinjection must be none, independent, or paired");
            } else if (key == "boundaryCondition")
                p.boundaryCondition = value;
            else if (key == "integrator") {
                if (value == "rk4")
                    p.integrator = IntegratorKind::rk4;
                else if (value == "dopri5")
                    p.integrator = IntegratorKind::dopri5;
                else
                    throw std::invalid_argument("integrator must be rk4 or dopri5");
            } else
                throw std::invalid_argument("unknown parameter: " + key);
            std::string trailing;
            if (fields >> trailing)
                throw std::invalid_argument("unexpected extra value: " + trailing);
        } catch (const std::exception &e) {
            throw std::runtime_error("parameter line " + std::to_string(lineNumber) + ": " +
                                     e.what());
        }
    }
    if (legacyNumSteps && !explicitEndTime)
        p.endTime = p.timeStep * static_cast<double>(*legacyNumSteps);
    p.validate();
    return p;
}
VortexSystem loadVortices(const std::string &filename) {
    std::ifstream input(filename);
    if (!input)
        throw std::runtime_error("cannot open initial-condition file: " + filename);
    VortexSystem vortices;
    std::string line;
    std::size_t lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        const auto comment = line.find('#');
        if (comment != std::string::npos)
            line.erase(comment);
        for (char &character : line)
            if (character == ',')
                character = ' ';
        std::istringstream fields(line);
        double x, y, circulation;
        if (!(fields >> x))
            continue;
        if (!(fields >> y >> circulation))
            throw std::runtime_error("invalid initial condition on line " +
                                     std::to_string(lineNumber));
        std::string trailing;
        if (fields >> trailing || !std::isfinite(x) || !std::isfinite(y) ||
            !std::isfinite(circulation))
            throw std::runtime_error("invalid initial condition on line " +
                                     std::to_string(lineNumber));
        vortices.x.push_back(x);
        vortices.y.push_back(y);
        vortices.circulation.push_back(circulation);
    }
    if (vortices.size() == 0)
        throw std::runtime_error("initial-condition file is empty");
    return vortices;
}
void initializeVortices(VortexSystem &vortices, double radius) {
    for (std::size_t i = 0; i < vortices.size(); ++i) {
        const double angle =
            2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(vortices.size());
        vortices.x[i] = radius * std::cos(angle);
        vortices.y[i] = radius * std::sin(angle);
        vortices.circulation[i] = (i % 2 == 0) ? 1.0 : -1.0;
    }
}

void initializePeriodicVortices(VortexSystem &vortices, double lengthX, double lengthY,
                                std::uint64_t seed) {
    if (vortices.size() % 2 != 0)
        throw std::invalid_argument("the built-in periodic initial condition requires an even N");

    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<double> xPosition(-0.5 * lengthX, 0.5 * lengthX);
    std::uniform_real_distribution<double> yPosition(-0.5 * lengthY, 0.5 * lengthY);
    for (std::size_t i = 0; i < vortices.size(); ++i) {
        vortices.x[i] = xPosition(generator);
        vortices.y[i] = yPosition(generator);
        // Equal numbers of positive and negative vortices satisfy periodic neutrality.
        vortices.circulation[i] = (i % 2 == 0) ? 1.0 : -1.0;
    }
}
