#include "checkpoint.h"
#include "compute.h"
#include "print.h"
#include "read.h"
#include "timestep.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct GeometryMetadata {
    double lengthX = 0.0;
    double lengthY = 0.0;
    int imageLayers = 0;
};

GeometryMetadata geometryMetadata(const SimParams &params) {
    if (params.boundaryCondition == "periodic")
        return {params.boxLengthX, params.boxLengthY, params.periodicImageLayers};
    if (params.boundaryCondition == "disk")
        return {params.diskRadius, 0.0, 0};
    return {};
}

std::unique_ptr<VelocityKernel> makeKernel(const SimParams &params) {
    if (params.boundaryCondition == "infinite")
        return std::make_unique<InfinitePlaneKernel>(params.coreRadius);
    if (params.boundaryCondition == "periodic")
        return std::make_unique<PeriodicBoxKernel>(params.boxLengthX, params.boxLengthY,
                                                   params.periodicImageLayers);
    return std::make_unique<DiskKernel>(params.diskRadius);
}

bool checkpointMatches(const Checkpoint &checkpoint, const SimParams &params,
                       const GeometryMetadata &geometry) {
    return checkpoint.coreRadius == params.coreRadius &&
           checkpoint.integrator == params.integrator &&
           checkpoint.boundaryCondition == params.boundaryCondition &&
           checkpoint.geometryLengthX == geometry.lengthX &&
           checkpoint.geometryLengthY == geometry.lengthY &&
           checkpoint.periodicImageLayers == geometry.imageLayers &&
           checkpoint.dipoleRemoval == params.dipoleRemoval &&
           (!params.dipoleRemoval ||
            (checkpoint.dipoleRemovalDistance == params.dipoleRemovalDistance &&
             checkpoint.dipoleReinjection == params.dipoleReinjection));
}

VortexSystem makeInitialState(const SimParams &params) {
    if (!params.initialConditionFile.empty())
        return loadVortices(params.initialConditionFile);

    VortexSystem vortices(params.vortexCount);
    if (params.boundaryCondition == "periodic") {
        initializePeriodicVortices(vortices, params.boxLengthX, params.boxLengthY,
                                   params.randomSeed);
        return vortices;
    }

    // Keep the demonstration ring away from the disk wall.
    double radius = 1.0;
    if (params.boundaryCondition == "disk")
        radius = 0.5 * params.diskRadius;
    initializeVortices(vortices, radius);
    return vortices;
}

} // namespace

int main(int argc, char **argv) {
    try {
        const std::string parameterFile = argc > 1 ? argv[1] : "params.txt";
        const SimParams params = loadParams(parameterFile);

#ifdef _OPENMP
        if (params.numThreads > 0)
            omp_set_num_threads(params.numThreads);
#endif

        const bool restarting = !params.restartFile.empty();
        Checkpoint restart;
        if (restarting)
            restart = loadCheckpoint(params.restartFile);

        const GeometryMetadata geometry = geometryMetadata(params);
        if (restarting && !checkpointMatches(restart, params, geometry))
            throw std::runtime_error(
                "checkpoint geometry or integrator settings do not match parameters");

        VortexSystem vortices = restarting ? std::move(restart.vortices) : makeInitialState(params);
        auto kernel = makeKernel(params);
        const Invariants initial =
            restarting ? restart.initialInvariants : computeInvariants(vortices, *kernel);
        DipoleManager dipoles =
            restarting ? DipoleManager(params, restart.dipoleState) : DipoleManager(params);
        if (!restarting)
            dipoles.process(vortices);
        Invariants segmentReference = restarting && restart.hasSegmentInvariants
                                          ? restart.segmentInvariants
                                          : computeInvariants(vortices, *kernel);
        RungeKuttaIntegrator integrator(vortices.size());
        VelocityField velocity(vortices.size());

        // Detect the common rerun/restart collision before opening and possibly truncating CSVs.
        const std::size_t firstCheckpointIndex = restarting ? restart.outputIndex + 1 : 0;
        if ((!restarting || restart.time < params.endTime) && !params.overwriteCheckpoints) {
            const auto firstCheckpoint =
                checkpointPath(params.checkpointDirectory, firstCheckpointIndex);
            if (std::filesystem::exists(firstCheckpoint))
                throw std::runtime_error("refusing to overwrite checkpoint: " +
                                         firstCheckpoint.string());
        }
        TrajectoryWriter trajectory(params.outputFile, params.overwriteOutput);
        DiagnosticsWriter diagnostics(params.diagnosticsFile, initial, params.overwriteOutput);

        double time = restarting ? restart.time : 0.0;
        double dt = restarting ? restart.suggestedTimeStep
                               : std::clamp(params.timeStep, params.minimumTimeStep,
                                            params.maximumTimeStep);
        double nextOutput = restarting ? restart.nextOutputTime : params.outputTime;
        std::size_t acceptedSteps = restarting ? restart.acceptedSteps : 0;
        std::size_t outputIndex = restarting ? restart.outputIndex : 0;

        if (time > params.endTime)
            throw std::runtime_error("checkpoint time is later than endTime");

        const auto writeFrame = [&] {
            kernel->evaluate(vortices, velocity);
            trajectory.write(time, vortices, velocity);
            const Invariants current = computeInvariants(vortices, *kernel);
            const DipoleEventState eventState = dipoles.state();
            diagnostics.write(time, current, segmentReference, eventState.removedPairs,
                              eventState.reinjectedPairs);
            printDiagnostics(time, acceptedSteps, current, initial, params.boundaryCondition,
                             segmentReference, eventState.removedPairs, eventState.reinjectedPairs);
        };
        const auto writeCurrentCheckpoint = [&] {
            writeCheckpoint(params.checkpointDirectory, vortices, initial, time, dt, nextOutput,
                            acceptedSteps, outputIndex, params.coreRadius, params.integrator,
                            params.boundaryCondition, geometry.lengthX, geometry.lengthY,
                            geometry.imageLayers, params.dipoleRemoval,
                            params.dipoleRemovalDistance, params.dipoleReinjection, dipoles.state(),
                            segmentReference, params.overwriteCheckpoints);
        };

        // A restarted branch records its starting frame but does not duplicate its source
        // checkpoint.
        writeFrame();
        if (!restarting)
            writeCurrentCheckpoint();

        while (time < params.endTime) {
            // Clipping to both boundaries makes CSV and checkpoint times exactly reproducible.
            const double stepSize = std::min({dt, params.endTime - time, nextOutput - time});
            double acceptedStep = stepSize;

            if (params.integrator == IntegratorKind::rk4) {
                integrator.rk4Step(vortices, stepSize, *kernel);
            } else {
                const StepResult result =
                    integrator.dopri5Step(vortices, stepSize, *kernel, params);
                acceptedStep = result.acceptedTimeStep;
                dt = result.suggestedTimeStep;
            }

            time += acceptedStep;
            ++acceptedSteps;
            if (dipoles.process(vortices) != 0)
                segmentReference = computeInvariants(vortices, *kernel);

            const double roundingSlack =
                16.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(time));
            if (time + roundingSlack >= nextOutput || time + roundingSlack >= params.endTime) {
                while (nextOutput <= time + roundingSlack)
                    nextOutput += params.outputTime;
                ++outputIndex;
                writeFrame();
                writeCurrentCheckpoint();
            }
        }

        return 0;
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
