#include "print.h"
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
TrajectoryWriter::TrajectoryWriter(const std::string &filename, bool overwrite) {
    if (!overwrite && std::filesystem::exists(filename))
        throw std::runtime_error("refusing to overwrite output file: " + filename);
    output_.open(filename, std::ios::trunc);
    if (!output_)
        throw std::runtime_error("cannot open output file: " + filename);
    output_ << "time,index,x,y,circulation,u,v\n";
    output_ << std::setprecision(17);
    output_.flush();
}
void TrajectoryWriter::write(double time, const VortexSystem &vortices,
                             const VelocityField &velocity) {
    if (velocity.x.size() != vortices.size() || velocity.y.size() != vortices.size())
        throw std::invalid_argument("velocity and vortex arrays have different lengths");
    for (std::size_t i = 0; i < vortices.size(); ++i)
        output_ << time << ',' << i << ',' << vortices.x[i] << ',' << vortices.y[i] << ','
                << vortices.circulation[i] << ',' << velocity.x[i] << ',' << velocity.y[i] << '\n';
    output_.flush();
    if (!output_)
        throw std::runtime_error("failed while writing trajectory output");
}
DiagnosticsWriter::DiagnosticsWriter(const std::string &filename, const Invariants &initial,
                                     bool overwrite)
    : initial_(initial) {
    if (!overwrite && std::filesystem::exists(filename))
        throw std::runtime_error("refusing to overwrite diagnostics file: " + filename);
    output_.open(filename, std::ios::trunc);
    if (!output_)
        throw std::runtime_error("cannot open diagnostics file: " + filename);
    output_ << "time,circulation,linear_impulse_x,linear_impulse_y,angular_impulse,hamiltonian,"
               "delta_circulation,delta_linear_impulse_x,delta_linear_impulse_y,"
               "delta_angular_impulse,delta_hamiltonian,segment_delta_circulation,"
               "segment_delta_linear_impulse_x,segment_delta_linear_impulse_y,"
               "segment_delta_angular_impulse,segment_delta_hamiltonian,removed_pairs,"
               "reinjected_pairs\n";
    output_ << std::setprecision(17);
    output_.flush();
}
void DiagnosticsWriter::write(double time, const Invariants &value,
                              const Invariants &segmentReference, std::size_t removedPairs,
                              std::size_t reinjectedPairs) {
    output_ << time << ',' << value.circulation << ',' << value.linearImpulseX << ','
            << value.linearImpulseY << ',' << value.angularImpulse << ',' << value.hamiltonian
            << ',' << value.circulation - initial_.circulation << ','
            << value.linearImpulseX - initial_.linearImpulseX << ','
            << value.linearImpulseY - initial_.linearImpulseY << ','
            << value.angularImpulse - initial_.angularImpulse << ','
            << value.hamiltonian - initial_.hamiltonian << ','
            << value.circulation - segmentReference.circulation << ','
            << value.linearImpulseX - segmentReference.linearImpulseX << ','
            << value.linearImpulseY - segmentReference.linearImpulseY << ','
            << value.angularImpulse - segmentReference.angularImpulse << ','
            << value.hamiltonian - segmentReference.hamiltonian << ',' << removedPairs << ','
            << reinjectedPairs << '\n';
    output_.flush();
    if (!output_)
        throw std::runtime_error("failed while writing diagnostics output");
}
void printDiagnostics(double time, std::size_t steps, const Invariants &value,
                      const Invariants &initial, const std::string &boundaryCondition,
                      const Invariants &segmentReference, std::size_t removedPairs,
                      std::size_t reinjectedPairs) {
    std::cout << std::setprecision(10) << "time=" << time << " steps=" << steps
              << " circulation=" << value.circulation
              << " dCirculation=" << value.circulation - initial.circulation
              << " segmentDCirculation=" << value.circulation - segmentReference.circulation
              << " H=" << value.hamiltonian << " dH=" << value.hamiltonian - initial.hamiltonian
              << " segmentDH=" << value.hamiltonian - segmentReference.hamiltonian;

    if (boundaryCondition == "infinite" || boundaryCondition == "periodic") {
        std::cout << " Ix=" << value.linearImpulseX
                  << " dIx=" << value.linearImpulseX - initial.linearImpulseX
                  << " segmentDIx=" << value.linearImpulseX - segmentReference.linearImpulseX
                  << " Iy=" << value.linearImpulseY
                  << " dIy=" << value.linearImpulseY - initial.linearImpulseY
                  << " segmentDIy=" << value.linearImpulseY - segmentReference.linearImpulseY;
    }
    if (boundaryCondition == "infinite" || boundaryCondition == "disk") {
        std::cout << " L=" << value.angularImpulse
                  << " dL=" << value.angularImpulse - initial.angularImpulse
                  << " segmentDL=" << value.angularImpulse - segmentReference.angularImpulse;
    }
    std::cout << " removedPairs=" << removedPairs << " reinjectedPairs=" << reinjectedPairs << '\n'
              << std::flush;
}
