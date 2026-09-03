#include "compute.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
int main(int argc, char **argv) {
    const std::size_t count = argc > 1 ? std::stoull(argv[1]) : 2000;
    const int repeats = argc > 2 ? std::stoi(argv[2]) : 5;
    VortexSystem vortices(count);
    std::mt19937_64 random(1234567);
    std::uniform_real_distribution<double> position(-1.0, 1.0);
    for (std::size_t i = 0; i < count; ++i) {
        vortices.x[i] = position(random);
        vortices.y[i] = position(random);
        vortices.circulation[i] = (i % 2 == 0) ? 1.0 : -1.0;
    }
    VelocityField velocity(count);
    InfinitePlaneKernel kernel;
    kernel.evaluate(vortices, velocity);
    const auto start = std::chrono::steady_clock::now();
    for (int repetition = 0; repetition < repeats; ++repetition)
        kernel.evaluate(vortices, velocity);
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    const double interactions = static_cast<double>(repeats) * count * (count - 1);
    std::cout << std::setprecision(6) << "N=" << count << " repeats=" << repeats
              << " seconds=" << seconds << " interactions_per_second=" << interactions / seconds
              << '\n';
}
