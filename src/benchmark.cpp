#include "compute.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
int main(int argc, char **argv) {
    try {
        if (argc > 3)
            throw std::invalid_argument("expected vortex count and repetition count");
        const auto positiveInteger = [](const char *text) {
            std::size_t parsed = 0;
            const std::string token(text);
            if (token.empty() || token.front() == '-')
                throw std::invalid_argument("benchmark arguments must be positive integers");
            const auto value = std::stoull(token, &parsed);
            if (parsed != token.size() || value == 0 ||
                value > static_cast<unsigned long long>(std::numeric_limits<int>::max()))
                throw std::invalid_argument(
                    "benchmark arguments must be positive integers no greater than INT_MAX");
            return static_cast<int>(value);
        };
        const std::size_t count = argc > 1 ? positiveInteger(argv[1]) : 2000;
        const int repeats = argc > 2 ? positiveInteger(argv[2]) : 5;
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
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
