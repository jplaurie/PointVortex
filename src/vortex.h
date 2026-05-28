#include <vector>
#include <cstddef>

struct VortexSystem {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> circ;
    std::vector<double> u;
    std::vector<double> v;
    
    explicit VortexSystem(std::size_t N):
        x(N),
        y(N),
        circ(N),
        u(N),
        v(N)
    {}

    std::size_t size() const {
        return x.size();
    }

};