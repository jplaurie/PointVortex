
std::vector<std::pair<double, double>> computeVelocities(
    const std::vector<Vortex>& vortices, double coreSize = 1e-5) {
    
    int N = vortices.size();
    std::vector<std::pair<double, double>> velocity(N, {0.0, 0.0});
    
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        double vx = 0.0, vy = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            double dx = vortices[i].x - vortices[j].x;
            double dy = vortices[i].y - vortices[j].y;
            double r2 = dx * dx + dy * dy + coreSize;

            double coeff = vortices[j].gamma / (2 * M_PI * r2);
            vx += coeff * (-dy);
            vy += coeff * dx;
        }
        velocity[i] = {vx, vy};
    }

    return velocity;
}




bool rkf45Step(std::vector<Vortex>& vortices, double& dt, double tol = 1e-6) {
    const int N = vortices.size();
    const double safety = 0.9, min_dt = 1e-10, max_dt = 1.0;

    std::vector<Vortex> y0 = vortices;

    auto computeStage = [&](const std::vector<Vortex>& state, double h) {
        return computeVelocities(state);
    };

    auto k1 = computeVelocities(y0);

    std::vector<Vortex> yk2 = y0;
    for (int i = 0; i < N; ++i) {
        yk2[i].x += dt * 0.25 * k1[i].first;
        yk2[i].y += dt * 0.25 * k1[i].second;
    }
    auto k2 = computeVelocities(yk2);

    std::vector<Vortex> yk3 = y0;
    for (int i = 0; i < N; ++i) {
        yk3[i].x += dt * (3.0/32.0 * k1[i].first + 9.0/32.0 * k2[i].first);
        yk3[i].y += dt * (3.0/32.0 * k1[i].second + 9.0/32.0 * k2[i].second);
    }
    auto k3 = computeVelocities(yk3);

    std::vector<Vortex> yk4 = y0;
    for (int i = 0; i < N; ++i) {
        yk4[i].x += dt * (1932.0/2197.0 * k1[i].first - 7200.0/2197.0 * k2[i].first + 7296.0/2197.0 * k3[i].first);
        yk4[i].y += dt * (1932.0/2197.0 * k1[i].second - 7200.0/2197.0 * k2[i].second + 7296.0/2197.0 * k3[i].second);
    }
    auto k4 = computeVelocities(yk4);

    std::vector<Vortex> yk5 = y0;
    for (int i = 0; i < N; ++i) {
        yk5[i].x += dt * (439.0/216.0 * k1[i].first - 8.0 * k2[i].first + 3680.0/513.0 * k3[i].first - 845.0/4104.0 * k4[i].first);
        yk5[i].y += dt * (439.0/216.0 * k1[i].second - 8.0 * k2[i].second + 3680.0/513.0 * k3[i].second - 845.0/4104.0 * k4[i].second);
    }
    auto k5 = computeVelocities(yk5);

    std::vector<Vortex> yk6 = y0;
    for (int i = 0; i < N; ++i) {
        yk6[i].x += dt * (-8.0/27.0 * k1[i].first + 2.0 * k2[i].first - 3544.0/2565.0 * k3[i].first
                        + 1859.0/4104.0 * k4[i].first - 11.0/40.0 * k5[i].first);
        yk6[i].y += dt * (-8.0/27.0 * k1[i].second + 2.0 * k2[i].second - 3544.0/2565.0 * k3[i].second
                        + 1859.0/4104.0 * k4[i].second - 11.0/40.0 * k5[i].second);
    }
    auto k6 = computeVelocities(yk6);

    // Estimate next state and error
    double maxError = 0.0;

    for (int i = 0; i < N; ++i) {
        double dx4 = dt * (25.0/216.0 * k1[i].first + 1408.0/2565.0 * k3[i].first + 2197.0/4104.0 * k4[i].first - 0.2 * k5[i].first);
        double dy4 = dt * (25.0/216.0 * k1[i].second + 1408.0/2565.0 * k3[i].second + 2197.0/4104.0 * k4[i].second - 0.2 * k5[i].second);

        double dx5 = dt * (16.0/135.0 * k1[i].first + 6656.0/12825.0 * k3[i].first + 28561.0/56430.0 * k4[i].first
                         - 9.0/50.0 * k5[i].first + 2.0/55.0 * k6[i].first);
        double dy5 = dt * (16.0/135.0 * k1[i].second + 6656.0/12825.0 * k3[i].second + 28561.0/56430.0 * k4[i].second
                         - 9.0/50.0 * k5[i].second + 2.0/55.0 * k6[i].second);

        double err = std::hypot(dx5 - dx4, dy5 - dy4);
        maxError = std::max(maxError, err);
    }

    // Adaptive timestep control
    if (maxError < tol) {
        // Accept step
        for (int i = 0; i < N; ++i) {
            vortices[i].x += dt * (16.0/135.0 * k1[i].first + 6656.0/12825.0 * k3[i].first
                                 + 28561.0/56430.0 * k4[i].first - 9.0/50.0 * k5[i].first + 2.0/55.0 * k6[i].first);
            vortices[i].y += dt * (16.0/135.0 * k1[i].second + 6656.0/12825.0 * k3[i].second
                                 + 28561.0/56430.0 * k4[i].second - 9.0/50.0 * k5[i].second + 2.0/55.0 * k6[i].second);
        }

        // Adjust dt
        dt = std::min(max_dt, safety * dt * std::pow(tol / maxError, 0.25));
        return true;
    } else {
        // Reject step, decrease dt
        dt = std::max(min_dt, safety * dt * std::pow(tol / maxError, 0.25));
        return false;
    }
}