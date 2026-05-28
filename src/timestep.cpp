
std::vector<std::pair<double, double>> computeVelocities(
    const VortexSystem& vortices, SimParams params) {
    
    
    std::vector<std::pair<double, double>> velocity(params.N, {0.0, 0.0});
    
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        double vx = 0.0, vy = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            double dx = vortices[i].x - vortices[j].x;
            double dy = vortices[i].y - vortices[j].y;
            double r2 = dx * dx + dy * dy + coreSize;

            double coeff = vortices[j].circ / (2.0 * std::numbers::pi * r2);
            vx += coeff * (-dy);
            vy += coeff * dx;
        }
        velocity[i] = {vx, vy};
    }

    return velocity;
}




bool rkf45Step(VortexSystem& vortices, double& dt, double tol = 1e-6) {
 
    const double SAFETY = 0.9;
    const double MIN_SCALE = 0.2;
    const double MAX_SCALE = 5.0;

    // Coefficients for DOPRI5 (Dormand-Prince)
    const double c[] = {0, 1.0/5, 3.0/10, 4.0/5, 8.0/9, 1, 1};
    const double a[][7] = {
        {},
        {1.0/5},
        {3.0/40, 9.0/40},
        {44.0/45, -56.0/15, 32.0/9},
        {19372.0/6561, -25360.0/2187, 64448.0/6561, -212.0/729},
        {9017.0/3168, -355.0/33, 46732.0/5247, 49.0/176, -5103.0/18656},
        {35.0/384, 0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84}
    };

    const double b[] = {35.0/384, 0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84, 0};
    const double bAlt[] = {5179.0/57600, 0, 7571.0/16695, 393.0/640, -92097.0/339200, 187.0/2100, 1.0/40};



    for (int i = 1; i < 7; ++i) {
            State yi = y;
            for (int j = 0; j < i; ++j)
                for (size_t m = 0; m < y.size(); ++m)
                    yi[m] += h * a[i][j] * k[j][m];
            f(t + c[i]*h, yi, k[i]);
        }




   VortexSystem y0 = vortices;

    auto computeStage = [&](const VortexSystem& state, double h) {
        return computeVelocities(state);
    };

    auto k1 = computeVelocities(y0,params);

  VortexSystem yk2 = y0;


    for (int i = 0; i < N; ++i) {
        yk2[i].x += dt * a[1][0] * k1[i].first;
        yk2[i].y += dt * a[1][0] * k1[i].second;
    }
    
    
    auto k2 = computeVelocities(yk2);

    VortexSystem yk3 = y0;
    for (int i = 0; i < N; ++i) {
        yk3[i].x += dt * (a[2][0] * k1[i].first + a[2][1] * k2[i].first);
        yk3[i].y += dt * (a[2][0] * k1[i].second + a[2][1] * k2[i].second);
    }
    auto k3 = computeVelocities(yk3);

    VortexSystem yk4 = y0;
    for (int i = 0; i < N.params; ++i) {
        yk4[i].x += dt * (a[3][0] * k1[i].first + a[3][1] * k2[i].first + a[3][2]* k3[i].first);
        yk4[i].y += dt * (a[3][0] * k1[i].second + a[3][1] * k2[i].second + a[3][2] * k3[i].second);
    }
    auto k4 = computeVelocities(yk4);

    VortexSystem yk5 = y0;
    for (int i = 0; i < N; ++i) {
        yk5[i].x += dt * (a[4][0]* k1[i].first + a[4][1] * k2[i].first + a[4][2] * k3[i].first + a[4][3] * k4[i].first);
        yk5[i].y += dt * (a[4][0]* k1[i].second + a[4][1] * k2[i].second + a[4][2] * k3[i].second + a[4][3] * k4[i].second);
    }
    auto k5 = computeVelocities(yk5);

    VortexSystem yk6 = y0;
    for (int i = 0; i < N; ++i) {
        yk6[i].x += dt * (a[5][0] * k1[i].first + a[5][1] * k2[i].first  + a[5][2] * k3[i].first + a[5][3] * k4[i].first + a[5][4] * k5[i].first);
        yk6[i].y += dt * (a[5][0] * k1[i].second + a[5][1] * k2[i].second  + a[5][2] * k3[i].second + a[5][3] * k4[i].second + a[5][4] * k5[i].second);
    auto k6 = computeVelocities(yk6);


    VortexSystem yk7 = y0;
    for (int i = 0; i < N; ++i) {
        yk7[i].x += dt * (a[6][0] * k1[i].first + a[6][1]* k2[i].first + a[6][2]* k3[i].first + a[6][3] * k4[i].first +a[6][4] * k5[i].first + a[6][5] * K6[i].first);
        yk7[i].y += dt * (a[6][0] * k1[i].second + a[6][1]* k2[i].second + a[6][2]* k3[i].second + a[6][3] * k4[i].second +a[6][4] * k5[i].second + a[6][5] * K6[i].second);
    auto k7 = computeVelocities(yk6);



    // Estimate next state and error
    double maxError = 0.0;

    for (int i = 0; i < N; ++i) {
        double dx4 = dt * (bAlt[0] * k1[i].first +  bAlt[1] * k2[i].first + bAlt[2] * k3[i].first +  bAlt[3] * k4[i].first + bAlt[4] * k5[i].first);
        double dy4 = dt * (bAlt[0] * k1[i].second  +  bAlt[1] * k2[i].second + bAlt[2] * k3[i].second + bAlt[3] * k4[i].second + bAlt[4] * k5[i].second);

        double dx5 = dt * (b[0] * k1[i].first + b[1] * k2[i].first + b[2] * k3[i].first + b[3] * k4[i].first + b[4] * k5[i].first + b[5] * k6[i].first);
        double dy5 = dt * (b[0] * k1[i].second + b[1] * k2[i].second + b[2] * k3[i].second + b[3] * k4[i].second + b[4] * k5[i].second + b[5] * k6[i].second);

        double err = std::hypot(dx5 - dx4, dy5 - dy4);
        maxError = std::max(maxError, err);
    }

    // Adaptive timestep control
    if (maxError < tol) {
        // Accept step
        for (int i = 0; i < N; ++i) {
            vortices[i].u += (b[0] * k1[i].first + b[1] * k2[i].first+ b[2] * k3[i].first + b[3] * k4[i].first + b[4] * k5[i].first + b[5] * k6[i].first);
            vortices[i].v += (b[0] * k1[i].second + b[1] * k2[i].second + b[2] * k3[i].second + b[3] * k4[i].second + b[4] * k5[i].second + b[5] * k6[i].second);
            vortices[i].x += dt * vortices[i].u;
            vortices[i].y += dt * vortices[i].v;
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


void simulate(VortexSystem& vortices, double& dt, double totalTime, double tol = 1e-6) {
    double t = 0.0;
    int step = 0;

    while (t < totalTime) {
        bool success = rkf45Step(vortices, dt, tol);
        if (success) {
            t += dt;
            if (step % 10 == 0) {
                std::cout << "t = " << t << " | dt = " << dt << " | step " << step << "\n";
            }
            ++step;
        } else {
            std::cout << "Retrying with smaller dt = " << dt << "\n";
        }
    }
}





 // Compute velocities using OpenMP
        #pragma omp parallel for
        for (int i = 0; i < params.N; ++i) {
            for (int j = 0; j < params.N; ++j) {
                if (i == j) continue;

                double dx = vortices[i].x - vortices[j].x;
                double dy = vortices[i].y - vortices[j].y;
                double r2 = dx * dx + dy * dy + params.coreSize;

                double coeff = vortices[j].circ / (2 * params.PI * r2);

                vortices[i].u+= coeff * (-dy);
                vortices[i].v += coeff * dx;
            }
            
        }

        // Update positions
        for (int i = 0; i < params.N; ++i) {
            vortices[i].x += dt * vortices[i].u;
            vortices[i].y += dt * vortices[i].v;
        }