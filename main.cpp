#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include <chrono> // timer

using namespace std;
using namespace Eigen;

// Function to apply periodic boundary conditions
double periodic(double dl, double bound) {
    return dl - bound * round(dl / bound);
}

// Phi1
//double fuel_f(double n, double k) {
//    return -k * n;
//}

// Fuel function Phi2
double fuel_f(double n, double k) {
    return -k * tanh(n);
}

// Phi3
//double fuel_f(double n, double k) {
//    // steady log(cosh(n))
//    double logcosh;
//    if (std::abs(n) < 20.0) {
//        logcosh = std::log(std::cosh(n));
//    } else {
//        // avoid explosion of cosh(n)
//        logcosh = std::abs(n) + std::log1p(std::exp(-2.0 * std::abs(n))) - std::log(2.0);
//    }
//
//    return -k * n + k * logcosh;
//}

// Phi4
// double fuel_f(double n, double k) {
//    return -k * n / (1 + n);
// }

// Gamma^(1/2) matrix computation
MatrixXd gamma_half_matrix(double theta, double v0, double zeta, double zeta_theta, double zeta_n) {
    double D1 = 1 / zeta_theta;
    double c = 1 / zeta_n;

    double sqrt_zeta = sqrt(zeta);
    double sqrt_D1 = sqrt(D1);
    double term1 = (1 + c * zeta - sqrt(1 - 2 * c * zeta + (c * c + v0 * v0) * zeta * zeta)) / (2 * zeta);
    double term2 = (1 + c * zeta + sqrt(1 - 2 * c * zeta + (c * c + v0 * v0) * zeta * zeta)) / (2 * zeta);

    MatrixXd D_half(4, 4);
    D_half << 1 / sqrt_zeta, 0, 0, 0,
            0, sqrt(term1), 0, 0,
            0, 0, sqrt(term2), 0,
            0, 0, 0, sqrt_D1;

    double cos_theta = cos(theta);
    double sin_theta = sin(theta);

    MatrixXd P(4, 4);
    P << -sin_theta, v0 * zeta * cos_theta * cos_theta / (-1 + c * zeta - sqrt(1 - 2 * c * zeta + (c * c + v0 * v0) * zeta * zeta)),
            v0 * zeta * cos_theta * cos_theta / (-1 + c * zeta + sqrt(1 - 2 * c * zeta + (c * c + v0 * v0) * zeta * zeta)), 0,
            cos_theta, v0 * zeta * cos_theta * sin_theta / (-1 + c * zeta - sqrt(1 - 2 * c * zeta + (c * c + v0 * v0) * zeta * zeta)),
            v0 * zeta * cos_theta * sin_theta / (-1 + c * zeta + sqrt(1 - 2 * c * zeta + (c * c + v0 * v0) * zeta * zeta)), 0,
            0, 0, 0, cos_theta,
            0, cos_theta, cos_theta, 0;

    MatrixXd P_inv(4, 4);
    P_inv << -sin_theta, cos_theta, 0, 0,
            -v0 * zeta / (2 * sqrt(c * c * zeta * zeta - 2 * c * zeta + v0 * v0 * zeta * zeta + 1)),
            -v0 * zeta * tan(theta) / (2 * sqrt(c * c * zeta * zeta - 2 * c * zeta + v0 * v0 * zeta * zeta + 1)),
            0,
            (-c * zeta + sqrt(c * c * zeta * zeta - 2 * c * zeta + v0 * v0 * zeta * zeta + 1) + 1) / (2 * sqrt(c * c * zeta * zeta - 2 * c * zeta + v0 * v0 * zeta * zeta + 1) * cos_theta),
            v0 * zeta / (2 * sqrt(c * c * zeta * zeta - 2 * c * zeta + v0 * v0 * zeta * zeta + 1)),
            v0 * zeta * tan(theta) / (2 * sqrt(c * c * zeta * zeta - 2 * c * zeta + v0 * v0 * zeta * zeta + 1)),
            0,
            (c * zeta + sqrt(c * c * zeta * zeta - 2 * c * zeta + v0 * v0 * zeta * zeta + 1) - 1) / (2 * sqrt(c * c * zeta * zeta - 2 * c * zeta + v0 * v0 * zeta * zeta + 1) * cos_theta),
            0, 0, 1 / cos_theta, 0;

    return P * D_half * P_inv;
}

// Main SDE integrator function
void integrator_SDE_fuel(int T, int N, double dt, double n0, double v0, double q_r, double q_n, double q_t, double zeta_theta, double J0, double R, double L, double zeta, double zeta_n, double k, int seed, int dump) {
    // Initialize random number generator
    random_device rd;
    //mt19937 gen(rd());
    std::mt19937 gen(seed); // fixed random seed

    normal_distribution<> noise_dist(0.0, 1.0);

    // Particle properties
    vector<double> x(N), y(N), x_unwrap(N), y_unwrap(N), phi(N), n(N, n0);
    vector<vector<double>> x_all, y_all, x_unwrap_all, y_unwrap_all, phi_all, n_all, pot_all;

    for (int i = 0; i < N; ++i) {
        x[i] = uniform_real_distribution<>(-L, L)(gen);
        y[i] = uniform_real_distribution<>(-L, L)(gen);
        phi[i] = uniform_real_distribution<>(0, 2 * M_PI)(gen);

        x_unwrap[i] = x[i];
        y_unwrap[i] = y[i];
    }

    auto start = std::chrono::high_resolution_clock::now();

    for (int t = 0; t < T; ++t) {
        vector<double> x_new(N), y_new(N), x_unwrap_new(N), y_unwrap_new(N), phi_new(N), n_new(N), pot_new(N);

        if (t % dump == 0) {
            auto now = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
            std::cout << "Step " << t << ", Elapsed time: " << elapsed << " seconds\n";
        }

        for (int i = 0; i < N; ++i) {
            double force_sum = 0.0;
            double Fi_x = 0.0; // repulsive band around R
            double Fi_y = 0.0;

            double pot_i = 0; // potential

            for (int j = 0; j < N; ++j) {
                if (i != j) {
                    double dx = periodic(x[i] - x[j], 2 * L);
                    double dy = periodic(y[i] - y[j], 2 * L);
                    double dist = sqrt(dx * dx + dy * dy);

                    double lambda_0 = 1e-20; // small lambda_0 indicates a narrowband
                    double dJ_dr = -J0 / (2 * lambda_0) * (1 - tanh((dist - 1) / lambda_0) * tanh((dist - 1) / lambda_0));
                    Fi_x += dJ_dr * (dx / dist) * cos(phi[i] - phi[j]);
                    Fi_y += dJ_dr * (dy / dist) * cos(phi[i] - phi[j]);

                    if (dist < R) {
                        force_sum += J0 * sin(phi[i] - phi[j]);
                        pot_i += -J0 * cos(phi[i] - phi[j]);
                    }
                }
            }

            MatrixXd gamma_half = gamma_half_matrix(phi[i], v0, zeta, zeta_theta, zeta_n);
            VectorXd random_noises(4);
            //for (int k0 = 0; k0 < 4; ++k0) random_noises[k0] = noise_dist(gen);
	    	random_noises[0] = sqrt(2.0 * q_t) * noise_dist(gen);
            random_noises[1] = sqrt(2.0 * q_t) * noise_dist(gen);
            random_noises[2] = sqrt(2.0 * q_r) * noise_dist(gen);
            random_noises[3] = sqrt(2.0 * q_n) * noise_dist(gen);

            VectorXd correlated_noises = gamma_half * random_noises;

            double noise_theta = correlated_noises[2];
            double Delta_phi = (-1.0 / zeta_theta) * force_sum * dt + sqrt(dt) * noise_theta;
            phi_new[i] = phi[i] + Delta_phi;

            double noise_n = correlated_noises[3];
            double Delta_n = fuel_f(n[i], k) * dt + sqrt(dt) * noise_n;
            n_new[i] = n[i] + Delta_n;

            double noise_r_x = correlated_noises[0];
            double noise_r_y = correlated_noises[1];

            double Delta_x = - fuel_f(n[i], k) * v0 * cos(phi[i]) * dt - v0 * cos(phi[i]) * sqrt(dt) * noise_n + sqrt(dt) * noise_r_x + 1.0 / zeta * Fi_x * dt;
            double Delta_y = - fuel_f(n[i], k) * v0 * sin(phi[i]) * dt - v0 * sin(phi[i]) * sqrt(dt) * noise_n + sqrt(dt) * noise_r_y + 1.0 / zeta * Fi_y * dt;
            x_new[i] = x[i] + Delta_x;
            y_new[i] = y[i] + Delta_y;

            x_new[i] = periodic(x_new[i], 2 * L);
            y_new[i] = periodic(y_new[i], 2 * L);

            x_unwrap_new[i] = x_unwrap[i] + Delta_x;
            y_unwrap_new[i] = y_unwrap[i] + Delta_y;

            pot_new[i] = pot_i;
        }

        if (t % dump == 0) {
            x_all.push_back(x);
            y_all.push_back(y);
            phi_all.push_back(phi);
            n_all.push_back(n);

            x_unwrap_all.push_back(x_unwrap);
            y_unwrap_all.push_back(y_unwrap);

            pot_all.push_back(pot_new);
        }

        x = x_new;
        y = y_new;
        phi = phi_new;
        n = n_new;

        x_unwrap = x_unwrap_new;
        y_unwrap = y_unwrap_new;
    }

    // The instantaneous angle
    int time_gap = 20;
    int T_dump = T / dump - time_gap;

    vector<vector<double>> theta_inst(N); // instantaneous angle: N by T_dump matrix

    for (int p = 0 ; p < N; ++p) {
        for (int t = 0; t < T_dump; ++t) {
            double inst_x = x_unwrap_all[t + time_gap][p] - x_unwrap_all[t][p];
            double inst_y = y_unwrap_all[t + time_gap][p] - y_unwrap_all[t][p];
            double theta_pt = atan2(inst_y, inst_x);
            theta_inst[p].push_back(theta_pt);
        }
    }

    // Generate filename based on parameters
    stringstream filename;
    filename << "./simulations/v0=" << v0 << "_qr=" << q_r << "_qn=" << q_n << "_qt=" << q_t
             << "_N=" << N << "_R=" << R << "_L=" << 2 * L << "_zeta=" << zeta << "_zeta_theta=" << zeta_theta
             << "_zeta_n=" << zeta_n << "_J0=" << J0 << "_n=" << n0 << "_k=" << k << "_seed=" << seed << ".dump";

    ofstream dump_file(filename.str());

    for (int t = 0; t < T_dump; ++t) {
    // for (int t = T_dump - 3000; t < T_dump; ++t) {
    	dump_file << "ITEM: TIMESTEP\n";
        dump_file << t << "\nITEM: NUMBER OF ATOMS\n" << N << "\nITEM: BOX BOUNDS pp pp pp\n";
        dump_file << -L << " " << L << "\n";
        dump_file << -L << " " << L << "\n";
        dump_file << -L << " " << L << "\n";
        dump_file << "ITEM: ATOMS id x y x_unwrap y_unwrap n phi cos_phi sin_phi inst_angle mux muy pot\n";

        for (int i = 0; i < N; ++i) {
            dump_file << i + 1 << " " << x_all[t][i] << " " << y_all[t][i] << " " << x_unwrap_all[t][i] << " " << y_unwrap_all[t][i]
                      << " " << n_all[t][i] << " " << phi_all[t][i] << " " << cos(phi_all[t][i]) << " " << sin(phi_all[t][i])
                      << " " << theta_inst[i][t] <<  " " << cos(theta_inst[i][t]) << " " << sin(theta_inst[i][t])
                      << " " << pot_all[t][i] << "\n";
        }

    }
    dump_file.close();
}

int main() {
    // parameters
    double dt = 0.01;
    int T = 302000;
    int N = 1000;
    double n0 = 4000.0;
    double v0 = 0.24;
    double q_r = 1.0 / 80.0;
    double q_n = 1.0 / 80.0;
    double q_t = 1.0 / 80.0;
    double zeta_theta = 125.0 / 144.0;
    double zeta = 8.0;
    double J0 = 175.0 / 4608.0;
    double R = 1.0;
    double L = 15.81 / 2.0;
    double zeta_n = 8.0;
    double k = 2.0;
    int seed = 10000;

    int dump = static_cast<int>(round(1.0 / dt));

    integrator_SDE_fuel(T, N, dt, n0, v0, q_r, q_n, q_t, zeta_theta, J0, R, L, zeta, zeta_n, k, seed, dump);

    return 0;
}
