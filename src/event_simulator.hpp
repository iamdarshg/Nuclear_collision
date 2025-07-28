#ifndef EVENT_SIMULATOR_HPP
#define EVENT_SIMULATOR_HPP

#include <iostream>
#include <random>
#include <cmath>
#include "matrix_elements.hpp"
#include "potential.hpp"
#include "kinematics.hpp"
#include "glauber.hpp"
#include "dglap.hpp"

namespace Event {

inline void run_simulation(int num_events, double E_per_nucleon_GeV) {
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::uniform_real_distribution<double> impact_b(0.0, 15.0); // fm
    std::uniform_real_distribution<double> Q2_dist(1.0, 100.0); // GeV²
    std::uniform_real_distribution<double> x_dist(0.01, 0.99);

    int collisions = 0;
    int broken_flux_tubes = 0;
    double total_weight = 0.0;

    for (int i = 0; i < num_events; ++i) {
        // --- Sample Q² and Bjorken x
        double Q2 = Q2_dist(rng);
        double x = x_dist(rng);

        // Evolve gluon PDF via DGLAP
        double f_g = DGLAP::gluon_pdf(x, Q2);

        // Compute sqrt(s) from beam energy
        double s = 2 * E_per_nucleon_GeV * 0.938; // p-n CM approx
        double x1 = x, x2 = x;  // parton fractions

        // Kinematics
        double rapidity = Kinematics::rapidity(E_per_nucleon_GeV, x * E_per_nucleon_GeV);
        double x_bj = Kinematics::bjorken_x(Q2, s);

        // Glauber factor from impact parameter
        double b = impact_b(rng);
        double coherence = Glauber::coherence_factor(b);

        // Check QCD scattering amplitude
        double amp = Matrix::compute_amplitude(x1, x2, Q2, s);
        double prob = f_g * f_g * amp * coherence;

        // Include saturation effect (suppress low x)
        if (x_bj < 0.01)
            prob *= std::pow(x_bj / 0.01, 0.3); // toy CGC-like model

        // Simulate string tension between nucleons (e.g. p–p)
        double r_fm = 0.5 + 2.0 * uniform(rng); // randomly 0.5–2.5 fm
        auto string_res = StrongForce::compute_string_potential(r_fm);

        if (uniform(rng) < prob) {
            collisions++;
            total_weight += prob;
        }

        if (string_res.flux_tube_breaks) {
            broken_flux_tubes++;
        }
    }

    std::cout << "Total events simulated: " << num_events << "\n";
    std::cout << "Estimated effective collisions: " << collisions << "\n";
    std::cout << "Avg. event weight: " << total_weight / num_events << "\n";
    std::cout << "Flux tubes broken: " << broken_flux_tubes << "\n";
}

} // namespace Event

#endif
