#ifndef POTENTIAL_HPP
#define POTENTIAL_HPP

#include <cmath>
#include <stdexcept>

namespace StrongForce {

// Constants (in natural units, GeV and fm)
constexpr double sigma = 0.9;          // String tension (GeV/fm)
constexpr double kappa = 0.25;         // Coulomb-like term coefficient
constexpr double r_break = 1.5;        // Flux string break threshold (fm)
constexpr double pair_energy_threshold = 2.0 * 0.33; // light quark pair mass

struct StringForceResult {
    double potential;
    bool flux_tube_breaks;
};

// Compute confining potential with flux tube breaking
inline StringForceResult compute_string_potential(double r_fm) {
    StringForceResult result;

    if (r_fm <= 0.01) {
        result.potential = 1e6; // Regularization at short distance
        result.flux_tube_breaks = false;
        return result;
    }

    double V = sigma * r_fm + kappa / r_fm;

    result.potential = V;
    result.flux_tube_breaks = (r_fm > r_break && V > pair_energy_threshold);
    return result;
}

} // namespace StrongForce

#endif
