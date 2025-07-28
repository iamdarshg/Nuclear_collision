#ifndef GLAUBER_HPP
#define GLAUBER_HPP

namespace Glauber {

// Simple Glauber attenuation factor (impact parameter b ~ 0–15 fm)
inline double coherence_factor(double b_fm) {
    constexpr double R_nucleus = 7.0; // Approx Hg nucleus radius in fm
    constexpr double alpha = 0.2;     // Coherence scaling

    if (b_fm > 2.0 * R_nucleus) return 0.0;

    double x = b_fm / R_nucleus;
    return std::exp(-alpha * x * x);
}

} // namespace Glauber

#endif
