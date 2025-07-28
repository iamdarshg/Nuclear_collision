#ifndef DGLAP_HPP
#define DGLAP_HPP

#include <cmath>

namespace DGLAP {

// Constants
constexpr double Q0 = 1.0;  // GeV
constexpr double lambda = 0.3;  // Evolution rate (Λ_QCD-like)
constexpr double n_q = 5.0;     // Quark valence power

// Gluon PDF: grows as Q² increases, drops at high-x
inline double gluon_pdf(double x, double Q2) {
    double evolution = std::log(Q2 / Q0 + 1.0);
    return (1.0 + lambda * evolution) * std::pow(1.0 - x, 5) / x;
}

// Quark PDF: valence-like (non-singlet) structure function
inline double quark_pdf(double x, double Q2) {
    double evolution = std::log(Q2 / Q0 + 1.0);
    return (1.0 + 0.5 * lambda * evolution) * std::pow(x, 0.5) * std::pow(1.0 - x, n_q);
}

// Sample PDF based on flavor
inline double get_pdf(const std::string& parton, double x, double Q2) {
    if (parton == "g") return gluon_pdf(x, Q2);
    if (parton == "u" || parton == "d" || parton == "s") return quark_pdf(x, Q2);
    return 0.0;
}

} // namespace DGLAP

#endif
