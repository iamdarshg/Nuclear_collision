#ifndef MATRIX_ELEMENTS_HPP
#define MATRIX_ELEMENTS_HPP

#include <cmath>
#include <array>
#include <random>

// Define type alias for 4-momentum (E, px, py, pz)
using FourVec = std::array<double, 4>;

namespace Physics {

constexpr double alpha_em_0 = 1.0 / 137.0;
constexpr double Lambda_QED2 = 0.04; // ~0.2 GeV^2
constexpr double alpha_s = 0.118; // QCD coupling at ~Z pole
constexpr double pi = 3.141592653589793;

// u from s, t
inline double u_from_st(double s, double t) {
    return -s - t;
}

inline double running_alpha_em(double Q2) {
    // One-loop vacuum polarization correction
    return alpha_em_0 / (1.0 - (alpha_em_0 / (3.0 * pi)) * std::log(Q2 / Lambda_QED2));
}

// Mandelstam variables
inline double s_hat(const FourVec& p1, const FourVec& p2) {
    return std::pow(p1[0] + p2[0], 2) - std::pow(p1[1] + p2[1], 2)
           - std::pow(p1[2] + p2[2], 2) - std::pow(p1[3] + p2[3], 2);
}

// Tree-level qq → qq via gluon exchange
inline double matrix_element_qq_to_qq(double s, double t) {
    return (4.0 / 9.0) * 4.0 * pi * alpha_s * alpha_s * (s * s + u_from_st(s, t) * u_from_st(s, t)) / (t * t);
}

// Tree-level gg → q q̄
inline double matrix_element_gg_to_qqbar(double s) {
    return (pi * alpha_s * alpha_s / (3.0 * s)) * std::log(s / 0.1);
}

// Tree-level γγ → e⁺e⁻
inline double matrix_element_photon_photon_to_fermions(double s) {
    double alpha = running_alpha_em(s);
    return 4.0 * pi * alpha * alpha / s;
}



}



#include <cmath>
#include "matrix_elements.hpp"
#include "pdf_evolution.hpp"

#include <cmath>
#include <iostream>
#include <array>

const double alpha_em = 1.0 / 137.0;
const double alpha_s = 0.118; // Running should be added later
const double hbar_c = 0.1973269804; // GeV·fm
const double pi = 3.14159265358979323846;

struct FourVector {
    double E, px, py, pz;
};

namespace Matrix {

    // Mandelstam variables
    double compute_s(const FourVector& p1, const FourVector& p2) {
        return std::pow(p1.E + p2.E, 2) - std::pow(p1.px + p2.px, 2)
             - std::pow(p1.py + p2.py, 2) - std::pow(p1.pz + p2.pz, 2);
    }

    double compute_t(const FourVector& p1, const FourVector& p3) {
        return std::pow(p1.E - p3.E, 2) - std::pow(p1.px - p3.px, 2)
             - std::pow(p1.py - p3.py, 2) - std::pow(p1.pz - p3.pz, 2);
    }

    double compute_u(const FourVector& p1, const FourVector& p4) {
        return std::pow(p1.E - p4.E, 2) - std::pow(p1.px - p4.px, 2)
             - std::pow(p1.py - p4.py, 2) - std::pow(p1.pz - p4.pz, 2);
    }

    // Placeholder for full multi-parton matrix element calculation
    double compute_qcd_amplitude(double s, double t, double u) {
        // Basic gluon exchange approximation (extend with full Feynman diagrams)
        double prefactor = (4.0 * pi * alpha_s);
        double denom = (t * t + u * u);
        if (denom == 0.0) denom = 1e-8; // Avoid div by zero
        return prefactor * (s * s + u * u) / denom;
    }

    double compute_qed_amplitude(double s, double t, double u) {
        // Electron-level scattering amplitude as a template
        double prefactor = (4.0 * pi * alpha_em);
        double denom = (t * t + u * u);
        if (denom == 0.0) denom = 1e-8;
        return prefactor * (s * s + u * u) / denom;
    }

    // Full corrected scattering amplitude using QCD + QED + PDF
    double compute_amplitude(double x1, double x2, double Q2, double s) {
        auto xf1 = PDF::get_flavor_xf(x1, Q2);
        auto xf2 = PDF::get_flavor_xf(x2, Q2);

        double amp = 0.0;

        // Loop over all flavor combinations (u,d,s,c,b,g)
        for (size_t i = 0; i < xf1.size(); ++i) {
            for (size_t j = 0; j < xf2.size(); ++j) {
                double partonic_s = x1 * x2 * s;
                double t = -Q2;
                double u = -partonic_s - t;

                double amp_qcd = compute_qcd_amplitude(partonic_s, t, u);
                double amp_qed = compute_qed_amplitude(partonic_s, t, u);

                // Include flavor weights from PDFs
                double weight = xf1[i] * xf2[j];
                amp += weight * (amp_qcd + amp_qed);
            }
        }

        return amp;
    }

    double compute_cross_section(double amplitude, double s) {
        return (1.0 / (16.0 * pi * s * s)) * std::abs(amplitude * amplitude) * hbar_c * hbar_c;
    }

}
// namespace Matrix

#endif
