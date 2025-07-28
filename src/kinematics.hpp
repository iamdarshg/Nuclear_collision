#ifndef KINEMATICS_HPP
#define KINEMATICS_HPP

#include <cmath>

namespace Kinematics {

// Rapidity from pz and E
inline double rapidity(double E, double pz) {
    return 0.5 * std::log((E + pz) / (E - pz + 1e-12));
}

// Bjorken-x approximation (for 2-nucleon frame)
inline double bjorken_x(double Q2, double s) {
    return Q2 / (s + 1e-12); // +eps for numerical safety
}

} // namespace Kinematics

#endif
