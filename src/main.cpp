#include <iostream>
#include "event_simulator.hpp"

int main() {
    constexpr int Z = 80;              // Atomic number of Hg
    constexpr int A = 199;             // Mass number of Hg
    constexpr int total_nucleons = A;  // per nucleus

    constexpr int num_collisions = 100000000;  // Total parton-parton samples
    constexpr double energy_per_nucleon = 0.005; // GeV/u

    std::cout << "=== Hg-199 + Hg-199 Collision Simulation ===\n";
    std::cout << "Using " << total_nucleons << " nucleons per nucleus\n";
    std::cout << "Energy per nucleon: " << energy_per_nucleon << " GeV\n";
    std::cout << "Simulating " << num_collisions << " sampled parton-parton collisions...\n\n";

    Event::run_simulation(num_collisions, energy_per_nucleon);

    return 0;
}
