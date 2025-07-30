#include <vector>
#ifndef PDF_EVOLUTION_HPP
#define PDF_EVOLUTION_HPP
#include "pdf_evolution.cpp" 
namespace PDF {
    void initialize();
    std::vector<double> get_flavor_xf(double x, double Q2);
}

#endif // PDF_EVOLUTION_HPP
