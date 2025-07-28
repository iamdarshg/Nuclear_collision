#include "pdf_evolution.hpp"
#include <apfel/apfelxx.h>
#include <memory>

static std::shared_ptr<apfel::GridPDF> pdf_grid;

void PDF::initialize() {
    apfel::InitializeQCD(5, {0.001, 1.3, 4.5, 100.0, 173.1});
    pdf_grid = apfel::BuildPDFs();
}

std::vector<double> PDF::get_flavor_xf(double x, double Q2) {
    auto xf = pdf_grid->EvaluatexQ2(x, Q2);
    return {
        xf[apfel::parton_id::up],
        xf[apfel::parton_id::down],
        xf[apfel::parton_id::strange],
        xf[apfel::parton_id::charm],
        xf[apfel::parton_id::bottom],
        xf[apfel::parton_id::gluon]
    };
}
