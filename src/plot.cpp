#include "Python.h"
#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "event_simulator.hpp"

PYBIND11_MODULE(plot, m) {
    m.doc() = "Event simulator for nuclear collisions";

    m.def("simulate_collision_event", &simulate_collision_event, "Simulate a collision event",
          py::arg("v_initial"), py::arg("impact_parameter_val"));
}
