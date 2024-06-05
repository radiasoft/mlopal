//
// Class CSR2DMLWakeFunction
//
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Solvers/CSR2DMLWakeFunction.h"
#include <pybind11/embed.h>
namespace py = pybind11;

// TODO(e-carlin): needed?
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "Utilities/Util.h"
#include <torch/torch.h>

CSR2DMLWakeFunction::CSR2DMLWakeFunction(const std::string& name,
                                 const unsigned int& N):
    WakeFunction(name, N)
{}

void CSR2DMLWakeFunction::apply(PartBunchBase<double, 3>* bunch) {
    Inform msg("MLWake ");
    // TODO(e-carlin): call model
    py::scoped_interpreter guard{}; // start the interpreter

    py::module sys = py::module::import("sys");
    sys.attr("path").attr("append")("."); // add the current directory to Python's sys.path

    py::module py_module = py::module::import("py_module"); // import the python module
    py::function py_function = py_module.attr("py_function"); // get the python function

    int a = 5, b = 10;
    py::object result = py_function(a, b); // call the python function

    std::cout << "xxxxxxxxxxxxxxxxxx The result is: " << result.cast<int>() << std::endl; // print the result
}

void CSR2DMLWakeFunction::initialize(const ElementBase* ref) {
    // TODO(e-carlin): del?
}

WakeType CSR2DMLWakeFunction::getType() const {
    return WakeType::CSR2DMLWakeFunction;
}
