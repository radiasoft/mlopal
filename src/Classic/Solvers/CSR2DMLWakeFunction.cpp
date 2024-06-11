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
#include <filesystem>
namespace py = pybind11;

// TODO(e-carlin): needed?
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "Utilities/Util.h"
#include <torch/torch.h>

CSR2DMLWakeFunction::CSR2DMLWakeFunction(const std::string& name, std::filesystem::path pyFilepath):
    // TODO(e-carlin): 0 is hard-coded because we don't use nBins_m
    WakeFunction(name, 0),
    pyFilepath_m(pyFilepath)
{}

void CSR2DMLWakeFunction::apply(PartBunchBase<double, 3>* bunch) {
    Inform msg("MLWake ");
    std::cout << "111111111111111111111111111111111" << pyFilepath_m << std::endl;
    py::scoped_interpreter guard{}; // start the interpreter
    // // TODO(e-carlin): call model
    std::cout << "22222222222222222222222222222222222222" << pyFilepath_m.parent_path() << std::endl;
    py::module sys = py::module::import("sys");
    sys.attr("path").attr("append")(pyFilepath_m.parent_path().string().c_str());
    std::cout << "33333333333333333333333333333333333" << pyFilepath_m.stem() << std::endl;
    py::module py_module = py::module::import(pyFilepath_m.stem().string().c_str()); // import the python module
    std::cout << "55555555555555555555555555555555555" << std::endl;
    // TODO(e-carlin): change the name of the func
    py::function py_function = py_module.attr("get_wake"); // get the python function
    std::cout << "66666666666666666666666666666666666" << std::endl;
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
