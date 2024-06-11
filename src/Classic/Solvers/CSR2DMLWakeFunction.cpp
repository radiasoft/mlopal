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

CSR2DMLWakeFunction::CSR2DMLWakeFunction(const std::string& name, std::filesystem::path pyFilepath):
    // TODO(e-carlin): 0 is hard-coded because we don't use nBins_m
    WakeFunction(name, 0),
    pyFilepath_m(pyFilepath)
{}

void CSR2DMLWakeFunction::apply(PartBunchBase<double, 3>* bunch) {
    Inform msg("MLWake ");
    py::scoped_interpreter guard{};
    py::module sys = py::module::import("sys");
    sys.attr("path").attr("append")(pyFilepath_m.parent_path().string().c_str());
    py::module py_module = py::module::import(pyFilepath_m.stem().string().c_str());
    py::function get_wake = py_module.attr("get_wake");
    py::object result = get_wake(bunch); // call the python function
    std::cout << "xxxxxxxxxxxxxxxxxx The result is: " << result.cast<int>() << std::endl; // print the result
}

void CSR2DMLWakeFunction::initialize(const ElementBase* ref) {
    // TODO(e-carlin): del?
}

WakeType CSR2DMLWakeFunction::getType() const {
    return WakeType::CSR2DMLWakeFunction;
}
