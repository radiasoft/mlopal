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
#include <filesystem>
#include <pybind11/embed.h>
namespace py = pybind11;

// TODO(e-carlin): needed?
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "Utilities/Util.h"

CSR2DMLWakeFunction::CSR2DMLWakeFunction(
    const std::string& name,
    std::filesystem::path pyFilepath,
    const unsigned int& N
    ):
    WakeFunction(name, N),
    planeDensity_m(),
    pyFilepath_m(pyFilepath)
{}

void CSR2DMLWakeFunction::apply(PartBunchBase<double, 3>* bunch) {
    Inform msg("MLWake ");
    // bunch->calcLineDensity(nBins_m, lineDensity_m, meshInfo);
    std::cout << "1111111111111111" << std::endl;
    std::pair<double, double> meshInfoX;
    std::pair<double, double> meshInfoY;
    // TODO(e-carlin): one nBins_m
    bunch->calcPlaneDensity(
        nBins_m,
        nBins_m,
        planeDensity_m,
        // bunch,
        meshInfoX,
        meshInfoY
    );
    // bunch->calcLineDensity(nBins_m, lineDensity_m, meshInfo);
    std::cout << "22222222222222222" << std::endl;
    py::scoped_interpreter guard{};
    py::object result = getWakeFn()("Im the bunch object");
    std::cout << "xxxxxxxxxxxxxxxxxx The result is: " << result.cast<int>() << std::endl;
}

py::function CSR2DMLWakeFunction::getWakeFn() {
    py::module::import("sys").attr("path").attr("append")(pyFilepath_m.parent_path().string().c_str());
    return py::module::import(pyFilepath_m.stem().string().c_str()).attr("get_wake");
}

void CSR2DMLWakeFunction::initialize(const ElementBase* ref) {
    // TODO(e-carlin): del?
}

WakeType CSR2DMLWakeFunction::getType() const {
    return WakeType::CSR2DMLWakeFunction;
}
