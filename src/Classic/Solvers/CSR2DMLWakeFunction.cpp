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
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/SBend.h"
#include "Algorithms/PartBunchBase.h"
#include "Solvers/CSR2DMLWakeFunction.h"
#include "AbstractObjects/OpalData.h"

#include <filesystem>
#include <pybind11/embed.h>
#include <pybind11/stl.h> // type conversion
namespace py = pybind11;

// TODO(e-carlin): needed?
#include "Utilities/Util.h"

CSR2DMLWakeFunction::CSR2DMLWakeFunction(
    const std::string& name,
    std::filesystem::path pyFilepath,
    const unsigned int& N
    ):
    WakeFunction(name, N),
    bendRadius_m(0.0),
    planeDensity_m(),
    pyFilepath_m(pyFilepath),
    totalBendAngle_m(0.0)
{}

void CSR2DMLWakeFunction::apply(PartBunchBase<double, 3>* bunch) {
    Inform msg("MLWake ");
    std::cout << "1111111111111111" << std::endl;
    std::pair<double, double> meshInfoX;
    std::pair<double, double> meshInfoY;
    // TODO(e-carlin): one nBins_m
    bunch->calcPlaneDensity(
        nBins_m,
        nBins_m,
        planeDensity_m,
        meshInfoX,
        meshInfoY
    );
    Vector_t smin, smax;
    bunch->get_bounds(smin, smax);
    std::cout << "22222222222222222" << std::endl;
    py::scoped_interpreter guard{};
    py::object result = getWakeFn()(
        planeDensity_m,
        bendRadius_m,
        totalBendAngle_m,
        smax(0) - smin(0),
        smax(2) - smin(2)
    );
    std::cout << "xxxxxxxxxxxxxxxxxx The result is: " << result.cast<int>() << std::endl;
}

py::function CSR2DMLWakeFunction::getWakeFn() {
    py::module::import("sys").attr("path").attr("append")(pyFilepath_m.parent_path().string().c_str());
    return py::module::import(pyFilepath_m.stem().string().c_str()).attr("get_wake");
}

void CSR2DMLWakeFunction::initialize(const ElementBase* ref) {
    // TODO(e-carlin): copied from CSRWakeFunction.cpp
    // TODO(e-carlin): do we need all of these values?
    if (ref->getType() == ElementType::RBEND ||
        ref->getType() == ElementType::SBEND) {

        const Bend2D *bend = static_cast<const Bend2D *>(ref);
        // TODO(e-carlin): used below. Needed?
        // double End;

        bendRadius_m = bend->getBendRadius();
        totalBendAngle_m = std::abs(bend->getBendAngle());
        // TODO(e-carlin): needed?
        // bend->getDimensions(Begin_m, End);
        // Length_m = bend->getEffectiveLength();
        // FieldBegin_m = bend->getEffectiveCenter() - Length_m / 2.0;
        // bendName_m = bend->getName();
    }
}

WakeType CSR2DMLWakeFunction::getType() const {
    return WakeType::CSR2DMLWakeFunction;
}
