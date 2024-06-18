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
    // py::object result = getWakeFn()(
    //     planeDensity_m,
    //     bendRadius_m,
    //     totalBendAngle_m,
    //     smax(0) - smin(0),
    //     smax(2) - smin(2)
    // );
    // std::cout << "The result is: " << result.cast<int>() << std::endl;
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
