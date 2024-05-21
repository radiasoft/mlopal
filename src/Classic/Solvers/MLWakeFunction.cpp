//
// Class MLWakeFunction
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
#include "Solvers/MLWakeFunction.h"

// TODO(e-carlin): needed?
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "Utilities/Util.h"

MLWakeFunction::MLWakeFunction(const std::string& name,
                                 const unsigned int& N):
    WakeFunction(name, N)
{}

void MLWakeFunction::apply(PartBunchBase<double, 3>* bunch) {
    Inform msg("MLWake ");

    // TODO(e-carlin): call model
}

void MLWakeFunction::initialize(const ElementBase* ref) {
    // TODO(e-carlin): del?
}

WakeType MLWakeFunction::getType() const {
    return WakeType::MLWakeFunction;
}
