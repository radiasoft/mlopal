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

// TODO(e-carlin): needed?
#ifndef MLWAKEFUNCTION_HH
#define MLWAKEFUNCTION_HH

#include "Solvers/WakeFunction.h"

#include <vector>

class Filter;
class ElementBase;

class MLWakeFunction: public WakeFunction {
public:
    MLWakeFunction(const std::string& name, const unsigned int& N);

    void apply(PartBunchBase<double, 3>* bunch) override;

    void initialize(const ElementBase* ref) override;

    virtual WakeType getType() const override;
};

#endif //MLWAKEFUNCTION_HH
