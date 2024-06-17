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

// TODO(e-carlin): needed?
#ifndef CSR2DMLWAKEFUNCTION_HH
#define CSR2DMLWAKEFUNCTION_HH

#include "Solvers/WakeFunction.h"

#include <vector>
#include <filesystem>
#include <pybind11/embed.h>
namespace py = pybind11;

class Filter;
class ElementBase;

class CSR2DMLWakeFunction: public WakeFunction {
public:
    CSR2DMLWakeFunction(const std::string& name, const std::filesystem::path pyFilepath, const unsigned int& N);

    void apply(PartBunchBase<double, 3>* bunch) override;


    void initialize(const ElementBase* ref) override;

    virtual WakeType getType() const override;
private:
    std::vector<std::vector<double>> planeDensity_m;
    /// python file that calls ml model
    std::filesystem::path pyFilepath_m;
    // Radius of curvature of effective dipole.
    double bendRadius_m;
    double totalBendAngle_m;
    py::function getWakeFn();
};

#endif //CSR2DMLWAKEFUNCTION_HH
