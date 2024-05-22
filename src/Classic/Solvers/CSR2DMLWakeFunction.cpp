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

// TODO(e-carlin): needed?
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "Utilities/Util.h"

CSR2DMLWakeFunction::CSR2DMLWakeFunction(const std::string& name,
                                 const unsigned int& N):
    WakeFunction(name, N)
{}

void CSR2DMLWakeFunction::apply(PartBunchBase<double, 3>* bunch) {
    Inform msg("MLWake ");

// std::pair<std::vector<std::vector<double>>, std::map<std::string, std::vector<double>>> convertToModelInput(PartBunchBase& partBunch) {
    // std::vector<std::vector<double>> lambda_distribution = partBunch.Q;  // Assuming Q is a 2D array

    // double x_min = std::numeric_limits<double>::max();
    // double x_max = std::numeric_limits<double>::min();
    // double z_min = std::numeric_limits<double>::max();
    // double z_max = std::numeric_limits<double>::min();

    // // Find spatial extent of the bunch
    // for (const auto& position : partBunch.R) {
    //     x_min = std::min(x_min, position.x);
    //     x_max = std::max(x_max, position.x);
    //     z_min = std::min(z_min, position.z);
    //     z_max = std::max(z_max, position.z);
    // }

    // // TODO(e-carlin): pretty unsure about these.
    // std::map<std::string, std::vector<double>> scalars = {
    //     {"s", {std::accumulate(partBunch.Q.begin(), partBunch.Q.end(), 0.0)}},
    //     {"Sx", {x_max - x_min}},
    //     {"Sz", {z_max - z_min}},
    //     {"rho_max", {*std::max_element(partBunch.Q.begin(), partBunch.Q.end())}}
    // };

    // return {lambda_distribution, scalars};
    // TODO(e-carlin): call model
}

void CSR2DMLWakeFunction::initialize(const ElementBase* ref) {
    // TODO(e-carlin): del?
}

WakeType CSR2DMLWakeFunction::getType() const {
    return WakeType::CSR2DMLWakeFunction;
}
