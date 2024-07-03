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

#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <vector>


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
    std::pair<double, double> meshInfoZ;
    // TODO(e-carlin): nBins for X and for Z
    bunch->calcPlaneDensity(
        nBins_m,
        nBins_m,
        planeDensity_m,
        meshInfoX,
        meshInfoZ
    );
    std::string planeDensityFile = writePlaneDensity(bunch, nBins_m, nBins_m);
    std::string command = "python " + pyFilepath_m.string() + ' ' + planeDensityFile;
    int ret = std::system(command.c_str());
    if (WEXITSTATUS(ret) != 0) {
        throw std::runtime_error("command=" + command + " exited with error=" + std::to_string(ret));
    }
    std::vector<std::vector<double>> wakeX;
    std::vector<std::vector<double>> wakeZ;
    std::tie(wakeX, wakeZ) = readWakes(planeDensityFile, nBins_m, nBins_m);

    // Apply wake to particles in bunch
    const double &meshSpacingX = meshInfoX.second;
    const double &meshSpacingZ = meshInfoZ.second;
    const double &meshOriginX = meshInfoX.first + 0.5 * meshSpacingX;
    const double &meshOriginZ = meshInfoZ.first + 0.5 * meshSpacingZ;
    for (unsigned int i = 0; i < bunch->getLocalNum(); ++i) {
        const Vector_t &R = bunch->R[i];
        double distanceToOriginX = (R(0) - meshOriginX) / meshSpacingX;
        double distanceToOriginZ = (R(2) - meshOriginZ) / meshSpacingZ;

        unsigned int indexX = (unsigned int)floor(distanceToOriginX);
        unsigned int indexZ = (unsigned int)floor(distanceToOriginZ);
        double leverX = distanceToOriginX - indexX;
        double leverZ = distanceToOriginZ - indexZ;

        PAssert_LT(indexX, planeDensity_m.size() - 1);
        PAssert_LT(indexZ, planeDensity_m[0].size() - 1);

        bunch->Ef[i](0) += (1. - leverX) * wakeX[indexX][indexZ] + leverX * wakeX[indexX + 1][indexZ];
        bunch->Ef[i](2) += (1. - leverZ) * wakeZ[indexX][indexZ] + leverZ * wakeZ[indexX + 1][indexZ];
    }
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


std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> CSR2DMLWakeFunction::readWakes(std::filesystem::path planeDensityFile, unsigned int nBinsX, unsigned int nBinsZ) {
    // TODO(e-carlin): clean up all repetition in here.
    std::vector<std::vector<double>> wakeX(nBinsX, std::vector<double>(nBinsZ));
    std::vector<std::vector<double>> wakeZ(nBinsX, std::vector<double>(nBinsZ));
    std::ifstream infile(planeDensityFile.parent_path() / "wake.bin", std::ios::in | std::ios::binary);
    int num_elements = nBinsX * nBinsZ;
    std::vector<double> wakeFlatX(num_elements);
    std::vector<double> waktFlatZ(num_elements);

    infile.read((char*)&waktFlatZ[0], num_elements*sizeof(double));
    // Reshape the flat vector into a 2D vector
    for (unsigned int i = 0; i < nBinsX; ++i) {
        for (unsigned int j = 0; j < nBinsZ; ++j) {
            wakeZ[i][j] = waktFlatZ[i * nBinsZ + j];
        }
    }

    infile.read((char*)&wakeFlatX[0], num_elements*sizeof(double));
    // Reshape the flat vector into a 2D vector
    for (unsigned int i = 0; i < nBinsX; ++i) {
        for (unsigned int j = 0; j < nBinsZ; ++j) {
            wakeX[i][j] = wakeFlatX[i * nBinsZ + j];
        }
    }


    infile.close();
    return std::make_tuple(wakeX, wakeZ);
}

std::string CSR2DMLWakeFunction::writePlaneDensity(PartBunchBase<double, 3>* bunch, unsigned int nBinsX, unsigned int nBinsZ) {
    Vector_t smin, smax;
    bunch->get_bounds(smin, smax);

    std::string dataFile = "data.bin";
    std::ofstream outfile(dataFile, std::ios::out | std::ios::binary);

    // Header with dimensions of array
    double rows = static_cast<double>(planeDensity_m.size());
    double cols = static_cast<double>(planeDensity_m[0].size());
    assert(rows == nBinsX);
    assert(cols == nBinsZ);
    outfile.write((char*)&rows, sizeof(double));
    outfile.write((char*)&cols, sizeof(double));

    // Write out array
    for (const auto &row : planeDensity_m) {
        outfile.write((char*)&row[0], row.size()*sizeof(double));
    }

    // Scalars
    outfile.write((char*)&bendRadius_m, sizeof(double));
    outfile.write((char*)&totalBendAngle_m, sizeof(double));
    double xSpan = smax(0) - smin(0);
    outfile.write((char*)&xSpan, sizeof(double));
    double zSpan = smax(2) - smin(2);
    outfile.write((char*)&zSpan, sizeof(double));
    outfile.close();
    return std::filesystem::absolute(dataFile);
}
