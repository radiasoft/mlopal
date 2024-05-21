//
// Class StatWriter
//   This class writes bunch statistics (*.stat).
//
// Copyright (c) 2019, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
//                     Christof Metzger-Kraus, Open Sourcerer
// All rights reserved
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
#include "StatWriter.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "Utilities/Timer.h"
#include "Physics/Units.h"

#include <sstream>

StatWriter::StatWriter(const std::string& fname, bool restart)
    : StatBaseWriter(fname, restart)
{ }


void StatWriter::fillHeader(const losses_t &losses) {

    if (this->hasColumns()) {
        return;
    }

    columns_m.addColumn("t", "double", "ns", "Time");
    columns_m.addColumn("s", "double", "m", "Path length");
    columns_m.addColumn("numParticles", "long", "1", "Number of Macro Particles");
    columns_m.addColumn("charge", "double", "1", "Bunch Charge");
    columns_m.addColumn("energy", "double", "MeV", "Mean Bunch Energy");

    columns_m.addColumn("rms_x", "double", "m", "RMS Beamsize in x");
    columns_m.addColumn("rms_y", "double", "m", "RMS Beamsize in y");
    columns_m.addColumn("rms_s", "double", "m", "RMS Beamsize in s");

    columns_m.addColumn("rms_px", "double", "1", "RMS Normalized Momenta in x");
    columns_m.addColumn("rms_py", "double", "1", "RMS Normalized Momenta in y");
    columns_m.addColumn("rms_ps", "double", "1", "RMS Normalized Momenta in s");

    columns_m.addColumn("emit_x", "double", "m", "Normalized Emittance x");
    columns_m.addColumn("emit_y", "double", "m", "Normalized Emittance y");
    columns_m.addColumn("emit_s", "double", "m", "Normalized Emittance s");

    columns_m.addColumn("mean_x", "double", "m", "Mean Beam Position in x");
    columns_m.addColumn("mean_y", "double", "m", "Mean Beam Position in y");
    columns_m.addColumn("mean_s", "double", "m", "Mean Beam Position in s");

    columns_m.addColumn("ref_x", "double", "m", "x coordinate of reference particle in lab cs");
    columns_m.addColumn("ref_y", "double", "m", "y coordinate of reference particle in lab cs");
    columns_m.addColumn("ref_z", "double", "m", "z coordinate of reference particle in lab cs");

    columns_m.addColumn("ref_px", "double", "1", "x momentum of reference particle in lab cs");
    columns_m.addColumn("ref_py", "double", "1", "y momentum of reference particle in lab cs");
    columns_m.addColumn("ref_pz", "double", "1", "z momentum of reference particle in lab cs");

    columns_m.addColumn("max_x", "double", "m", "Max Beamsize in x");
    columns_m.addColumn("max_y", "double", "m", "Max Beamsize in y");
    columns_m.addColumn("max_s", "double", "m", "Max Beamsize in s");

    columns_m.addColumn("xpx", "double", "1", "Correlation xpx");
    columns_m.addColumn("ypy", "double", "1", "Correlation ypy");
    columns_m.addColumn("zpz", "double", "1", "Correlation zpz");

    columns_m.addColumn("Dx", "double", "m", "Dispersion in x");
    columns_m.addColumn("DDx", "double", "1", "Derivative of dispersion in x");
    columns_m.addColumn("Dy", "double", "m", "Dispersion in y");
    columns_m.addColumn("DDy", "double", "1", "Derivative of dispersion in y");

    columns_m.addColumn("Bx_ref", "double", "T", "Bx-Field component of ref particle");
    columns_m.addColumn("By_ref", "double", "T", "By-Field component of ref particle");
    columns_m.addColumn("Bz_ref", "double", "T", "Bz-Field component of ref particle");

    columns_m.addColumn("Ex_ref", "double", "MV/m", "Ex-Field component of ref particle");
    columns_m.addColumn("Ey_ref", "double", "MV/m", "Ey-Field component of ref particle");
    columns_m.addColumn("Ez_ref", "double", "MV/m", "Ez-Field component of ref particle");

    columns_m.addColumn("dE", "double", "MeV", "energy spread of the beam");
    columns_m.addColumn("dt", "double", "ns", "time step size");
    columns_m.addColumn("partsOutside", "double",  "1", "outside n*sigma of the beam");
    columns_m.addColumn("DebyeLength", "double",  "m", "Debye length in the boosted frame");
    columns_m.addColumn("plasmaParameter", "double",  "1", "Plasma parameter that gives no. of particles in a Debye sphere");
    columns_m.addColumn("temperature", "double",  "K", "Temperature of the beam");
    columns_m.addColumn("rmsDensity", "double",  "1", "RMS number density of the beam");

    if (Options::computePercentiles) {
        columns_m.addColumn("68_Percentile_x", "double", "m",
                            "68.27 percentile (1 sigma of normal distribution) of x-component of position");
        columns_m.addColumn("68_Percentile_y", "double", "m",
                            "68.27 percentile (1 sigma of normal distribution) of y-component of position");
        columns_m.addColumn("68_Percentile_z", "double", "m",
                            "68.27 percentile (1 sigma of normal distribution) of z-component of position");

        columns_m.addColumn("95_Percentile_x", "double", "m",
                            "95.45 percentile (2 sigma of normal distribution) of x-component of position");
        columns_m.addColumn("95_Percentile_y", "double", "m",
                            "95.45 percentile (2 sigma of normal distribution) of y-component of position");
        columns_m.addColumn("95_Percentile_z", "double", "m",
                            "95.45 percentile (2 sigma of normal distribution) of z-component of position");

        columns_m.addColumn("99_Percentile_x", "double", "m",
                            "99.73 percentile (3 sigma of normal distribution) of x-component of position");
        columns_m.addColumn("99_Percentile_y", "double", "m",
                            "99.73 percentile (3 sigma of normal distribution) of y-component of position");
        columns_m.addColumn("99_Percentile_z", "double", "m",
                            "99.73 percentile (3 sigma of normal distribution) of z-component of position");

        columns_m.addColumn("99_99_Percentile_x", "double", "m",
                            "99.994 percentile (4 sigma of normal distribution) of x-component of position");
        columns_m.addColumn("99_99_Percentile_y", "double", "m",
                            "99.994 percentile (4 sigma of normal distribution) of y-component of position");
        columns_m.addColumn("99_99_Percentile_z", "double", "m",
                            "99.994 percentile (4 sigma of normal distribution) of z-component of position");

        columns_m.addColumn("normalizedEps68Percentile_x", "double", "m",
                            "x-component of normalized emittance at 68 percentile (1 sigma of normal distribution)");
        columns_m.addColumn("normalizedEps68Percentile_y", "double", "m",
                            "y-component of normalized emittance at 68 percentile (1 sigma of normal distribution)");
        columns_m.addColumn("normalizedEps68Percentile_z", "double", "m",
                            "z-component of normalized emittance at 68 percentile (1 sigma of normal distribution)");

        columns_m.addColumn("normalizedEps95Percentile_x", "double", "m",
                            "x-component of normalized emittance at 95 percentile (2 sigma of normal distribution)");
        columns_m.addColumn("normalizedEps95Percentile_y", "double", "m",
                            "y-component of normalized emittance at 95 percentile (2 sigma of normal distribution)");
        columns_m.addColumn("normalizedEps95Percentile_z", "double", "m",
                            "z-component of normalized emittance at 95 percentile (2 sigma of normal distribution)");

        columns_m.addColumn("normalizedEps99Percentile_x", "double", "m",
                            "x-component of normalized emittance at 99 percentile (3 sigma of normal distribution)");
        columns_m.addColumn("normalizedEps99Percentile_y", "double", "m",
                            "y-component of normalized emittance at 99 percentile (3 sigma of normal distribution)");
        columns_m.addColumn("normalizedEps99Percentile_z", "double", "m",
                            "z-component of normalized emittance at 99 percentile (3 sigma of normal distribution)");

        columns_m.addColumn("normalizedEps99_99Percentile_x", "double", "m",
                            "x-component of normalized emittance at 99.99 percentile (4 sigma of normal distribution)");
        columns_m.addColumn("normalizedEps99_99Percentile_y", "double", "m",
                            "y-component of normalized emittance at 99.99 percentile (4 sigma of normal distribution)");
        columns_m.addColumn("normalizedEps99_99Percentile_z", "double", "m",
                            "z-component of normalized emittance at 99.99 percentile (4 sigma of normal distribution)");
    }
    if (Options::dumpBeamMatrix) {
        columns_m.addColumn("S11","double","m**2","element 1,1 of 6D beam matrix");
        columns_m.addColumn("S12","double","m","element 1,2 of 6D beam matrix");
        columns_m.addColumn("S13","double","m**2","element 1,3 of 6D beam matrix");
        columns_m.addColumn("S14","double","m","element 1,4 of 6D beam matrix");
        columns_m.addColumn("S15","double","m**2","element 1,5 of 6D beam matrix");
        columns_m.addColumn("S16","double","m","element 1,6 of 6D beam matrix");
        columns_m.addColumn("S22","double","1","element 2,2 of 6D beam matrix");
        columns_m.addColumn("S23","double","m","element 2,3 of 6D beam matrix");
        columns_m.addColumn("S24","double","1","element 2,4 of 6D beam matrix");
        columns_m.addColumn("S25","double","m","element 2,5 of 6D beam matrix");
        columns_m.addColumn("S26","double","1","element 2,6 of 6D beam matrix");
        columns_m.addColumn("S33","double","m**2","element 3,3 of 6D beam matrix");
        columns_m.addColumn("S34","double","m","element 3,4 of 6D beam matrix");
        columns_m.addColumn("S35","double","m**2","element 3,5 of 6D beam matrix");
        columns_m.addColumn("S36","double","m","element 3,6 of 6D beam matrix");
        columns_m.addColumn("S44","double","1","element 4,4 of 6D beam matrix");
        columns_m.addColumn("S45","double","m","element 4,5 of 6D beam matrix");
        columns_m.addColumn("S46","double","1","element 4,6 of 6D beam matrix");
        columns_m.addColumn("S55","double","m**2","element 5,5 of 6D beam matrix");
        columns_m.addColumn("S56","double","m","element 5,6 of 6D beam matrix");
        columns_m.addColumn("S66","double","1","element 6,6 of 6D beam matrix");
    }
    if (OpalData::getInstance()->isInOPALCyclMode() &&
        Ippl::getNodes() == 1) {
        columns_m.addColumn("R0_x", "double", "m", "R0 Particle position in x");
        columns_m.addColumn("R0_y", "double", "m", "R0 Particle position in y");
        columns_m.addColumn("R0_s", "double", "m", "R0 Particle position in z");

        columns_m.addColumn("P0_x", "double", "1", "R0 Particle momentum in x");
        columns_m.addColumn("P0_y", "double", "1", "R0 Particle momentum in y");
        columns_m.addColumn("P0_s", "double", "1", "R0 Particle momentum in z");
    }

    if (OpalData::getInstance()->isInOPALCyclMode()) {
        columns_m.addColumn("halo_x", "double", "1", "Halo in x");
        columns_m.addColumn("halo_y", "double", "1", "Halo in y");
        columns_m.addColumn("halo_z", "double", "1", "Halo in z");

        columns_m.addColumn("azimuth", "double", "deg",
                            "Azimuth in global coordinates");
    }

    for (size_t i = 0; i < losses.size(); ++ i) {
        columns_m.addColumn(losses[i].first, "long", "1",
                            "Number of lost particles in element");
    }

    if ( mode_m == std::ios::app )
        return;

    OPALTimer::Timer simtimer;
    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    std::stringstream ss;
    ss << "Statistics data '"
       << OpalData::getInstance()->getInputFn()
       << "' " << dateStr << " " << timeStr;

    this->addDescription(ss.str(), "stat parameters");

    this->addDefaultParameters();


    this->addInfo("ascii", 1);
}


void StatWriter::write(const PartBunchBase<double, 3> *beam, Vector_t FDext[],
                       const losses_t &losses, const double& azimuth,
                       const size_t npOutside)
{
    double Ekin = beam->get_meanKineticEnergy();

    double  pathLength = beam->get_sPos();

    /// Write data to files. If this is the first write to the beam statistics file, write SDDS
    /// header information.

    double Q = beam->getCharge();

    if (Ippl::myNode() != 0) {
        return;
    }

    fillHeader(losses);

    this->open();

    this->writeHeader();

    columns_m.addColumnValue("t", beam->getT() * Units::s2ns);   // 1
    columns_m.addColumnValue("s", pathLength);                     // 2
    columns_m.addColumnValue("numParticles", beam->getTotalNum()); // 3
    columns_m.addColumnValue("charge", Q);                         // 4
    columns_m.addColumnValue("energy", Ekin);                      // 5

    columns_m.addColumnValue("rms_x", beam->get_rrms()(0));        // 6
    columns_m.addColumnValue("rms_y", beam->get_rrms()(1));        // 7
    columns_m.addColumnValue("rms_s", beam->get_rrms()(2));        // 8

    columns_m.addColumnValue("rms_px", beam->get_prms()(0));       // 9
    columns_m.addColumnValue("rms_py", beam->get_prms()(1));       // 10
    columns_m.addColumnValue("rms_ps", beam->get_prms()(2));       // 11

    columns_m.addColumnValue("emit_x", beam->get_norm_emit()(0));  // 12
    columns_m.addColumnValue("emit_y", beam->get_norm_emit()(1));  // 13
    columns_m.addColumnValue("emit_s", beam->get_norm_emit()(2));  // 14

    columns_m.addColumnValue("mean_x", beam->get_rmean()(0) );     // 15
    columns_m.addColumnValue("mean_y", beam->get_rmean()(1) );     // 16
    columns_m.addColumnValue("mean_s", beam->get_rmean()(2) );     // 17

    columns_m.addColumnValue("ref_x", beam->RefPartR_m(0));        // 18
    columns_m.addColumnValue("ref_y", beam->RefPartR_m(1));        // 19
    columns_m.addColumnValue("ref_z", beam->RefPartR_m(2));        // 20

    columns_m.addColumnValue("ref_px", beam->RefPartP_m(0));       // 21
    columns_m.addColumnValue("ref_py", beam->RefPartP_m(1));       // 22
    columns_m.addColumnValue("ref_pz", beam->RefPartP_m(2));       // 23

    columns_m.addColumnValue("max_x", beam->get_maxExtent()(0));   // 24
    columns_m.addColumnValue("max_y", beam->get_maxExtent()(1));   // 25
    columns_m.addColumnValue("max_s", beam->get_maxExtent()(2));   // 26

    // Write out Courant Snyder parameters.
    columns_m.addColumnValue("xpx", beam->get_rprms()(0));         // 27
    columns_m.addColumnValue("ypy", beam->get_rprms()(1));         // 28
    columns_m.addColumnValue("zpz", beam->get_rprms()(2));         // 29

    // Write out dispersion.
    columns_m.addColumnValue("Dx",  beam->get_Dx());               // 30
    columns_m.addColumnValue("DDx", beam->get_DDx());              // 31
    columns_m.addColumnValue("Dy",  beam->get_Dy());               // 32
    columns_m.addColumnValue("DDy", beam->get_DDy());              // 33

    // Write head/reference particle/tail field information.
    columns_m.addColumnValue("Bx_ref", FDext[0](0));               // 34 B-ref x
    columns_m.addColumnValue("By_ref", FDext[0](1));               // 35 B-ref y
    columns_m.addColumnValue("Bz_ref", FDext[0](2));               // 36 B-ref z

    columns_m.addColumnValue("Ex_ref", FDext[1](0));               // 37 E-ref x
    columns_m.addColumnValue("Ey_ref", FDext[1](1));               // 38 E-ref y
    columns_m.addColumnValue("Ez_ref", FDext[1](2));               // 39 E-ref z

    columns_m.addColumnValue("dE", beam->getdE());                 // 40 dE energy spread
    columns_m.addColumnValue("dt", beam->getdT() * Units::s2ns); // 41 dt time step size
    columns_m.addColumnValue("partsOutside", npOutside);           // 42 number of particles outside n*sigma

    columns_m.addColumnValue("DebyeLength", beam->get_debyeLength()); // 43 Debye length in the boosted frame
    columns_m.addColumnValue("plasmaParameter", beam->get_plasmaParameter()); // 43 plasma parameter
    columns_m.addColumnValue("temperature", beam->get_temperature()); // 44 Temperature 
    columns_m.addColumnValue("rmsDensity", beam->get_rmsDensity()); // 45 RMS number density

    if (Options::computePercentiles) {
        columns_m.addColumnValue("68_Percentile_x", beam->get_68Percentile()[0]);
        columns_m.addColumnValue("68_Percentile_y", beam->get_68Percentile()[1]);
        columns_m.addColumnValue("68_Percentile_z", beam->get_68Percentile()[2]);

        columns_m.addColumnValue("95_Percentile_x", beam->get_95Percentile()[0]);
        columns_m.addColumnValue("95_Percentile_y", beam->get_95Percentile()[1]);
        columns_m.addColumnValue("95_Percentile_z", beam->get_95Percentile()[2]);

        columns_m.addColumnValue("99_Percentile_x", beam->get_99Percentile()[0]);
        columns_m.addColumnValue("99_Percentile_y", beam->get_99Percentile()[1]);
        columns_m.addColumnValue("99_Percentile_z", beam->get_99Percentile()[2]);

        columns_m.addColumnValue("99_99_Percentile_x", beam->get_99_99Percentile()[0]);
        columns_m.addColumnValue("99_99_Percentile_y", beam->get_99_99Percentile()[1]);
        columns_m.addColumnValue("99_99_Percentile_z", beam->get_99_99Percentile()[2]);

        columns_m.addColumnValue("normalizedEps68Percentile_x", beam->get_normalizedEps_68Percentile()[0]);
        columns_m.addColumnValue("normalizedEps68Percentile_y", beam->get_normalizedEps_68Percentile()[1]);
        columns_m.addColumnValue("normalizedEps68Percentile_z", beam->get_normalizedEps_68Percentile()[2]);

        columns_m.addColumnValue("normalizedEps95Percentile_x", beam->get_normalizedEps_95Percentile()[0]);
        columns_m.addColumnValue("normalizedEps95Percentile_y", beam->get_normalizedEps_95Percentile()[1]);
        columns_m.addColumnValue("normalizedEps95Percentile_z", beam->get_normalizedEps_95Percentile()[2]);

        columns_m.addColumnValue("normalizedEps99Percentile_x", beam->get_normalizedEps_99Percentile()[0]);
        columns_m.addColumnValue("normalizedEps99Percentile_y", beam->get_normalizedEps_99Percentile()[1]);
        columns_m.addColumnValue("normalizedEps99Percentile_z", beam->get_normalizedEps_99Percentile()[2]);

        columns_m.addColumnValue("normalizedEps99_99Percentile_x", beam->get_normalizedEps_99_99Percentile()[0]);
        columns_m.addColumnValue("normalizedEps99_99Percentile_y", beam->get_normalizedEps_99_99Percentile()[1]);
        columns_m.addColumnValue("normalizedEps99_99Percentile_z", beam->get_normalizedEps_99_99Percentile()[2]);
    }
    if (Options::dumpBeamMatrix) {
        columns_m.addColumnValue("S11", beam->getSigmaMatrix()[0][0]);
        columns_m.addColumnValue("S12", beam->getSigmaMatrix()[0][1]);
        columns_m.addColumnValue("S13", beam->getSigmaMatrix()[0][2]);
        columns_m.addColumnValue("S14", beam->getSigmaMatrix()[0][3]);
        columns_m.addColumnValue("S15", beam->getSigmaMatrix()[0][4]);
        columns_m.addColumnValue("S16", beam->getSigmaMatrix()[0][5]);
        columns_m.addColumnValue("S22", beam->getSigmaMatrix()[1][1]);
        columns_m.addColumnValue("S23", beam->getSigmaMatrix()[1][2]);
        columns_m.addColumnValue("S24", beam->getSigmaMatrix()[1][3]);
        columns_m.addColumnValue("S25", beam->getSigmaMatrix()[1][4]);
        columns_m.addColumnValue("S26", beam->getSigmaMatrix()[1][5]);
        columns_m.addColumnValue("S33", beam->getSigmaMatrix()[2][2]);
        columns_m.addColumnValue("S34", beam->getSigmaMatrix()[2][3]);
        columns_m.addColumnValue("S35", beam->getSigmaMatrix()[2][4]);
        columns_m.addColumnValue("S36", beam->getSigmaMatrix()[2][5]);
        columns_m.addColumnValue("S44", beam->getSigmaMatrix()[3][3]);
        columns_m.addColumnValue("S45", beam->getSigmaMatrix()[3][4]);
        columns_m.addColumnValue("S46", beam->getSigmaMatrix()[3][5]);
        columns_m.addColumnValue("S55", beam->getSigmaMatrix()[4][4]);
        columns_m.addColumnValue("S56", beam->getSigmaMatrix()[4][5]);
        columns_m.addColumnValue("S66", beam->getSigmaMatrix()[5][5]);
    }
    if (OpalData::getInstance()->isInOPALCyclMode()) {
        if (Ippl::getNodes() == 1) {
            if (beam->getLocalNum() > 0) {
                columns_m.addColumnValue("R0_x", beam->R[0](0));
                columns_m.addColumnValue("R0_y", beam->R[0](1));
                columns_m.addColumnValue("R0_s", beam->R[0](2));
                columns_m.addColumnValue("P0_x", beam->P[0](0));
                columns_m.addColumnValue("P0_y", beam->P[0](1));
                columns_m.addColumnValue("P0_s", beam->P[0](2));
            } else {
                columns_m.addColumnValue("R0_x", 0.0);
                columns_m.addColumnValue("R0_y", 0.0);
                columns_m.addColumnValue("R0_s", 0.0);
                columns_m.addColumnValue("P0_x", 0.0);
                columns_m.addColumnValue("P0_y", 0.0);
                columns_m.addColumnValue("P0_s", 0.0);
            }
        }
        Vector_t halo = beam->get_halo();
        columns_m.addColumnValue("halo_x", halo(0));
        columns_m.addColumnValue("halo_y", halo(1));
        columns_m.addColumnValue("halo_z", halo(2));

        columns_m.addColumnValue("azimuth", azimuth);
    }

    for(size_t i = 0; i < losses.size(); ++ i) {
        long unsigned int loss = losses[i].second;
        columns_m.addColumnValue(losses[i].first, loss);
    }

    this->writeRow();

    this->close();
}
