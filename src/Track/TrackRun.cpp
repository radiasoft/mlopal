//
// Class TrackRun
//   The RUN command.
//
// Copyright (c) 200x - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Track/TrackRun.h"

#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/ObjectFunction.h"
#include "AbstractObjects/OpalData.h"

#include "Algorithms/ParallelCyclotronTracker.h"
#include "Algorithms/ParallelTTracker.h"
#include "Algorithms/ThickTracker.h"
#include "Algorithms/Tracker.h"

#include "Attributes/Attributes.h"

#include "Beamlines/TBeamline.h"

#include "BasicActions/Option.h"

#include "Distribution/Distribution.h"

#include "Physics/Physics.h"
#include "Physics/Units.h"

#include "Track/Track.h"

#include "Utilities/OpalException.h"

#include "Structure/Beam.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/DataSink.h"
#include "Structure/FieldSolver.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/H5PartWrapperForPC.h"
#include "Structure/H5PartWrapperForPT.h"

#include "OPALconfig.h"
#include "changes.h"

#include <boost/assign.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>


extern Inform *gmsg;

std::shared_ptr<Tracker> TrackRun::itsTracker = nullptr;

namespace {
    // The attributes of class TrackRun.
    enum {
        METHOD,           // Tracking method to use.
        TURNS,            // The number of turns to be tracked.
        MBMODE,           // The working way for multi-bunch mode for OPAL-cycl: "FORCE" or "AUTO"
        PARAMB,           // The control parameter for "AUTO" mode of multi-bunch,
        MB_ETA,           // The scale parameter for binning in multi-bunch mode
        MB_BINNING,       // The binning type in multi-bunch mode
        BEAM,             // The beam to track
        FIELDSOLVER,      // The field solver attached
        BOUNDARYGEOMETRY, // The boundary geometry
        DISTRIBUTION,     // The particle distribution
        TRACKBACK,
        SIZE
    };
}

const std::string TrackRun::defaultDistribution("DISTRIBUTION");

const boost::bimap<TrackRun::RunMethod, std::string> TrackRun::stringMethod_s =
    boost::assign::list_of<const boost::bimap<TrackRun::RunMethod, std::string>::relation>
    (RunMethod::PARALLELT,  "PARALLEL-T")
    (RunMethod::CYCLOTRONT, "CYCLOTRON-T")
    (RunMethod::THICK,      "THICK");


TrackRun::TrackRun():
    Action(SIZE, "RUN",
           "The \"RUN\" sub-command tracks the defined particles through "
           "the given lattice."),
    dist(nullptr),
    fs(nullptr),
    ds(nullptr),
    phaseSpaceSink_m(nullptr),
    isFollowupTrack_m(false),
    method_m(RunMethod::NONE),
    macromass_m(0.0),
    macrocharge_m(0.0) {
    itsAttr[METHOD] = Attributes::makePredefinedString
                      ("METHOD", "Name of tracking algorithm to use.",
                       {"THICK", "OPAL-T", "PARALLEL-T", "OPAL-CYCL", "CYCLOTRON-T"});
    
    itsAttr[TURNS] = Attributes::makeReal
        ("TURNS", "Number of turns to be tracked; Number of neighboring bunches to be tracked in cyclotron.", 1.0);

    itsAttr[MBMODE] = Attributes::makePredefinedString
        ("MBMODE", "The working way for multi-bunch mode for OPAL-cycl.",
         {"FORCE", "AUTO"}, "FORCE");

    itsAttr[PARAMB] = Attributes::makeReal
        ("PARAMB", "Control parameter to define when to start multi-bunch mode, only available in \"AUTO\" mode.", 5.0);

    itsAttr[MB_ETA] = Attributes::makeReal
        ("MB_ETA", "The scale parameter for binning in multi-bunch mode.", 0.01);

    itsAttr[MB_BINNING] = Attributes::makePredefinedString
        ("MB_BINNING", "Type of energy binning in multi-bunch mode.",
         {"GAMMA_BINNING", "BUNCH_BINNING"}, "GAMMA_BINNING");

    itsAttr[BEAM] = Attributes::makeString
        ("BEAM", "Name of beam.");

    itsAttr[FIELDSOLVER] = Attributes::makeString
        ("FIELDSOLVER", "Field solver to be used.");

    itsAttr[BOUNDARYGEOMETRY] = Attributes::makeString
        ("BOUNDARYGEOMETRY", "Boundary geometry to be used NONE (default).", "NONE");

    itsAttr[DISTRIBUTION] = Attributes::makeStringArray
        ("DISTRIBUTION", "List of particle distributions to be used.");

    itsAttr[TRACKBACK] = Attributes::makeBool
        ("TRACKBACK", "Track in reverse direction, default: false.", false);

    registerOwnership(AttributeHandler::SUB_COMMAND);
    opal = OpalData::getInstance();
}


TrackRun::TrackRun(const std::string& name, TrackRun* parent):
    Action(name, parent),
    dist(nullptr),
    fs(nullptr),
    ds(nullptr),
    phaseSpaceSink_m(nullptr),
    isFollowupTrack_m(false),
    method_m(RunMethod::NONE),
    macromass_m(0.0),
    macrocharge_m(0.0) {
    opal = OpalData::getInstance();
}


TrackRun::~TrackRun() {
    delete phaseSpaceSink_m;
}


TrackRun* TrackRun::clone(const std::string& name) {
    return new TrackRun(name, this);
}


void TrackRun::execute() {
    const int currentVersion = ((OPAL_VERSION_MAJOR * 100) + OPAL_VERSION_MINOR) * 100;
    if (Options::version < currentVersion) {
        unsigned int fileVersion = Options::version / 100;
        bool newerChanges = false;
        for (auto it = Versions::changes.begin(); it != Versions::changes.end(); ++ it) {
            if (it->first > fileVersion) {
                newerChanges = true;
                break;
            }
        }
        if (newerChanges) {
            Inform errorMsg("Error");
            errorMsg << "\n******************** V E R S I O N   M I S M A T C H ***********************\n" << endl;
            for (auto it = Versions::changes.begin(); it != Versions::changes.end(); ++ it) {
                if (it->first > fileVersion) {
                    errorMsg << it->second << endl;
                }
            }
            errorMsg << "\n"
                     << "* Make sure you do understand these changes and adjust your input file \n"
                     << "* accordingly. Then add\n"
                     << "* OPTION, VERSION = " << currentVersion << ";\n"
                     << "* to your input file. " << endl;
            errorMsg << "\n****************************************************************************\n" << endl;
            throw OpalException("TrackRun::execute", "Version mismatch");
        }
    }

    isFollowupTrack_m = opal->hasBunchAllocated();
    if (!itsAttr[DISTRIBUTION] && !isFollowupTrack_m) {
        throw OpalException("TrackRun::execute",
                            "\"DISTRIBUTION\" must be set in \"RUN\" command.");
    }
    if (!itsAttr[FIELDSOLVER]) {
        throw OpalException("TrackRun::execute",
                            "\"FIELDSOLVER\" must be set in \"RUN\" command.");
    }
    if (!itsAttr[BEAM]) {
        throw OpalException("TrackRun::execute",
                            "\"BEAM\" must be set in \"RUN\" command.");
    }

    // Get algorithm to use.
    setRunMethod();
    switch (method_m) {
       case RunMethod::THICK: {
            setupThickTracker();
            break;
        }
        case RunMethod::PARALLELT: {
            setupTTracker();
            break;
        }
        case RunMethod::CYCLOTRONT: {
            setupCyclotronTracker();
            break;
        }
        default: {
            throw OpalException("TrackRun::execute",
                                "Unknown \"METHOD\" for the \"RUN\" command");
        }
    }

    if (method_m == RunMethod::THICK) {
       	int turns = int(std::round(Attributes::getReal(itsAttr[TURNS])));

        // Track for the all but last turn.
        for (int turn = 1; turn < turns; ++turn) {
            itsTracker->execute();
        }
        // Track the last turn.
        itsTracker->execute();

    } else {
        itsTracker->execute();

        opal->setRestartRun(false);
    }

    opal->bunchIsAllocated();
}

void TrackRun::setRunMethod() {
    if (!itsAttr[METHOD]) {
        throw OpalException("TrackRun::setRunMethod",
                            "The attribute \"METHOD\" isn't set for the \"RUN\" command");
    } else {
        auto it = stringMethod_s.right.find(Attributes::getString(itsAttr[METHOD]));
        if (it != stringMethod_s.right.end()) {
            method_m = it->second;
        }
    }
}

std::string TrackRun::getRunMethodName() const {
    return stringMethod_s.left.at(method_m);
}

void TrackRun::setupThickTracker() {
    if (isFollowupTrack_m) {
        Track::block->bunch->setLocalTrackStep(0);
    }

    Beam* beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

    setBoundaryGeometry();

    setupFieldsolver();

    if (opal->inRestartRun()) {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  opal->getRestartStep(),
                                                  OpalData::getInstance()->getRestartFileName(),
                                                  H5_O_WRONLY);
    } else if (isFollowupTrack_m) {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  -1,
                                                  opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    } else {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    }

    macrocharge_m = setDistributionParallelT(beam);

    *gmsg << *this  << endl;

    Track::block->bunch->setdT(Track::block->dT.front());
    Track::block->bunch->dtScInit_m = Track::block->dtScInit;
    Track::block->bunch->deltaTau_m = Track::block->deltaTau;

    if (!isFollowupTrack_m && !opal->inRestartRun()) {
        Track::block->bunch->setT(Track::block->t0_m);
    }

    if (Track::block->bunch->getIfBeamEmitting()) {
        Track::block->bunch->setChargeZeroPart(macrocharge_m);
    } else {
        Track::block->bunch->setCharge(macrocharge_m);
    }

    // set coupling constant
    double coefE = 1.0 / (4 * Physics::pi * Physics::epsilon_0);
    Track::block->bunch->setCouplingConstant(coefE);

    // statistical data are calculated (rms, eps etc.)
    Track::block->bunch->calcBeamParameters();

    initDataSink();

    if (!isFollowupTrack_m) {
        *gmsg << *dist << endl;
    }

    if (Track::block->bunch->getTotalNum() > 0) {
        double spos = /*Track::block->bunch->get_sPos() +*/ Track::block->zstart;
        auto &zstop = Track::block->zstop;
        auto &timeStep = Track::block->localTimeSteps;
        auto &dT = Track::block->dT;

        unsigned int i = 0;
        while (i + 1 < zstop.size() && zstop[i + 1] < spos) {
            ++ i;
        }

        zstop.erase(zstop.begin(), zstop.begin() + i);
        timeStep.erase(timeStep.begin(), timeStep.begin() + i);
        dT.erase(dT.begin(), dT.begin() + i);

        Track::block->bunch->setdT(dT.front());
    } else {
        Track::block->zstart = 0.0;
    }

    *gmsg << *beam << endl;
    *gmsg << *fs   << endl;

    itsTracker.reset(new ThickTracker(*Track::block->use->fetchLine(),
                                  Track::block->bunch, *beam, *ds, Track::block->reference,
                                  false, false, Track::block->localTimeSteps,
                                  Track::block->zstart, Track::block->zstop, Track::block->dT,
                                  Track::block->truncOrder));
}


void TrackRun::setupTTracker(){
    OpalData::getInstance()->setInOPALTMode();

    if (isFollowupTrack_m) {
        Track::block->bunch->setLocalTrackStep(0);
    }

    Beam* beam = Beam::find(Attributes::getString(itsAttr[BEAM]));
    Track::block->bunch->setBeamFrequency(beam->getFrequency() * Units::MHz2Hz);
    Track::block->bunch->setPType(beam->getParticleName());

    setBoundaryGeometry();

    setupFieldsolver();

    if (opal->inRestartRun()) {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  opal->getRestartStep(),
                                                  OpalData::getInstance()->getRestartFileName(),
                                                  H5_O_WRONLY);
    } else if (isFollowupTrack_m) {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  -1,
                                                  opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    } else {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    }

    macrocharge_m = setDistributionParallelT(beam);
    macromass_m   = beam->getMassPerParticle();

    *gmsg << *this  << endl;

    Track::block->bunch->setdT(Track::block->dT.front());
    Track::block->bunch->dtScInit_m = Track::block->dtScInit;
    Track::block->bunch->deltaTau_m = Track::block->deltaTau;

    if (!isFollowupTrack_m && !opal->inRestartRun()) {
        Track::block->bunch->setT(Track::block->t0_m);
    }

    if (Track::block->bunch->getIfBeamEmitting()) {
        Track::block->bunch->setChargeZeroPart(macrocharge_m);
        Track::block->bunch->setMassZeroPart(macromass_m);
    } else {
        Track::block->bunch->setCharge(macrocharge_m);
        Track::block->bunch->setMass(macromass_m);
    }
    // set coupling constant
    double coefE = 1.0 / (4 * Physics::pi * Physics::epsilon_0);
    Track::block->bunch->setCouplingConstant(coefE);

    // statistical data are calculated (rms, eps etc.)
    Track::block->bunch->calcBeamParameters();

    initDataSink();

    if (!isFollowupTrack_m) {
        *gmsg << std::scientific;
        *gmsg << *dist << endl;
    }

    if (Track::block->bunch->getTotalNum() > 0) {
        double spos = Track::block->zstart;
        auto& zstop = Track::block->zstop;
        auto it = Track::block->dT.begin();

        unsigned int i = 0;
        while (i + 1 < zstop.size() && zstop[i + 1] < spos) {
            ++ i;
            ++ it;
        }

        Track::block->bunch->setdT(*it);
    } else {
        Track::block->zstart = 0.0;
    }

    *gmsg << *beam << endl;
    *gmsg << *fs   << endl;

    // findPhasesForMaxEnergy();
    
    itsTracker.reset(new ParallelTTracker(*Track::block->use->fetchLine(),
                                      Track::block->bunch,
                                      *ds,
                                      Track::block->reference,
                                      false,
                                      Attributes::getBool(itsAttr[TRACKBACK]),
                                      Track::block->localTimeSteps,
                                      Track::block->zstart,
                                      Track::block->zstop,
                                      Track::block->dT));
}

void TrackRun::setupCyclotronTracker(){

    OpalData::getInstance()->setInOPALCyclMode();
    Beam* beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

    setBoundaryGeometry();

    setupFieldsolver();

    Track::block->bunch->setPType(beam->getParticleName());
    Track::block->bunch->POrigin = ParticleOrigin::REGULAR;

    std::vector<std::string> distr_str = Attributes::getStringArray(itsAttr[DISTRIBUTION]);
    if (distr_str.size() == 0) {
        dist = Distribution::find(defaultDistribution);
    } else {
        dist = Distribution::find(distr_str.at(0));
    }

    // multi-bunch parameters
    const int specifiedNumBunch = int(std::abs(std::round(Attributes::getReal(itsAttr[TURNS]))));
    const double mbPara         = Attributes::getReal(itsAttr[PARAMB]);
    const std::string mbMode    = Attributes::getString(itsAttr[MBMODE]);
    const double mbEta          = Attributes::getReal(itsAttr[MB_ETA]);
    const std::string mbBinning = Attributes::getString(itsAttr[MB_BINNING]);

    if (opal->inRestartRun()) {
        phaseSpaceSink_m = new H5PartWrapperForPC(opal->getInputBasename() + std::string(".h5"),
                                                  opal->getRestartStep(),
                                                  OpalData::getInstance()->getRestartFileName(),
                                                  H5_O_WRONLY);
    } else if (isFollowupTrack_m) {
        phaseSpaceSink_m = new H5PartWrapperForPC(opal->getInputBasename() + std::string(".h5"),
                                                  -1,
                                                  opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    } else {
        phaseSpaceSink_m = new H5PartWrapperForPC(opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    }

    if (beam->getNumberOfParticles() < 3 || beam->getCurrent() == 0.0) {
        macrocharge_m = beam->getCharge() * Physics::q_e;
        macromass_m   = beam->getMass();
        Track::block->bunch->setDistribution(dist,
                                             beam->getNumberOfParticles(),
                                             beam->getCurrent(),
                                             *Track::block->use->fetchLine());

    } else {
        /**
           getFrequency() gets RF frequency [MHz], NOT isochronous revolution frequency of particle!
           getCurrent() gets beamcurrent [A]
        */
        macrocharge_m = beam->getChargePerParticle();
        macromass_m   = beam->getMassPerParticle();

        if (!isFollowupTrack_m) {
            if (!opal->inRestartRun()) {
                Track::block->bunch->setDistribution(dist,
                                                     beam->getNumberOfParticles(),
                                                     beam->getCurrent(),
                                                     *Track::block->use->fetchLine());

            } else {
                dist->doRestartOpalCycl(Track::block->bunch,
                                        beam->getNumberOfParticles(),
                                        opal->getRestartStep(),
                                        specifiedNumBunch,
                                        phaseSpaceSink_m);
            }
        }
    }
    Track::block->bunch->setMass(macromass_m); // set the Mass per macro-particle, [GeV/c^2]
    Track::block->bunch->setCharge(macrocharge_m);  // set the charge per macro-particle, [C]

    Track::block->bunch->setdT(1.0 / (Track::block->stepsPerTurn * beam->getFrequency() * Units::MHz2Hz));
    Track::block->bunch->setStepsPerTurn(Track::block->stepsPerTurn);

    // set coupling constant
    double coefE = 1.0 / (4 * Physics::pi * Physics::epsilon_0);
    Track::block->bunch->setCouplingConstant(coefE);

    // statistical data are calculated (rms, eps etc.)
    Track::block->bunch->calcBeamParameters();

    initDataSink(specifiedNumBunch);

    itsTracker.reset(new ParallelCyclotronTracker(*Track::block->use->fetchLine(),
                                              Track::block->bunch, *ds, Track::block->reference,
                                              false, false, Track::block->localTimeSteps.front(),
                                              Track::block->timeIntegrator,
                                              specifiedNumBunch, mbEta, mbPara, mbMode, mbBinning));

    ParallelCyclotronTracker* cyclTracker = dynamic_cast<ParallelCyclotronTracker*>(itsTracker.get());

    if (opal->inRestartRun()) {
        H5PartWrapperForPC *h5pw = static_cast<H5PartWrapperForPC*>(phaseSpaceSink_m);
        cyclTracker->setBeGa(h5pw->getMeanMomentum());

        cyclTracker->setPr(h5pw->getReferencePr());
        cyclTracker->setPt(h5pw->getReferencePt());
        cyclTracker->setPz(h5pw->getReferencePz());

        cyclTracker->setR(h5pw->getReferenceR());
        cyclTracker->setTheta(h5pw->getReferenceT());
        cyclTracker->setZ(h5pw->getReferenceZ());

        // The following is for restarts in local frame
        cyclTracker->setPhi(h5pw->getAzimuth());
        cyclTracker->setPsi(h5pw->getElevation());
        cyclTracker->setPreviousH5Local(h5pw->getPreviousH5Local());

        if ( specifiedNumBunch > 1 ) {
            cyclTracker->setLastDumpedStep(opal->getRestartStep());
        }
    }

    // statistical data are calculated (rms, eps etc.)
    Track::block->bunch->calcBeamParameters();

    *gmsg << *this  << endl;
    *gmsg << *dist << endl;
    *gmsg << *beam << endl;
    *gmsg << *fs   << endl;
}

void TrackRun::setupFieldsolver() {
    fs = FieldSolver::find(Attributes::getString(itsAttr[FIELDSOLVER]));

    if (fs->getFieldSolverType() != FieldSolverType::NONE) {
        size_t numGridPoints = fs->getMX()*fs->getMY()*fs->getMT(); // total number of gridpoints
        Beam* beam = Beam::find(Attributes::getString(itsAttr[BEAM]));
        size_t numParticles = beam->getNumberOfParticles();

        if (!opal->inRestartRun() && numParticles < numGridPoints
            && fs->getFieldSolverType() != FieldSolverType::SAAMG // in SPIRAL/SAAMG we're meshing the whole domain -DW
            && fs->getFieldSolverType() != FieldSolverType::P3M //In P3M with one-one mapping grid points can be less than particles
            && !Options::amr)
        {
            throw OpalException("TrackRun::setupFieldsolver()",
                                "The number of simulation particles (" + std::to_string(numParticles) + ") \n" +
                                "is smaller than the number of gridpoints (" + std::to_string(numGridPoints) + ").\n" +
                                "Please increase the number of particles or reduce the size of the mesh.\n");
        }

        OpalData::getInstance()->addProblemCharacteristicValue("MX", fs->getMX());
        OpalData::getInstance()->addProblemCharacteristicValue("MY", fs->getMY());
        OpalData::getInstance()->addProblemCharacteristicValue("MT", fs->getMT());
    }

    fs->initCartesianFields();
    Track::block->bunch->setSolver(fs);
    if (fs->hasPeriodicZ()) {
        Track::block->bunch->setBCForDCBeam();
    } else {
        Track::block->bunch->setBCAllOpen();
    }
}


void TrackRun::initDataSink(const int& numBunch) {
    if (!opal->inRestartRun()) {
        if (!opal->hasDataSinkAllocated()) {
            opal->setDataSink(new DataSink(phaseSpaceSink_m, false, numBunch));
        } else {
            ds = opal->getDataSink();
            ds->changeH5Wrapper(phaseSpaceSink_m);
        }
    } else {
        opal->setDataSink(new DataSink(phaseSpaceSink_m, true, numBunch));
    }
    ds = opal->getDataSink();
}

void TrackRun::setBoundaryGeometry() {
    if (Attributes::getString(itsAttr[BOUNDARYGEOMETRY]) != "NONE") {
        // Ask the dictionary if BoundaryGeometry is allocated.
        // If it is allocated use the allocated BoundaryGeometry
        if (!OpalData::getInstance()->hasGlobalGeometry()) {
            const std::string geomDescriptor = Attributes::getString(itsAttr[BOUNDARYGEOMETRY]);
            BoundaryGeometry* bg = BoundaryGeometry::find(geomDescriptor)->clone(geomDescriptor);
            OpalData::getInstance()->setGlobalGeometry(bg);
        }
    }
}


double TrackRun::setDistributionParallelT(Beam* beam) {
    /*
     * Distribution(s) can be set via a single distribution or a list
     * (array) of distributions. If an array is defined the first in the
     * list is treated as the primary distribution. All others are added to
     * it to create the full distribution.
     */
    std::vector<std::string> distributionArray
        = Attributes::getStringArray(itsAttr[DISTRIBUTION]);
    const size_t numberOfDistributions = distributionArray.size();

    if (numberOfDistributions == 0) {
        dist = Distribution::find(defaultDistribution);
    } else {
        dist = Distribution::find(distributionArray.at(0));
        dist->setNumberOfDistributions(numberOfDistributions);

        if (numberOfDistributions > 1) {
            *gmsg << endl
                  << "---------------------------------" << endl
                  << "Found more than one distribution:" << endl << endl;
            *gmsg << "Main Distribution" << endl
                  << "---------------------------------" << endl
                  << distributionArray.at(0) << endl << endl
                  << "Secondary Distribution(s)" << endl
                  << "---------------------------------" << endl;

            for (size_t i = 1; i < numberOfDistributions; ++ i) {
                Distribution *distribution = Distribution::find(distributionArray.at(i));
                distribution->setNumberOfDistributions(numberOfDistributions);
                distrs_m.push_back(distribution);

                *gmsg << distributionArray.at(i) << endl;
            }
            *gmsg << endl
                  << "---------------------------------" << endl << endl;
        }
    }

    /*
     * Initialize distributions.
     */
    size_t numberOfParticles = beam->getNumberOfParticles();
    if (!isFollowupTrack_m) {
        if (!opal->inRestartRun()) {
            /*
             * Here we are not doing a restart run
             * and we do not have a bunch already allocated.
             */
            Track::block->bunch->setDistribution(dist,
                                                 distrs_m,
                                                 numberOfParticles);

            /*
             * If this is an injected beam (rather than an emitted beam), we
             * make sure it doesn't have any particles at z < 0.
             */

            opal->setGlobalPhaseShift(0.5 * dist->getTEmission() + dist->getEmissionTimeShift());
        } else {
            /*
             * Read in beam from restart file.
             */
            dist->doRestartOpalT(Track::block->bunch, numberOfParticles, opal->getRestartStep(), phaseSpaceSink_m);
        }
    }

    // Return charge per macroparticle.
    return beam->getChargePerParticle();
}

Inform& TrackRun::print(Inform& os) const {
    os << endl;
    os << "* ************* T R A C K  R U N *************************************************** " << endl;
    if (!isFollowupTrack_m) {
        os << "* Selected Tracking Method == " << getRunMethodName() << ", NEW TRACK" << '\n'
           << "* ********************************************************************************** " << '\n';
    } else {
        os << "* Selected Tracking Method == " << getRunMethodName() << ", FOLLOWUP TRACK" << '\n'
           << "* ********************************************************************************** " << '\n';
    }
    os << "* Phase space dump frequency    = " << Options::psDumpFreq << '\n'
       << "* Statistics dump frequency     = " << Options::statDumpFreq << " w.r.t. the time step." << '\n'
       << "* DT                            = " << Track::block->dT.front() << " [s]\n"
       << "* MAXSTEPS                      = " << Track::block->localTimeSteps.front() << '\n'
       << "* Mass of simulation particle   = " << macromass_m << " [GeV/c^2]" << '\n'
       << "* Charge of simulation particle = " << macrocharge_m << " [C]" << '\n';
    if (method_m == RunMethod::CYCLOTRONT) {
        os << "* Number of neighbour bunches   = " << int(std::abs(std::round(Attributes::getReal(itsAttr[TURNS])))) << '\n'
           << "* STEPSPERTURN                  = " << Track::block->stepsPerTurn << '\n';
    }
    os << "* ********************************************************************************** ";
    return os;
}

std::shared_ptr<Tracker> TrackRun::getTracker() {
    return itsTracker;
}
