//
// Class ParallelCyclotronTracker
//   Tracker for OPAL-Cycl
//
// Copyright (c) 2007 - 2014, Jianjun Yang, Andreas Adelmann and Matthias Toggweiler,
//                            Paul Scherrer Institut, Villigen PSI, Switzerland
// Copyright (c) 2014,        Daniel Winklehner, MIT, Cambridge, MA, USA
// Copyright (c) 2012 - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects"
// and the paper
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects:
// Model, implementation, and application"
// (https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.13.064201)
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
#include "Algorithms/ParallelCyclotronTracker.h"

#include "AbsBeamline/CCollimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Cyclotron.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/MultipoleT.h"
#include "AbsBeamline/MultipoleTBase.h"
#include "AbsBeamline/MultipoleTCurvedConstRadius.h"
#include "AbsBeamline/MultipoleTCurvedVarRadius.h"
#include "AbsBeamline/MultipoleTStraight.h"
#include "AbsBeamline/Offset.h"
#include "AbsBeamline/PluginElement.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/Ring.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/SBend3D.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/Stripper.h"
#include "AbsBeamline/Vacuum.h"
#include "AbsBeamline/VariableRFCavity.h"
#include "AbsBeamline/VariableRFCavityFringeField.h"
#include "AbsBeamline/VerticalFFAMagnet.h"

#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"

#include "Algorithms/Ctunes.h"
#include "Algorithms/PolynomialTimeDependence.h"

#include "BasicActions/DumpEMFields.h"
#include "BasicActions/DumpFields.h"

#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedBeamline.h"

#include "Distribution/Distribution.h"

#include "Elements/OpalBeamline.h"

#include "Physics/Physics.h"
#include "Physics/Units.h"

#include "Structure/BoundaryGeometry.h"
#include "Structure/DataSink.h"
#include "Structure/LossDataSink.h"

#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <boost/filesystem.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>

constexpr double c_mmtns = Physics::c * Units::m2mm / Units::s2ns;

Vector_t const ParallelCyclotronTracker::xaxis = Vector_t(1.0, 0.0, 0.0);
Vector_t const ParallelCyclotronTracker::yaxis = Vector_t(0.0, 1.0, 0.0);
Vector_t const ParallelCyclotronTracker::zaxis = Vector_t(0.0, 0.0, 1.0);

extern Inform *gmsg;

/**
 * Constructor ParallelCyclotronTracker
 *
 * @param beamline
 * @param bunch
 * @param ds
 * @param reference
 * @param revBeam
 * @param revTrack
 * @param maxSTEPS
 * @param timeIntegrator
 */
ParallelCyclotronTracker::ParallelCyclotronTracker(const Beamline& beamline,
                                                   PartBunchBase<double, 3>* bunch,
                                                   DataSink& ds,
                                                   const PartData& reference,
                                                   bool revBeam, bool revTrack,
                                                   int maxSTEPS,
                                                   Steppers::TimeIntegrator timeintegrator,
                                                   const int& numBunch,
                                                   const double& mbEta,
                                                   const double& mbPara,
                                                   const std::string& mbMode,
                                                   const std::string& mbBinning)
    : Tracker(beamline, bunch, reference, revBeam, revTrack)
    , bgf_m(nullptr)
    , maxSteps_m(maxSTEPS)
    , lastDumpedStep_m(0)
    , myNode_m(Ippl::myNode())
    , initialLocalNum_m(bunch->getLocalNum())
    , initialTotalNum_m(bunch->getTotalNum())
    , opalRing_m(nullptr)
    , itsStepper_mp(nullptr)
    , mode_m(TrackingMode::UNDEFINED)
    , stepper_m(timeintegrator)
{
    itsBeamline = dynamic_cast<Beamline*>(beamline.clone());
    itsDataSink = &ds;

    if ( numBunch > 1 ) {
        mbHandler_m = std::unique_ptr<MultiBunchHandler>(
            new MultiBunchHandler(bunch, numBunch, mbEta,
                                  mbPara, mbMode, mbBinning)
        );
    }

    IntegrationTimer_m = IpplTimings::getTimer("Integration");
    TransformTimer_m   = IpplTimings::getTimer("Frametransform");
    DumpTimer_m        = IpplTimings::getTimer("Dump");
    BinRepartTimer_m   = IpplTimings::getTimer("Binaryrepart");
    PluginElemTimer_m  = IpplTimings::getTimer("PluginElements");
    DelParticleTimer_m = IpplTimings::getTimer("DeleteParticles");

    setTrackingMode();
}

/**
 * Destructor ParallelCyclotronTracker
 *
 */
ParallelCyclotronTracker::~ParallelCyclotronTracker() {
    if (bgf_m)
        lossDs_m->save();
    for (Component* component : myElements) {
        delete(component);
    }
    for (auto fd : FieldDimensions) {
        delete(fd);
    }
    delete itsBeamline;
    // delete opalRing_m;
}


void ParallelCyclotronTracker::bgf_main_collision_test() {
    if (!bgf_m) return;

    Inform msg("bgf_main_collision_test ");

    /**
     *Here we check if a particle is outside the domain, flag it for deletion
     */

    Vector_t intecoords = 0.0;

    // This has to match the dT in the rk4 pusher
    double dtime = itsBunch_m->getdT() * getHarmonicNumber();

    int triId = 0;
    for (size_t i = 0; i < itsBunch_m->getLocalNum(); i++) {
        int res = bgf_m->partInside(itsBunch_m->R[i], itsBunch_m->P[i],
                                    dtime, intecoords, triId);
        if (res >= 0) {
            lossDs_m->addParticle(OpalParticle(itsBunch_m->ID[i],
                                               itsBunch_m->R[i], itsBunch_m->P[i],
                                               itsBunch_m->getT(),
                                               itsBunch_m->Q[i], itsBunch_m->M[i]),
                                  std::make_pair(turnnumber_m, itsBunch_m->bunchNum[i]));
            itsBunch_m->Bin[i] = -1;
            Inform gmsgALL("OPAL", INFORM_ALL_NODES);
            gmsgALL << level4 << "* Particle " << itsBunch_m->ID[i]
                    << " lost on boundary geometry" << endl;
        }
    }
}

// only used for dumping into stat file
void ParallelCyclotronTracker::dumpAngle(const double& theta,
                                         double& prevAzimuth,
                                         double& azimuth,
                                         const short& bunchNr) {
    if ( prevAzimuth < 0.0 ) { // only at first occurrence
        double plus = 0.0;
        if ( OpalData::getInstance()->inRestartRun() ) {
            plus = 360.0 * (turnnumber_m - bunchNr);
        }
        azimuth = theta + plus;
    } else {
        double dtheta = theta - prevAzimuth;
        if ( dtheta < 0 ) {
            dtheta += 360.0;
        }
        if ( dtheta > 180 ) { // rotating clockwise, reduce angle
            dtheta -= 360;
        }
        azimuth += dtheta;
    }
    prevAzimuth = theta;
}


double ParallelCyclotronTracker::computeRadius(const Vector_t& meanR) const {
    return Units::m2mm * std::sqrt(meanR(0) * meanR(0) + meanR(1) * meanR(1));
}


void ParallelCyclotronTracker::computePathLengthUpdate(std::vector<double>& dl,
                                                       const double& dt)
{
    // the last element in dotP is the dot-product over all particles
    std::vector<double> dotP(dl.size());
    if ( Options::psDumpFrame == DumpFrame::BUNCH_MEAN || isMultiBunch()) {

        for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {
            dotP[itsBunch_m->bunchNum[i]] += dot(itsBunch_m->P[i], itsBunch_m->P[i]);
        }

        allreduce(dotP.data(), dotP.size(), std::plus<double>());

        // dot-product over all particles
        double sum = std::accumulate(dotP.begin(), dotP.end() - 1, 0);
        dotP.back() = sum / double(itsBunch_m->getTotalNum());

        // bunch specific --> multi-bunches only
        for (short b = 0; b < (short)dotP.size() - 1; ++b) {
            dotP[b] /= double(itsBunch_m->getTotalNumPerBunch(b));
        }

    } else if ( itsBunch_m->getLocalNum() == 0 ) {
        // here we are in DumpFrame::GLOBAL mode
        dotP[0] = 0.0;
    } else {
        // here we are in DumpFrame::GLOBAL mode
        dotP[0] = dot(itsBunch_m->P[0], itsBunch_m->P[0]);
    }

    for (size_t i = 0; i < dotP.size(); ++i) {
        double const gamma = std::sqrt(1.0 + dotP[i]);
        double const beta  = std::sqrt(dotP[i]) / gamma;
        dl[i] = c_mmtns * dt * Units::mm2m * beta;

    }
}


/**
 *
 *
 * @param fn Base file name
 */
void ParallelCyclotronTracker::openFiles(size_t numFiles, std::string SfileName) {

    for (unsigned int i=0; i<numFiles; i++) {
        std::string SfileName2 = SfileName;
        if (i<=2) {
            SfileName2 += std::string("-Angle" + std::to_string(int(i)) + ".dat");
        }
        else {
            // for single Particle Mode, output after each turn, to define matched initial phase ellipse.
            SfileName2 += std::string("-afterEachTurn.dat");
        }

        outfTheta_m.emplace_back(new std::ofstream(SfileName2.c_str()));
        outfTheta_m.back()->precision(8);
        outfTheta_m.back()->setf(std::ios::scientific, std::ios::floatfield);
        *outfTheta_m.back() << "# r [mm]        beta_r*gamma       "
                            << "theta [deg]     beta_theta*gamma        "
                            << "z [mm]          beta_z*gamma" << std::endl;
    }
}

/**
 * Close all files related to
 * special output in the Cyclotron
 * mode.
 */
void ParallelCyclotronTracker::closeFiles() {
    for (auto & file : outfTheta_m) {
        file->close();
    }
}


/**
 *
 *
 * @param cycl
 */
void ParallelCyclotronTracker::visitCyclotron(const Cyclotron& cycl) {
    *gmsg << "* ----------------------------- Cyclotron -------------------------------- *" << endl;

    cycl_m = dynamic_cast<Cyclotron*>(cycl.clone());
    myElements.push_back(cycl_m);

    // Is this a Spiral Inflector Simulation? If yes, we'll give the user some
    // useful information
    spiral_flag = cycl_m->getSpiralFlag();
    if (spiral_flag) {
        *gmsg << endl << "* This is a Spiral Inflector Simulation! This means the following:" << endl;
        *gmsg         << "* 1.) It is up to the user to provide appropriate geometry, electric and magnetic fields!" << endl;
        *gmsg         << "*     (Use BANDRF type cyclotron and use RFMAPFN to load both magnetic" << endl;
        *gmsg         << "*     and electric fields, setting SUPERPOSE to an array of TRUE values.)" << endl;
        *gmsg         << "* 2.) For high currents it is strongly recommended to use the SAAMG fieldsolver," << endl;
        *gmsg         << "*     FFT does not give the correct results (boundary conditions are missing)." << endl;
        *gmsg         << "* 3.) The whole geometry will be meshed and used for the fieldsolver." << endl;
        *gmsg         << "*     There will be no transformations of the bunch into a local frame und consequently," << endl;
        *gmsg         << "*     the problem will be treated non-relativistically!" << endl;
        *gmsg         << "*     (This is not an issue for spiral inflectors as they are typically < 100 keV/amu.)" << endl;
        *gmsg << endl << "* Note: For now, multi-bunch mode (MBM) needs to be de-activated for spiral inflector" << endl;
        *gmsg         << "* and space charge needs to be solved every time-step. numBunch_m and scSolveFreq are reset." << endl;
        if (isMultiBunch()) {
            mbHandler_m = nullptr;
        }
    }

    // Fresh run (no restart):
    if (!OpalData::getInstance()->inRestartRun()) {

        // Get reference values from cyclotron element
        // For now, these are still stored in mm. should be the only ones. -DW
        referenceR     = cycl_m->getRinit();
        referenceTheta = cycl_m->getPHIinit();
        referenceZ     = cycl_m->getZinit();
        referencePr    = cycl_m->getPRinit();
        referencePz    = cycl_m->getPZinit();

        if (referenceTheta <= -180.0 || referenceTheta > 180.0) {
            throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron",
                                "PHIINIT is out of [-180, 180)!");
        }

        referencePtot =  itsReference.getGamma() * itsReference.getBeta();

        // Calculate reference azimuthal (tangential) momentum from total-, z- and radial momentum:
        float insqrt = referencePtot * referencePtot - \
            referencePr * referencePr - referencePz * referencePz;

        if (insqrt < 0) {
            if (insqrt > -1.0e-10) {
                referencePt = 0.0;
            } else {
                throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron",
                                    "Pt imaginary!");
            }
        } else {
            referencePt = std::sqrt(insqrt);
        }

        if (referencePtot < 0.0) {
            referencePt *= -1.0;
        }
        // End calculate referencePt

        // Restart a run:
    } else {

        // If the user wants to save the restarted run in local frame,
        // make sure the previous h5 file was local too
        if (Options::psDumpFrame != DumpFrame::GLOBAL) {
            if (!previousH5Local) {
                throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron",
                                    "You are trying a local restart from a global h5 file!");
            }
            // Else, if the user wants to save the restarted run in global frame,
            // make sure the previous h5 file was global too
        } else {
            if (previousH5Local) {
                throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron",
                                    "You are trying a global restart from a local h5 file!");
            }
        }

        // Adjust some of the reference variables from the h5 file
        referencePhi *= Units::deg2rad;
        referencePsi *= Units::deg2rad;
        referencePtot = bega;
        if (referenceTheta <= -180.0 || referenceTheta > 180.0) {
            throw OpalException("Error in ParallelCyclotronTracker::visitCyclotron",
                                "PHIINIT is out of [-180, 180)!");
        }
    }

    sinRefTheta_m = std::sin(referenceTheta * Units::deg2rad);
    cosRefTheta_m = std::cos(referenceTheta * Units::deg2rad);

    *gmsg << endl;
    *gmsg << "* Bunch global starting position:" << endl;
    *gmsg << "* RINIT   = " << referenceR  << " [mm]" << endl;
    *gmsg << "* PHIINIT = " << referenceTheta << " [deg]" << endl;
    *gmsg << "* ZINIT   = " << referenceZ << " [mm]" << endl;
    *gmsg << endl;
    *gmsg << "* Bunch global starting momenta:" << endl;
    *gmsg << "* Initial gamma = " << itsReference.getGamma() << endl;
    *gmsg << "* Initial beta  = " << itsReference.getBeta() << endl;
    *gmsg << "* Reference total momentum          = " << referencePtot << " [beta gamma]" << endl;
    *gmsg << "* Reference azimuthal momentum (Pt) = " << referencePt << " [beta gamma]" << endl;
    *gmsg << "* Reference radial momentum (Pr)    = " << referencePr  << " [beta gamma]" << endl;
    *gmsg << "* Reference axial momentum (Pz)     = " << referencePz << " [beta gamma]" << endl;
    *gmsg << endl;

    double sym = cycl_m->getSymmetry();
    *gmsg << "* " << sym << "-fold field symmetry " << endl;

    // ckr: this just returned the default value as defined in Component.h
    // double rff = cycl_m->getRfFrequ(0);
    // *gmsg << "* Rf frequency= " << rff << " [MHz]" << endl;

    std::string fmfn = cycl_m->getFieldMapFN();
    *gmsg << "* Field map file       = '" << fmfn << "'" << endl;

    std::string type = cycl_m->getCyclotronType();
    *gmsg << "* Type of cyclotron    = " << type << " " << endl;

    double rmin = cycl_m->getMinR();
    double rmax = cycl_m->getMaxR();
    *gmsg << "* Radial aperture      = " << rmin << " ... " << rmax<<" [m] "<< endl;

    double zmin = cycl_m->getMinZ();
    double zmax = cycl_m->getMaxZ();
    *gmsg << "* Vertical aperture    = " << zmin << " ... " << zmax<<" [m]"<< endl;

    double h = cycl_m->getCyclHarm();
    *gmsg << "* Number of trimcoils  = " << cycl_m->getNumberOfTrimcoils() << endl;
    *gmsg << "* Harmonic number h    = " << h << " " << endl;

    std::vector<double> rffrequ = cycl_m->getRfFrequ();
    *gmsg << "* RF frequency         = " << Util::doubleVectorToString(rffrequ) << " [MHz]" << endl;

    cycl_m->setBFieldType();
    if (cycl_m->getBFieldType() == Cyclotron::BFieldType::BANDRF) {
        std::vector<double> rfphi = cycl_m->getRfPhi();
        for (size_t i = 0; i < rfphi.size(); ++i) {
            rfphi[i] = rfphi[i] * Units::rad2deg;
        }
        *gmsg << "* RF inital phase      = " << Util::doubleVectorToString(rfphi) << " [deg]" << endl;

        std::vector<double> escale = cycl_m->getEScale();
        *gmsg << "* E-field scale factor = " << Util::doubleVectorToString(escale) << endl;

        std::vector<bool> superpose = cycl_m->getSuperpose();
        *gmsg << "* Superpose electric field maps -> " << Util::boolVectorToUpperString(superpose) << endl;
    }

    // Read in cyclotron field maps
    cycl_m->initialise(itsBunch_m, cycl_m->getBScale());

    double BcParameter[8] = {};
    BcParameter[0] = Units::mm2m * cycl_m->getRmin();
    BcParameter[1] = Units::mm2m * cycl_m->getRmax();

    // Store inner radius and outer radius of cyclotron field map in the list
    buildupFieldList(BcParameter, ElementType::CYCLOTRON, cycl_m);
}

/**
 *
 *
 * @param coll
 */
void ParallelCyclotronTracker::visitCCollimator(const CCollimator& coll) {
    *gmsg << "* ----------------------------- Collimator ------------------------------- *" << endl;

    CCollimator* elptr = dynamic_cast<CCollimator*>(coll.clone());
    myElements.push_back(elptr);

    *gmsg << "* Name    = " << elptr->getName() << endl;

    double xstart = elptr->getXStart();
    *gmsg << "* XStart  = " << xstart << " [m]" << endl;

    double xend = elptr->getXEnd();
    *gmsg << "* XEnd    = " << xend << " [m]" << endl;

    double ystart = elptr->getYStart();
    *gmsg << "* YStart  = " << ystart << " [m]" << endl;

    double yend = elptr->getYEnd();
    *gmsg << "* YEnd    = " << yend << " [m]" << endl;

    double zstart = elptr->getZStart();
    *gmsg << "* ZStart  = " << zstart << " [m]" << endl;

    double zend = elptr->getZEnd();
    *gmsg << "* ZEnd    = " << zend << " [m]" << endl;

    double width = elptr->getWidth();
    *gmsg << "* Width   = " << width << " [m]" << endl;

    elptr->initialise(itsBunch_m);

    double BcParameter[8] = {};

    BcParameter[0] = xstart;
    BcParameter[1] = xend;
    BcParameter[2] = ystart;
    BcParameter[3] = yend;
    BcParameter[4] = width;

    buildupFieldList(BcParameter, ElementType::CCOLLIMATOR, elptr);
}

/**
 *
 *
 * @param corr
 */
void ParallelCyclotronTracker::visitCorrector(const Corrector& corr) {
    *gmsg << "In Corrector; L= " << corr.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Corrector*>(corr.clone()));
}

/**
 *
 *
 * @param degrader
 */
void ParallelCyclotronTracker::visitDegrader(const Degrader& deg) {
    *gmsg << "In Degrader; L= " << deg.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Degrader*>(deg.clone()));

}

/**
 *
 *
 * @param drift
 */
void ParallelCyclotronTracker::visitDrift(const Drift& drift) {
    *gmsg << "In drift L= " << drift.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Drift*>(drift.clone()));
}

/**
 *
 *
 *  @param
 */
void ParallelCyclotronTracker::visitFlexibleCollimator(const FlexibleCollimator&) {

}

/**
 *
 *
 * @param off
 */
void ParallelCyclotronTracker::visitOffset(const Offset& off) {
    if (opalRing_m == nullptr)
        throw OpalException(
                            "ParallelCylcotronTracker::visitOffset",
                            "Attempt to place an offset when Ring not defined");
    Offset* offNonConst = const_cast<Offset*>(&off);
    offNonConst->updateGeometry(opalRing_m->getNextPosition(),
                                opalRing_m->getNextNormal());
    opalRing_m->appendElement(off);
}

/**
 *
 *
 * @param marker
 */
void ParallelCyclotronTracker::visitMarker(const Marker& marker) {
    //   *gmsg << "In Marker; L= " << marker.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Marker*>(marker.clone()));
    // Do nothing.
}

/**
 *
 *
 * @param corr
 */
void ParallelCyclotronTracker::visitMonitor(const Monitor& corr) {
    //   *gmsg << "In Monitor; L= " << corr.getElementLength() << endl;
    myElements.push_back(dynamic_cast<Monitor*>(corr.clone()));
    //   applyDrift(flip_s * corr.getElementLength());
}

/**
 *
 *
 * @param mult
 */
void ParallelCyclotronTracker::visitMultipole(const Multipole& mult) {
    *gmsg << "In Multipole; L= " << mult.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<Multipole*>(mult.clone()));
}

/**
 *
 *
 * @param multT
 */
void ParallelCyclotronTracker::visitMultipoleT(const MultipoleT& multT) {
    *gmsg << "Adding MultipoleT" << endl;
    if (opalRing_m != nullptr) {
        opalRing_m->appendElement(multT);
    } else {
        throw OpalException("ParallelCyclotronTracker::visitMultipoleT",
                            "Need to define a RINGDEFINITION to use MultipoleT element");
    }
    myElements.push_back(dynamic_cast<MultipoleT*>(multT.clone()));
}

/**
 *
 *
 * @param multTstraight
 */
void ParallelCyclotronTracker::visitMultipoleTStraight(const MultipoleTStraight& multTstraight) {
    *gmsg << "Adding MultipoleTStraight" << endl;
    if (opalRing_m != nullptr) {
        opalRing_m->appendElement(multTstraight);
    } else {
        throw OpalException("ParallelCyclotronTracker::visitMultipoleTStraight",
                            "Need to define a RINGDEFINITION to use MultipoleTStraight element");
    }
    myElements.push_back(dynamic_cast<MultipoleTStraight*>(multTstraight.clone()));
}

/**
 *
 *
 * @param multTccurv
 */
void ParallelCyclotronTracker::visitMultipoleTCurvedConstRadius(const MultipoleTCurvedConstRadius& multTccurv) {
    *gmsg << "Adding MultipoleTCurvedConstRadius" << endl;
    if (opalRing_m != nullptr) {
        opalRing_m->appendElement(multTccurv);
    } else {
        throw OpalException("ParallelCyclotronTracker::visitMultipoleTCurvedConstRadius",
                            "Need to define a RINGDEFINITION to use MultipoleTCurvedConstRadius element");
    }
    myElements.push_back(dynamic_cast<MultipoleTCurvedConstRadius*>(multTccurv.clone()));
}

/**
 *
 *
 * @param multTvcurv
 */
void ParallelCyclotronTracker::visitMultipoleTCurvedVarRadius(const MultipoleTCurvedVarRadius& multTvcurv) {
    *gmsg << "Adding MultipoleTCurvedVarRadius" << endl;
    if (opalRing_m != nullptr) {
        opalRing_m->appendElement(multTvcurv);
    } else {
        throw OpalException("ParallelCyclotronTracker::visitMultipoleTCurvedVarRadius",
                            "Need to define a RINGDEFINITION to use MultipoleTCurvedVarRadius element");
    }
    myElements.push_back(dynamic_cast<MultipoleTCurvedVarRadius*>(multTvcurv.clone()));
}

/**
 *
 *
 * @param prob
 */
void ParallelCyclotronTracker::visitProbe(const Probe& prob) {
    *gmsg << "* ----------------------------- Probe ------------------------------------ *" << endl;

    Probe* elptr = dynamic_cast<Probe*>(prob.clone());
    myElements.push_back(elptr);

    *gmsg << "* Name    = " << elptr->getName() << endl;

    double xstart = elptr->getXStart();
    *gmsg << "* XStart  = " << xstart << " [m]" << endl;

    double xend = elptr->getXEnd();
    *gmsg << "* XEnd    = " << xend << " [m]" << endl;

    double ystart = elptr->getYStart();
    *gmsg << "* YStart  = " << ystart << " [m]" << endl;

    double yend = elptr->getYEnd();
    *gmsg << "* YEnd    = " << yend << " [m]" << endl;

    // initialise, do nothing
    elptr->initialise(itsBunch_m);

    double BcParameter[8] = {};

    BcParameter[0] = xstart;
    BcParameter[1] = xend;
    BcParameter[2] = ystart;
    BcParameter[3] = yend;
    BcParameter[4] = 1 ; // width

    buildupFieldList(BcParameter, ElementType::PROBE, elptr);
}

/**
 *
 *
 * @param bend
 */
void ParallelCyclotronTracker::visitRBend(const RBend& bend) {
    *gmsg << "In RBend; L= " << bend.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<RBend*>(bend.clone()));
}

/**
 *
 *
 * @param as
 */
void ParallelCyclotronTracker::visitRFCavity(const RFCavity& as) {
    *gmsg << "* ----------------------------- RFCavity --------------------------------- * " << endl;

    RFCavity* elptr = dynamic_cast<RFCavity*>(as.clone());
    myElements.push_back(elptr);

    if ( elptr->getCavityType() != CavityType::SGSW ) {
        throw OpalException("ParallelCyclotronTracker::visitRFCavity",
                            "\"" + elptr->getCavityTypeString() + "\" is not valid \"TYPE\" for RFCavity.\n"
                            "The ParallelCyclotronTracker can only play with cyclotron type RF system currently...");
    }

    *gmsg << "* Name                      = " << elptr->getName() << endl;

    std::string fmfn = elptr->getFieldMapFN();
    *gmsg << "* RF Field map file         = '" << fmfn << "'" << endl;

    double rmin = elptr->getRmin();
    *gmsg << "* Minimal radius of cavity  = " << rmin << " [mm]" << endl;

    double rmax = elptr->getRmax();
    *gmsg << "* Maximal radius of cavity  = " << rmax << " [mm]" << endl;

    double rff = elptr->getCycFrequency();
    *gmsg << "* RF frequency (2*pi*f)     = " << rff << " [rad/s]" << endl;

    double angle = elptr->getAzimuth();
    *gmsg << "* Cavity azimuth position   = " << angle << " [deg] " << endl;

    double gap = elptr->getGapWidth();
    *gmsg << "* Cavity gap width          = " << gap << " [mm] " << endl;

    double pdis = elptr->getPerpenDistance();
    *gmsg << "* Cavity Shift distance     = " << pdis << " [mm] " << endl;

    double phi0 = elptr->getPhi0();
    *gmsg << "* Initial RF phase (t=0)    = " << phi0 << " [deg] " << endl;

    /*
      Setup time dependence and in case of no
      timedependence use a polynom with  a_0 = 1 and a_k = 0, k = 1,2,3.
    */

    std::shared_ptr<AbstractTimeDependence> freq_atd = nullptr;
    std::shared_ptr<AbstractTimeDependence> ampl_atd = nullptr;
    std::shared_ptr<AbstractTimeDependence> phase_atd = nullptr;

    dvector_t  unityVec;
    unityVec.push_back(1.);
    unityVec.push_back(0.);
    unityVec.push_back(0.);
    unityVec.push_back(0.);

    std::string frequencyModelName = elptr->getFrequencyModelName();
    std::string amplitudeModelName = elptr->getAmplitudeModelName();
    std::string phaseModelName = elptr->getPhaseModelName();

    if (!frequencyModelName.empty()) {
        freq_atd = AbstractTimeDependence::getTimeDependence(frequencyModelName);
        *gmsg << "* Variable frequency RF Model name " << frequencyModelName << endl;
    } else {
        freq_atd = std::shared_ptr<AbstractTimeDependence>(new PolynomialTimeDependence(unityVec));
    }

    if (!amplitudeModelName.empty()) {
        ampl_atd = AbstractTimeDependence::getTimeDependence(amplitudeModelName);
        *gmsg << "* Variable amplitude RF Model name " << amplitudeModelName << endl;
    } else {
        ampl_atd = std::shared_ptr<AbstractTimeDependence>(new PolynomialTimeDependence(unityVec));
    }

    if (!phaseModelName.empty()) {
        phase_atd = AbstractTimeDependence::getTimeDependence(phaseModelName);
        *gmsg << "* Variable phase RF Model name " << phaseModelName << endl;
    } else {
        phase_atd = std::shared_ptr<AbstractTimeDependence>(new PolynomialTimeDependence(unityVec));
    }
    // read cavity voltage profile data from file.
    elptr->initialise(itsBunch_m, freq_atd, ampl_atd, phase_atd);

    double BcParameter[8] = {};

    BcParameter[0] = Units::mm2m * rmin;
    BcParameter[1] = Units::mm2m * rmax;
    BcParameter[2] = Units::mm2m * pdis;
    BcParameter[3] = angle;

    buildupFieldList(BcParameter, ElementType::RFCAVITY, elptr);
}

/**
 *
 * @param ring
 */
void ParallelCyclotronTracker::visitRing(const Ring& ring) {
    *gmsg << "* ----------------------------- Ring ------------------------------------- *" << endl;

    delete opalRing_m;

    opalRing_m = dynamic_cast<Ring*>(ring.clone());

    myElements.push_back(opalRing_m);

    opalRing_m->initialise(itsBunch_m);

    referenceR = opalRing_m->getBeamRInit();
    referencePr = opalRing_m->getBeamPRInit();
    referenceTheta = opalRing_m->getBeamPhiInit();

    if (referenceTheta <= -180.0 || referenceTheta > 180.0) {
        throw OpalException("Error in ParallelCyclotronTracker::visitRing",
                            "PHIINIT is out of [-180, 180)!");
    }

    referenceZ = 0.0;
    referencePz = 0.0;

    referencePtot = itsReference.getGamma() * itsReference.getBeta();
    referencePt = std::sqrt(referencePtot * referencePtot - referencePr * referencePr);

    if (referencePtot < 0.0)
        referencePt *= -1.0;

    sinRefTheta_m = std::sin(referenceTheta * Units::deg2rad);
    cosRefTheta_m = std::cos(referenceTheta * Units::deg2rad);

    double BcParameter[8] = {}; // zero initialise array

    buildupFieldList(BcParameter, ElementType::RING, opalRing_m);

    // Finally print some diagnostic
    *gmsg << "* Initial beam radius = " << referenceR << " [mm] " << endl;
    *gmsg << "* Initial gamma = " << itsReference.getGamma() << endl;
    *gmsg << "* Initial beta  = " << itsReference.getBeta() << endl;
    *gmsg << "* Total reference momentum      = " << referencePtot << " [beta gamma]" << endl;
    *gmsg << "* Reference azimuthal momentum  = " << referencePt << " [beta gamma]" << endl;
    *gmsg << "* Reference radial momentum     = " << referencePr << " [beta gamma]" << endl;
    *gmsg << "* " << opalRing_m->getSymmetry() << " fold field symmetry " << endl;
    *gmsg << "* Harmonic number h = " << opalRing_m->getHarmonicNumber() << " " << endl;
}

/**
 *
 *
 * @param bend
 */
void ParallelCyclotronTracker::visitSBend(const SBend& bend) {
    *gmsg << "In SBend; L = " << bend.getElementLength() << " however the element is missing " << endl;
    myElements.push_back(dynamic_cast<SBend*>(bend.clone()));
}

void ParallelCyclotronTracker::visitSBend3D(const SBend3D& bend) {
    *gmsg << "Adding SBend3D" << endl;
    if (opalRing_m != nullptr)
        opalRing_m->appendElement(bend);
    else
        throw OpalException("ParallelCyclotronTracker::visitSBend3D",
                            "Need to define a RINGDEFINITION to use SBend3D element");
}

void ParallelCyclotronTracker::visitScalingFFAMagnet(const ScalingFFAMagnet& bend) {
    *gmsg << "Adding ScalingFFAMagnet" << endl;
    if (opalRing_m != nullptr) {
        ScalingFFAMagnet* newBend = bend.clone(); // setup the end field, if required
        newBend->setupEndField();
        opalRing_m->appendElement(*newBend);
    } else {
        throw OpalException("ParallelCyclotronTracker::visitScalingFFAMagnet",
                            "Need to define a RINGDEFINITION to use ScalingFFAMagnet element");
    }
}

/**
 *
 *
 * @param sept
 */
void ParallelCyclotronTracker::visitSeptum(const Septum& sept) {
    *gmsg << "* ----------------------------- Septum ----------------------------------- *" << endl;

    Septum* elptr = dynamic_cast<Septum*>(sept.clone());
    myElements.push_back(elptr);

    *gmsg << "* Name    = " << elptr->getName() << endl;

    double xstart = elptr->getXStart();
    *gmsg << "* XStart  = " << xstart << " [m]" << endl;

    double xend = elptr->getXEnd();
    *gmsg << "* XEnd    = " << xend << " [m]" << endl;

    double ystart = elptr->getYStart();
    *gmsg << "* YStart  = " << ystart << " [m]" << endl;

    double yend = elptr->getYEnd();
    *gmsg << "* YEnd    = " << yend << " [m]" << endl;

    double width = elptr->getWidth();
    *gmsg << "* Width   = " << width << " [m]" << endl;

    // initialise, do nothing
    elptr->initialise(itsBunch_m);

    double BcParameter[8] = {};

    BcParameter[0] = xstart;
    BcParameter[1] = xend;
    BcParameter[2] = ystart;
    BcParameter[3] = yend;
    BcParameter[4] = width;

    buildupFieldList(BcParameter, ElementType::SEPTUM, elptr);
}

/**
 *
 *
 * @param solenoid
 */
void ParallelCyclotronTracker::visitSolenoid(const Solenoid& solenoid) {
    myElements.push_back(dynamic_cast<Solenoid*>(solenoid.clone()));
    Component* elptr = *(--myElements.end());
    if (!elptr->hasAttribute("ELEMEDGE")) {
        *gmsg << "Solenoid: no position of the element given!" << endl;
        return;
    }
}

/**
 *
 *
 * @param stripper
 */
void ParallelCyclotronTracker::visitStripper(const Stripper& stripper) {
    *gmsg << "* ----------------------------- Stripper --------------------------------- *" << endl;

    Stripper* elptr = dynamic_cast<Stripper*>(stripper.clone());
    myElements.push_back(elptr);

    *gmsg << "* Name    = " << elptr->getName() << endl;

    double xstart = elptr->getXStart();
    *gmsg << "* XStart  = " << xstart << " [m]" << endl;

    double xend = elptr->getXEnd();
    *gmsg << "* XEnd    = " << xend << " [m]" << endl;

    double ystart = elptr->getYStart();
    *gmsg << "* YStart  = " << ystart << " [m]" << endl;

    double yend = elptr->getYEnd();
    *gmsg << "* YEnd    = " << yend << " [m]" << endl;

    double opcharge = elptr->getOPCharge();
    *gmsg << "* Charge of outcoming particle = +e * " << opcharge << endl;

    double opmass = elptr->getOPMass();
    *gmsg << "* Mass of the outcoming particle = " << opmass << " [GeV/c^2]" << endl;

    bool stop = elptr->getStop();
    *gmsg << std::boolalpha
          << "* Particles stripped will be deleted after interaction -> "
          << stop << endl;

    elptr->initialise(itsBunch_m);

    double BcParameter[8] = {};

    BcParameter[0] = xstart;
    BcParameter[1] = xend;
    BcParameter[2] = ystart;
    BcParameter[3] = yend;
    BcParameter[4] = 1; // width
    BcParameter[5] = opcharge;
    BcParameter[6] = opmass;

    buildupFieldList(BcParameter, ElementType::STRIPPER, elptr);
}

/**
 *
 *
 * @param vac
 */
void ParallelCyclotronTracker::visitVacuum(const Vacuum& vac) {
    *gmsg << "* ----------------------------- Vacuum ----------------------------------- *" << endl;

    Vacuum* elptr = dynamic_cast<Vacuum*>(vac.clone());
    myElements.push_back(elptr);

    double BcParameter[8] = {};

    std::string gas = elptr->getResidualGasName();
    *gmsg << "* Residual gas = " << gas << endl;

    double pressure = elptr->getPressure();
    if ( boost::filesystem::exists(elptr->getPressureMapFN()) ) {
        std::string pmfn = elptr->getPressureMapFN();
        *gmsg << "* Pressure field map file = '" << pmfn << "'" << endl;
        *gmsg << "* Default pressure = " << pressure << " [mbar]" << endl;
    } else  {
        *gmsg << "* Pressure     = " << pressure << " [mbar]" << endl;
    }
    double pscale = elptr->getPScale();

    double temperature = elptr->getTemperature();
    *gmsg << "* Temperature  = " << temperature << " [K]" << endl;

    bool stop = elptr->getStop();
    *gmsg << std::boolalpha
          << "* Particles stripped will be deleted after interaction -> "
          << stop << endl;

    elptr->initialise(itsBunch_m);

    BcParameter[0] = pressure;
    BcParameter[1] = pscale;
    BcParameter[2] = temperature;

    buildupFieldList(BcParameter, ElementType::VACUUM, elptr);
}

/**
 *
 *
 * @param cav
 */
void ParallelCyclotronTracker::visitVariableRFCavity(const VariableRFCavity& cav) {
    *gmsg << "Adding Variable RF Cavity" << endl;
    if (opalRing_m != nullptr)
        opalRing_m->appendElement(cav);
    else
        throw OpalException("ParallelCyclotronTracker::visitVariableRFCavity",
                            "Need to define a RINGDEFINITION to use VariableRFCavity element");
}

/**
 *
 *
 * @param cav
 */
void ParallelCyclotronTracker::visitVariableRFCavityFringeField
                                  (const VariableRFCavityFringeField& cav) {
    *gmsg << "Adding Variable RF Cavity with Fringe Field" << endl;
    if (opalRing_m != nullptr)
        opalRing_m->appendElement(cav);
    else
        throw OpalException("ParallelCyclotronTracker::visitVariableRFCavityFringeField",
                            "Need to define a RINGDEFINITION to use VariableRFCavity element");
}

/**
 *
 *
 * @param mag
 */
void ParallelCyclotronTracker::visitVerticalFFAMagnet(const VerticalFFAMagnet& mag) {
    *gmsg << "Adding Vertical FFA Magnet" << endl;
    if (opalRing_m != nullptr)
        opalRing_m->appendElement(mag);
    else
        throw OpalException("ParallelCyclotronTracker::visitVerticalFFAMagnet",
                            "Need to define a RINGDEFINITION to use VerticalFFAMagnet element");
}


/**
 *
 *
 * @param BcParameter
 * @param ElementType
 * @param elptr
 */
void ParallelCyclotronTracker::buildupFieldList(double BcParameter[], ElementType elementType, Component *elptr) {
    beamline_list::iterator sindex;

    type_pair *localpair = new type_pair();
    localpair->first = elementType;

    for (int i = 0; i < 8; i++)
        *(((localpair->second).first) + i) = *(BcParameter + i);

    (localpair->second).second = elptr;

    // always put cyclotron as the first element in the list.
    if (elementType == ElementType::RING || elementType == ElementType::CYCLOTRON) {
        sindex = FieldDimensions.begin();
    } else {
        sindex = FieldDimensions.end();
    }
    FieldDimensions.insert(sindex, localpair);
}

/**
 *
 *
 * @param bl
 */
void ParallelCyclotronTracker::visitBeamline(const Beamline& bl) {
    const FlaggedBeamline* fbl = static_cast<const FlaggedBeamline*>(&bl);
    fbl->iterate(*this, false);//*dynamic_cast<BeamlineVisitor *>(this), false);
}

void ParallelCyclotronTracker::checkNumPart(std::string s) {
    int nlp = itsBunch_m->getLocalNum();
    int minnlp = 0;
    int maxnlp = 111111;
    reduce(nlp, minnlp, OpMinAssign());
    reduce(nlp, maxnlp, OpMaxAssign());
    *gmsg << s << "min local particle number: "<< minnlp << endl;
    *gmsg << "*                     max local particle number: " << maxnlp << endl;
}


void ParallelCyclotronTracker::execute() {
    /*
      Initialize common variables and structures
      for the integrators
    */
    if (OpalData::getInstance()->inRestartRun()) {
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::APPEND);
    }

    step_m          = 0;
    restartStep0_m  = 0;
    turnnumber_m    = 1;
    azimuth_m       = -1.0;
    prevAzimuth_m   = -1.0;

    // Record how many bunches have already been injected. ONLY FOR MPM
    if (isMultiBunch())
        mbHandler_m->setNumBunch(itsBunch_m->getNumBunch());

    itsBeamline->accept(*this);
    if (opalRing_m != nullptr)
        opalRing_m->lockRing();

    // Display the selected elements
    *gmsg << "* ------------------------------------------------------------------------ *" << endl;
    *gmsg << "* The selected Beam line elements are :" << endl;
    for (auto fd : FieldDimensions) {
        ElementType type = fd->first;
        *gmsg << "* -> " <<  ElementBase::getTypeString(type) << endl;
        if (type == ElementType::RFCAVITY) {
            RFCavity* cav = static_cast<RFCavity*>((fd->second).second);
            CavityCrossData ccd = {cav, cav->getSinAzimuth(), cav->getCosAzimuth(), cav->getPerpenDistance() * Units::mm2m};
            cavCrossDatas_m.push_back(ccd);
        } else if ( type == ElementType::CCOLLIMATOR ||
                    type == ElementType::PROBE       ||
                    type == ElementType::SEPTUM      ||
                    type == ElementType::STRIPPER) {
            PluginElement* element = static_cast<PluginElement*>((fd->second).second);
            pluginElements_m.push_back(element);
        }
    }
    *gmsg << "* ------------------------------------------------------------------------ *" << endl << endl;

    // Get BoundaryGeometry that is already initialized
    bgf_m = OpalData::getInstance()->getGlobalGeometry();
    if (bgf_m)
        lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(bgf_m->getOpalName(),!Options::asciidump));

    // External field arrays for dumping
    for (int k = 0; k < 2; k++)
        FDext_m[k] = Vector_t(0.0, 0.0, 0.0);

    extE_m = Vector_t(0.0, 0.0, 0.0);
    extB_m = Vector_t(0.0, 0.0, 0.0);
    DumpFields::writeFields((((*FieldDimensions.begin())->second).second));
    DumpEMFields::writeFields((((*FieldDimensions.begin())->second).second));

    function_t func = std::bind(&ParallelCyclotronTracker::getFieldsAtPoint,
                                this,
                                std::placeholders::_1,
                                std::placeholders::_2,
                                std::placeholders::_3,
                                std::placeholders::_4);

    switch ( stepper_m ) {
        case Steppers::TimeIntegrator::LF2: {
            *gmsg << "* 2nd order Leap-Frog integrator" << endl;
            itsStepper_mp.reset(new LF2<function_t>(func));
            break;
        }
        case Steppers::TimeIntegrator::MTS: {
            *gmsg << "* Multiple time stepping (MTS) integrator" << endl;
            break;
        }
        case Steppers::TimeIntegrator::RK4:
        default: {
            *gmsg << "* 4th order Runge-Kutta integrator" << endl;
            itsStepper_mp.reset(new RK4<function_t>(func));
            break;
        }
    }

    if ( stepper_m == Steppers::TimeIntegrator::MTS) {
        MtsTracker();
    } else {
        GenericTracker();
    }

    *gmsg << "* ------------------------------------------------------------------------ *" << endl;
    *gmsg << "* Finalizing i.e. write data and close files :" << endl;
    for (auto fd : FieldDimensions) {
        ((fd->second).second)->finalise();
    }
    *gmsg << "* ------------------------------------------------------------------------ *" << endl;
}


void ParallelCyclotronTracker::MtsTracker() {
    /*
     * variable             unit        meaning
     * ------------------------------------------------
     * t                    [ns]        time
     * dt                   [ns]        time step
     * oldReferenceTheta    [rad]       azimuth angle
     * itsBunch_m->R        [m]         particle position
     *
     */

    double t = 0, dt = 0, oldReferenceTheta = 0;
    std::tie(t, dt, oldReferenceTheta) = initializeTracking_m();

    int const numSubsteps = std::max(Options::mtsSubsteps, 1);
    *gmsg << "* MTS: Number of substeps per step is " << numSubsteps << endl;

    double const dt_inner = dt / double(numSubsteps);
    *gmsg << "* MTS: The inner time step is therefore " << dt_inner << " [ns]" << endl;

//     int SteptoLastInj = itsBunch_m->getSteptoLastInj();

    bool flagTransition = false; // flag to determine when to transit from single-bunch to multi-bunches mode

    *gmsg << "* ---------------------- Start tracking ---------------------------------- *" << endl;

    if ( itsBunch_m->hasFieldSolver() )
        computeSpaceChargeFields_m();

    for (; (step_m < maxSteps_m) && (itsBunch_m->getTotalNum()>0); step_m++) {

        bool finishedTurn = false;

        if (step_m % Options::sptDumpFreq == 0) {
            singleParticleDump();
        }

        Ippl::Comm->barrier();

        // First half kick from space charge force
        if (itsBunch_m->hasFieldSolver()) {
            kick(0.5 * dt);
        }

        // Substeps for external field integration
        for (int n = 0; n < numSubsteps; ++n)
            borisExternalFields(dt_inner);

        // bunch injection
        // TODO: Where is correct location for this piece of code? Beginning/end of step? Before field solve?
        injectBunch(flagTransition);

        if ( itsBunch_m->hasFieldSolver() ) {
            computeSpaceChargeFields_m();
        } else {
            // if field solver is not available , only update bunch, to transfer particles between nodes if needed,
            // reset parameters such as LocalNum, initialTotalNum_m.
            // INFOMSG("No space charge Effects are included!"<<endl;);
            if ((step_m % Options::repartFreq * 100) == 0) { //TODO: why * 100?
                Vector_t const meanP = calcMeanP();
                double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * Physics::pi;
                Vector_t const meanR = calcMeanR();
                globalToLocal(itsBunch_m->R, phi, meanR);
                itsBunch_m->updateNumTotal();
                repartition();
                localToGlobal(itsBunch_m->R, phi, meanR);
            }
        }

        // Second half kick from space charge force
        if (itsBunch_m->hasFieldSolver())
            kick(0.5 * dt);

        // recalculate bingamma and reset the BinID for each particles according to its current gamma
        if (isMultiBunch() && (step_m % Options::rebinFreq == 0)) {
            mbHandler_m->updateParticleBins(itsBunch_m);
        }

        // dump some data after one push in single particle tracking
        if ( mode_m == TrackingMode::SINGLE ) {

            unsigned int i = 0;

            double temp_meanTheta = calculateAngle2(itsBunch_m->R[i](0),
                                                    itsBunch_m->R[i](1)); // [-pi, pi]

            dumpThetaEachTurn_m(t, itsBunch_m->R[i], itsBunch_m->P[i],
                                temp_meanTheta, finishedTurn);

            dumpAzimuthAngles_m(t, itsBunch_m->R[i], itsBunch_m->P[i],
                                oldReferenceTheta, temp_meanTheta);

            oldReferenceTheta = temp_meanTheta;
        } else if ( mode_m == TrackingMode::BUNCH ) {
            // both for single bunch and multi-bunch
            // avoid dump at the first step
            // finishedTurn has not been changed in first push
            if ( isTurnDone() ) {
                ++turnnumber_m;
                finishedTurn = true;

                *gmsg << endl;
                *gmsg << "*** Finished turn " << turnnumber_m - 1
                      << ", Total number of live particles: "
                      << itsBunch_m->getTotalNum() << endl;
            }

            // Recalculate bingamma and reset the BinID for each particles according to its current gamma
            if (isMultiBunch() && (step_m % Options::rebinFreq == 0)) {
                mbHandler_m->updateParticleBins(itsBunch_m);
            }
        }

        // printing + updating bunch parameters + updating t
        update_m(t, dt, finishedTurn);

    }

    // Some post-integration stuff
    *gmsg << endl;
    *gmsg << "* ---------------------- DONE TRACKING PARTICLES ------------------------- *" << endl;


    //FIXME
    dvector_t Ttime = dvector_t();
    dvector_t Tdeltr = dvector_t();
    dvector_t Tdeltz = dvector_t();
    ivector_t TturnNumber = ivector_t();

    finalizeTracking_m(Ttime, Tdeltr, Tdeltz, TturnNumber);
}


void ParallelCyclotronTracker::GenericTracker() {
    /*
     * variable             unit        meaning
     * ------------------------------------------------
     * t                    [ns]        time
     * dt                   [ns]        time step
     * oldReferenceTheta    [rad]       azimuth angle
     * itsBunch_m->R        [m]         particle position
     *
     */
    // Generic Tracker that has three modes defined by timeIntegrator_m:
    // 0 --- RK-4 (default)
    // 1 --- LF-2
    // (2 --- MTS ... not yet implemented in generic tracker)
    // mbHandler_m->getNumBunch() determines the number of bunches in multibunch mode (MBM, 1 for OFF)
    // Total number of particles determines single particle mode (SPM, 1 particle) or
    // Static Equilibrium Orbit Mode (SEO, 2 particles)

    double t = 0, dt = 0, oldReferenceTheta = 0;
    std::tie(t, dt, oldReferenceTheta) = initializeTracking_m();

    // If initialTotalNum_m = 2, trigger SEO mode and prepare for transverse tuning calculation
    // Where is the IF here? -DW
    dvector_t Ttime, Tdeltr, Tdeltz;
    ivector_t TturnNumber;

    // Apply the plugin elements: probe, collimator, stripper, septum once before first step
    bool flagNeedUpdate = applyPluginElements(dt);

    // Destroy particles if they are marked as Bin = -1 in the plugin elements
    // or out of global aperture
    deleteParticle(flagNeedUpdate);

    /********************************
     *     Main integration loop    *
     ********************************/
    *gmsg << endl;
    *gmsg << "* ---------------------- Start tracking ---------------------------------- *" << endl;

    for (; (step_m < maxSteps_m) && (itsBunch_m->getTotalNum()>0); step_m++) {

        bool finishedTurn = false;

        switch (mode_m) {
            case TrackingMode::SEO: {
                // initialTotalNum_m == 2
                seoMode_m(t, dt, finishedTurn, Ttime, Tdeltr, Tdeltz, TturnNumber);
                break;
            }
            case TrackingMode::SINGLE: {
                // initialTotalNum_m == 1
                singleMode_m(t, dt, finishedTurn, oldReferenceTheta);
                break;
            }
            case TrackingMode::BUNCH: {
                // initialTotalNum_m > 2
                // Start Tracking for number of particles > 2 (i.e. not single and not SEO mode)
                bunchMode_m(t, dt, finishedTurn);
                break;
            }
            case TrackingMode::UNDEFINED: {
            default:
                throw OpalException("ParallelCyclotronTracker::GenericTracker()",
                                    "No such tracking mode.");
            }
        }
        // Update bunch and some parameters and output some info
        update_m(t, dt, finishedTurn);

    } // end for: the integration is DONE after maxSteps_m steps or if all particles are lost!

    // Some post-integration stuff
    *gmsg << endl;
    *gmsg << "* ---------------------- DONE TRACKING PARTICLES ------------------------- *" << endl;

    finalizeTracking_m(Ttime, Tdeltr, Tdeltz, TturnNumber);
}

bool ParallelCyclotronTracker::getFieldsAtPoint(const double& t, const size_t& Pindex, Vector_t& Efield, Vector_t& Bfield) {

    bool outOfBound = this->computeExternalFields_m(Pindex, t, Efield, Bfield);

    // For runs without space charge effects, override this step to save time
    if (itsBunch_m->hasFieldSolver()) {

        // Don't do for reference particle
        if (itsBunch_m->ID[Pindex] != 0) {

            // add external Field and self space charge field
            Efield += itsBunch_m->Ef[Pindex];
            Bfield += itsBunch_m->Bf[Pindex];
        }
    }

    return outOfBound;
}


/**
 *
 *
 * @param Rold
 * @param Rnew
 * @param elptr
 * @param Dold
 *
 * @return
 */

bool ParallelCyclotronTracker::checkGapCross(Vector_t Rold, Vector_t Rnew,
                                             RFCavity * rfcavity, double& Dold)
{
    bool flag = false;
    double sinx = rfcavity->getSinAzimuth();
    double cosx = rfcavity->getCosAzimuth();
    // TODO: Presumably this is still in mm, so for now, change to m -DW
    double PerpenDistance = Units::mm2m * rfcavity->getPerpenDistance();
    double distNew = (Rnew[0] * sinx - Rnew[1] * cosx) - PerpenDistance;
    double distOld = (Rold[0] * sinx - Rold[1] * cosx) - PerpenDistance;
    if (distOld > 0.0 && distNew <= 0.0) flag = true;
    // This parameter is used correct cavity phase
    Dold = Units::m2mm * distOld;
    return flag;
}

bool ParallelCyclotronTracker::RFkick(RFCavity * rfcavity, const double t, const double dt, const int Pindex){
    // For OPAL 2.0: As long as the RFCavity is in mm, we have to change R to mm here -DW
    double radius = std::sqrt(std::pow(Units::m2mm * itsBunch_m->R[Pindex](0), 2.0) + std::pow(Units::m2mm * itsBunch_m->R[Pindex](1), 2.0)
                         - std::pow(rfcavity->getPerpenDistance() , 2.0));
    double rmin = rfcavity->getRmin();
    double rmax = rfcavity->getRmax();
    double nomalRadius = (radius - rmin) / (rmax - rmin);
    double tempP[3];
    if (nomalRadius <= 1.0 && nomalRadius >= 0.0) {

        for (int j = 0; j < 3; j++)
            tempP[j] = itsBunch_m->P[Pindex](j);  //[px,py,pz]  units: dimensionless

        // here evaluate voltage and conduct momenta kicking;
        rfcavity->getMomentaKick(nomalRadius, tempP, t, dt, itsBunch_m->ID[Pindex], itsBunch_m->getM(), itsBunch_m->getQ()); // t : ns

        for (int j = 0; j < 3; j++)
            itsBunch_m->P[Pindex](j) = tempP[j];
        return true;
    }
    return false;
}


struct adder {
    adder() : sum(0) {}
    double sum;
    void operator()(double x) { sum += x; }
};

/**
 *
 *
 * @param t
 * @param r
 * @param z
 * @param lastTurn
 * @param nur
 * @param nuz
 *
 * @return
 */
bool ParallelCyclotronTracker::getTunes(dvector_t& t, dvector_t& r, dvector_t& z,
                                        int lastTurn, double& /*nur*/, double& /*nuz*/) {
    TUNE_class *tune;

    int Ndat = t.size();

    /*
      remove mean
    */
    double rsum =  for_each(r.begin(), r.end(), adder()).sum;

    for (int i = 0; i < Ndat; i++)
        r[i] -= rsum;

    double zsum =  for_each(z.begin(), z.end(), adder()).sum;

    for (int i = 0; i < Ndat; i++)
        z[i] -= zsum;
    double ti = *(t.begin());
    double tf = t[t.size()-1];
    double T = (tf - ti);

    t.clear();
    for (int i = 0; i < Ndat; i++) {
        t.push_back(i);
    }

    T = t[Ndat-1];

    *gmsg << endl;
    *gmsg << "* ************************************* nuR ******************************************* *" << endl;
    *gmsg << endl << "* ===> " << Ndat << " data points  Ti=" << ti << " Tf= " << tf << " -> T= " << T << endl;

    int nhis_lomb = 10;
    int  stat = 0;
    // book tune class
    tune = new TUNE_class();
    stat = tune->lombAnalysis(t, r, nhis_lomb, T / lastTurn);
    if (stat != 0)
        *gmsg << "* TUNE: Lomb analysis failed" << endl;
    *gmsg << "* ************************************************************************************* *" << endl;

    delete tune;
    tune = nullptr;
    // FIXME: FixMe: need to come from the inputfile
    nhis_lomb = 100;

    if (zsum != 0.0) {

        *gmsg << endl;
        *gmsg << "* ************************************* nuZ ******************************************* *" << endl;
        *gmsg << endl << "* ===> " << Ndat << " data points  Ti=" << ti << " Tf= " << tf << " -> T= " << T << endl;

        // book tune class
        tune = new TUNE_class();
        stat = tune->lombAnalysis(t, z, nhis_lomb, T / lastTurn);
        if (stat != 0)
            *gmsg << "* TUNE: Lomb analysis failed" << endl;
        *gmsg << "* ************************************************************************************* *" << endl;

        delete tune;
        tune = nullptr;
    }
    return true;
}


double ParallelCyclotronTracker::getHarmonicNumber() const {
    if (opalRing_m != nullptr)
        return opalRing_m->getHarmonicNumber();
    Cyclotron* elcycl = dynamic_cast<Cyclotron*>(((*FieldDimensions.begin())->second).second);
    if (elcycl != nullptr)
        return elcycl->getCyclHarm();
    throw OpalException("ParallelCyclotronTracker::getHarmonicNumber()",
                        std::string("The first item in the FieldDimensions list does not ")
                        +std::string("seem to be a Ring or a Cyclotron element"));
}


Vector_t ParallelCyclotronTracker::calcMeanR(short bunchNr) const {
    Vector_t meanR(0.0, 0.0, 0.0);

    for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {
        // take all particles if bunchNr <= -1
        if ( bunchNr > -1 && itsBunch_m->bunchNum[i] != bunchNr)
            continue;

        for (int d = 0; d < 3; ++d) {
            meanR(d) += itsBunch_m->R[i](d);
        }
    }

    reduce(meanR, meanR, OpAddAssign());

    size_t n = itsBunch_m->getTotalNum();

    if ( bunchNr > -1 )
        n = itsBunch_m->getTotalNumPerBunch(bunchNr);

    return meanR / Vector_t(n);
}

Vector_t ParallelCyclotronTracker::calcMeanP() const {
    Vector_t meanP(0.0, 0.0, 0.0);

    for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {
        for (int d = 0; d < 3; ++d) {
            meanP(d) += itsBunch_m->P[i](d);
        }
    }

    reduce(meanP, meanP, OpAddAssign());
    return meanP / Vector_t(itsBunch_m->getTotalNum());
}

void ParallelCyclotronTracker::repartition() {
    if ((step_m % Options::repartFreq) == 0) {
        IpplTimings::startTimer(BinRepartTimer_m);
        itsBunch_m->do_binaryRepart();
        Ippl::Comm->barrier();
        IpplTimings::stopTimer(BinRepartTimer_m);
    }
}

void ParallelCyclotronTracker::globalToLocal(ParticleAttrib<Vector_t>& particleVectors,
                                             double phi, Vector_t const translationToGlobal) {
    IpplTimings::startTimer(TransformTimer_m);
    particleVectors -= translationToGlobal;

    Tenzor<double, 3> const rotation( std::cos(phi), std::sin(phi), 0,
                                      -std::sin(phi), std::cos(phi), 0,
                                      0,        0, 1); // clockwise rotation

    for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {

        particleVectors[i] = dot(rotation, particleVectors[i]);
    }
    IpplTimings::stopTimer(TransformTimer_m);
}

void ParallelCyclotronTracker::localToGlobal(ParticleAttrib<Vector_t>& particleVectors,
                                             double phi, Vector_t const translationToGlobal) {
    IpplTimings::startTimer(TransformTimer_m);
    Tenzor<double, 3> const rotation(std::cos(phi), -std::sin(phi), 0,
                                     std::sin(phi),  std::cos(phi), 0,
                                     0,         0, 1); // counter-clockwise rotation

    for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {

        particleVectors[i] = dot(rotation, particleVectors[i]);
    }

    particleVectors += translationToGlobal;
    IpplTimings::stopTimer(TransformTimer_m);
}


inline void ParallelCyclotronTracker::globalToLocal(ParticleAttrib<Vector_t>& particleVectors,
                                                    Quaternion_t const quaternion,
                                                    Vector_t const meanR) {
    IpplTimings::startTimer(TransformTimer_m);

    // Translation from global to local
    particleVectors -= meanR;

    // Rotation to align P_mean with x-axis
    rotateWithQuaternion(particleVectors, quaternion);
    IpplTimings::stopTimer(TransformTimer_m);
}

inline void ParallelCyclotronTracker::localToGlobal(ParticleAttrib<Vector_t>& particleVectors,
                                                    Quaternion_t const quaternion,
                                                    Vector_t const meanR) {
    IpplTimings::startTimer(TransformTimer_m);
    // Reverse the quaternion by multiplying the axis components (x,y,z) with -1
    Quaternion_t reverseQuaternion = quaternion * -1.0;
    reverseQuaternion(0) *= -1.0;

    // Rotation back to original P_mean direction
    rotateWithQuaternion(particleVectors, reverseQuaternion);

    // Translation from local to global
    particleVectors += meanR;
    IpplTimings::stopTimer(TransformTimer_m);
}

inline void ParallelCyclotronTracker::globalToLocal(ParticleAttrib<Vector_t>& particleVectors,
                                                    double const phi,
                                                    double const psi,
                                                    Vector_t const meanR) {

    IpplTimings::startTimer(TransformTimer_m);
    //double const tolerance = 1.0e-4; // TODO What is a good angle threshold? -DW

    // Translation from global to local
    particleVectors -= meanR;

    // Rotation to align P_mean with x-axis
    rotateAroundZ(particleVectors, phi);

    // If theta is large enough (i.e. there is significant momentum in z direction)
    // rotate around x-axis next
    //if (fabs(psi) > tolerance)
    rotateAroundX(particleVectors, psi);
    IpplTimings::stopTimer(TransformTimer_m);
}

inline void ParallelCyclotronTracker::globalToLocal(Vector_t& myVector,
                                                    double const phi,
                                                    double const psi,
                                                    Vector_t const meanR) {

    IpplTimings::startTimer(TransformTimer_m);
    //double const tolerance = 1.0e-4; // TODO What is a good angle threshold? -DW

    // Translation from global to local
    myVector -= meanR;

    // Rotation to align P_mean with x-axis
    rotateAroundZ(myVector, phi);

    // If theta is large enough (i.e. there is significant momentum in z direction)
    // rotate around x-axis next
    //if (fabs(psi) > tolerance)
    rotateAroundX(myVector, psi);
    IpplTimings::stopTimer(TransformTimer_m);
}

inline void ParallelCyclotronTracker::localToGlobal(ParticleAttrib<Vector_t>& particleVectors,
                                                    double const phi,
                                                    double const psi,
                                                    Vector_t const meanR) {

    IpplTimings::startTimer(TransformTimer_m);
    //double const tolerance = 1.0e-4; // TODO What is a good angle threshold? -DW

    // If theta is large enough (i.e. there is significant momentum in z direction)
    // rotate back around x-axis next
    //if (fabs(psi) > tolerance)
    rotateAroundX(particleVectors, -psi);

    // Rotation to align P_mean with x-axis
    rotateAroundZ(particleVectors, -phi);

    // Translation from local to global
    particleVectors += meanR;
    IpplTimings::stopTimer(TransformTimer_m);
}

inline void ParallelCyclotronTracker::localToGlobal(Vector_t& myVector,
                                                    double const phi,
                                                    double const psi,
                                                    Vector_t const meanR) {

    IpplTimings::startTimer(TransformTimer_m);
    //double const tolerance = 1.0e-4; // TODO What is a good angle threshold? -DW

    // If theta is large enough (i.e. there is significant momentum in z direction)
    // rotate back around x-axis next
    //if (fabs(psi) > tolerance)
    rotateAroundX(myVector, -psi);

    // Rotation to align P_mean with x-axis
    rotateAroundZ(myVector, -phi);

    // Translation from local to global
    myVector += meanR;
    IpplTimings::stopTimer(TransformTimer_m);
}

inline void ParallelCyclotronTracker::rotateWithQuaternion(ParticleAttrib<Vector_t>& particleVectors, Quaternion_t const quaternion) {

    Vector_t const quaternionVectorComponent = Vector_t(quaternion(1), quaternion(2), quaternion(3));
    double const quaternionScalarComponent = quaternion(0);

    for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {

        particleVectors[i] = 2.0f * dot(quaternionVectorComponent, particleVectors[i]) * quaternionVectorComponent +
            (quaternionScalarComponent * quaternionScalarComponent -
             dot(quaternionVectorComponent, quaternionVectorComponent)) * particleVectors[i] + 2.0f *
            quaternionScalarComponent * cross(quaternionVectorComponent, particleVectors[i]);
    }
}

inline void ParallelCyclotronTracker::normalizeQuaternion(Quaternion_t& quaternion){

    double tolerance = 1.0e-10;
    double length2 = dot(quaternion, quaternion);

    if (std::abs(length2) > tolerance && std::abs(length2 - 1.0f) > tolerance) {

        double length = std::sqrt(length2);
        quaternion /= length;
    }
}

inline void ParallelCyclotronTracker::normalizeVector(Vector_t& vector) {

    double tolerance = 1.0e-10;
    double length2 = dot(vector, vector);

    if (std::abs(length2) > tolerance && std::abs(length2 - 1.0f) > tolerance) {

        double length = std::sqrt(length2);
        vector /= length;
    }
}

inline void ParallelCyclotronTracker::rotateAroundZ(ParticleAttrib<Vector_t>& particleVectors, double const phi) {
    // Clockwise rotation of particles 'particleVectors' by 'phi' around Z axis

    Tenzor<double, 3> const rotation( std::cos(phi), std::sin(phi), 0,
                                      -std::sin(phi), std::cos(phi), 0,
                                      0,        0, 1);

    for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {
        particleVectors[i] = dot(rotation, particleVectors[i]);
    }
}

inline void ParallelCyclotronTracker::rotateAroundZ(Vector_t& myVector, double const phi) {
    // Clockwise rotation of single vector 'myVector' by 'phi' around Z axis

    Tenzor<double, 3> const rotation( std::cos(phi), std::sin(phi), 0,
                                      -std::sin(phi), std::cos(phi), 0,
                                      0,        0, 1);

    myVector = dot(rotation, myVector);
}

inline void ParallelCyclotronTracker::rotateAroundX(ParticleAttrib<Vector_t>& particleVectors, double const psi) {
    // Clockwise rotation of particles 'particleVectors' by 'psi' around X axis

    Tenzor<double, 3> const rotation(1,  0,          0,
                                     0,  std::cos(psi), std::sin(psi),
                                     0, -std::sin(psi), std::cos(psi));

    for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {

        particleVectors[i] = dot(rotation, particleVectors[i]);
    }
}

inline void ParallelCyclotronTracker::rotateAroundX(Vector_t& myVector, double const psi) {
    // Clockwise rotation of single vector 'myVector' by 'psi' around X axis

    Tenzor<double, 3> const rotation(1,  0,          0,
                                     0,  std::cos(psi), std::sin(psi),
                                     0, -std::sin(psi), std::cos(psi));

    myVector = dot(rotation, myVector);
}

inline void ParallelCyclotronTracker::getQuaternionTwoVectors(Vector_t u, Vector_t v, Quaternion_t& quaternion) {
    // four vector (w,x,y,z) of the quaternion of P_mean with the positive x-axis

    normalizeVector(u);
    normalizeVector(v);

    double k_cos_theta = dot(u, v);
    double k = std::sqrt(dot(u, u) * dot(v, v));
    double tolerance1 = 1.0e-5;
    double tolerance2 = 1.0e-8;
    Vector_t resultVectorComponent;

    if (std::abs(k_cos_theta / k + 1.0) < tolerance1) {
        // u and v are almost exactly antiparallel so we need to do
        // 180 degree rotation around any vector orthogonal to u

        resultVectorComponent = cross(u, xaxis);

        // If by chance u is parallel to xaxis, use zaxis instead
        if (dot(resultVectorComponent, resultVectorComponent) < tolerance2) {
            resultVectorComponent = cross(u, zaxis);
        }

        double halfAngle = 0.5 * Physics::pi;
        double sinHalfAngle = std::sin(halfAngle);

        resultVectorComponent *= sinHalfAngle;

        k = 0.0;
        k_cos_theta = std::cos(halfAngle);

    } else {
        resultVectorComponent = cross(u, v);
    }

    quaternion(0) = k_cos_theta + k;
    quaternion(1) = resultVectorComponent(0);
    quaternion(2) = resultVectorComponent(1);
    quaternion(3) = resultVectorComponent(2);

    normalizeQuaternion(quaternion);
}


bool ParallelCyclotronTracker::push(double h) {
    /* h   [ns]
     * dt1 [ns]
     * dt2 [ns]
     */
    IpplTimings::startTimer(IntegrationTimer_m);

    h *= Units::ns2s;

    bool flagNeedUpdate = false;
    for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {
        Vector_t const oldR = itsBunch_m->R[i];
        double const gamma = std::sqrt(1.0 + dot(itsBunch_m->P[i], itsBunch_m->P[i]));
        double const c_gamma = Physics::c / gamma;
        Vector_t const v = itsBunch_m->P[i] * c_gamma;
        itsBunch_m->R[i] += h * v;
        for (const auto & ccd : cavCrossDatas_m) {
            double const distNew = (itsBunch_m->R[i][0] * ccd.sinAzimuth - itsBunch_m->R[i][1] * ccd.cosAzimuth) - ccd.perpenDistance;
            bool tagCrossing = false;
            double distOld;
            if (distNew <= 0.0) {
                distOld = (oldR[0] * ccd.sinAzimuth - oldR[1] * ccd.cosAzimuth) - ccd.perpenDistance;
                if (distOld > 0.0) tagCrossing = true;
            }
            if (tagCrossing) {
                double const dt1 = distOld / std::sqrt(dot(v, v));
                double const dt2 = h - dt1;

                // Retrack particle from the old postion to cavity gap point
                itsBunch_m->R[i] = oldR + dt1 * v;

                // Momentum kick
                //itsBunch_m->R[i] *= 1000.0; // RFkick uses [itsBunch_m->R] == mm
                RFkick(ccd.cavity, itsBunch_m->getT() * Units::s2ns, dt1, i);
                //itsBunch_m->R[i] *= 0.001;

                itsBunch_m->R[i] += dt2 * itsBunch_m->P[i] * c_gamma;
            }
        }
        flagNeedUpdate |= (itsBunch_m->Bin[i] < 0);
    }

    IpplTimings::stopTimer(IntegrationTimer_m);
    return flagNeedUpdate;
}


bool ParallelCyclotronTracker::kick(double h) {
    IpplTimings::startTimer(IntegrationTimer_m);

    bool flagNeedUpdate = false;
    BorisPusher pusher;
    double const q = itsBunch_m->Q[0] / Physics::q_e; // For now all particles have the same charge
    double const M = itsBunch_m->M[0] * Units::GeV2eV; // For now all particles have the same rest energy

    for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {

        pusher.kick(itsBunch_m->R[i], itsBunch_m->P[i],
                    itsBunch_m->Ef[i], itsBunch_m->Bf[i],
                    h * Units::ns2s, M, q);
        flagNeedUpdate |= (itsBunch_m->Bin[i] < 0);
    }
    IpplTimings::stopTimer(IntegrationTimer_m);
    return flagNeedUpdate;
}


void ParallelCyclotronTracker::borisExternalFields(double h) {
    // h in [ns]

    // push particles for first half step
    bool flagNeedUpdate = push(0.5 * h);

    // Evaluate external fields
    IpplTimings::startTimer(IntegrationTimer_m);
    for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); ++i) {

        itsBunch_m->Ef[i] = Vector_t(0.0, 0.0, 0.0);
        itsBunch_m->Bf[i] = Vector_t(0.0, 0.0, 0.0);

        computeExternalFields_m(i, itsBunch_m->getT() * Units::s2ns,
                                itsBunch_m->Ef[i], itsBunch_m->Bf[i]);
    }
    IpplTimings::stopTimer(IntegrationTimer_m);

    // Kick particles for full step
    flagNeedUpdate |= kick(h);

    // push particles for second half step
    flagNeedUpdate |= push(0.5 * h);

    // apply the plugin elements: probe, collimator, stripper, septum
    flagNeedUpdate |= applyPluginElements(h);
    // destroy particles if they are marked as Bin=-1 in the plugin elements or out of global apeture
    deleteParticle(flagNeedUpdate);
}


bool ParallelCyclotronTracker::applyPluginElements(const double dt) {
    IpplTimings::startTimer(PluginElemTimer_m);

    for (beamline_list::iterator sindex = ++(FieldDimensions.begin()); sindex != FieldDimensions.end(); ++sindex) {
        if (((*sindex)->first) == ElementType::VACUUM) {
            Vacuum* vac = static_cast<Vacuum*>(((*sindex)->second).second);
            vac->checkVacuum(itsBunch_m, cycl_m);
        }
    }

    bool flag = false;
    for (PluginElement* element : pluginElements_m) {
        bool tmp = element->check(itsBunch_m,
                                  turnnumber_m,
                                  itsBunch_m->getT(),
                                  dt);
        flag |= tmp;

        if ( tmp ) {
            itsBunch_m->updateNumTotal();
            *gmsg << "* Total number of particles after PluginElement= "
                  << itsBunch_m->getTotalNum() << endl;
        }
    }

    IpplTimings::stopTimer(PluginElemTimer_m);
    return flag;
}

bool ParallelCyclotronTracker::deleteParticle(bool flagNeedUpdate) {
    IpplTimings::startTimer(DelParticleTimer_m);
    // Update immediately if any particles are lost during this step

    if ((step_m + 1) % Options::delPartFreq != 0) {
        IpplTimings::stopTimer(DelParticleTimer_m);
        return false;
    }

    allreduce(flagNeedUpdate, 1, std::logical_or<bool>());

    if (flagNeedUpdate) {
        short bunchCount = itsBunch_m->getNumBunch();
        std::vector<size_t> locLostParticleNum(bunchCount, 0);

        const int leb = itsBunch_m->getLastemittedBin();
        std::unique_ptr<size_t[]> localBinCount;

        if ( isMultiBunch() ) {
            localBinCount = std::unique_ptr<size_t[]>(new size_t[leb]);
            for (int i = 0; i < leb; ++i)
                localBinCount[i] = 0;
        }

        for (unsigned int i = 0; i < itsBunch_m->getLocalNum(); i++) {
            if (itsBunch_m->Bin[i] < 0) {
                ++locLostParticleNum[itsBunch_m->bunchNum[i]];
                itsBunch_m->destroy(1, i);
            } else if ( isMultiBunch() ) {
                /* we need to count the local number of particles
                 * per energy bin
                 */
                ++localBinCount[itsBunch_m->Bin[i]];
            }
        }

        if ( isMultiBunch() ) {
            // set the local bin count
            for (int i = 0; i < leb; ++i) {
                itsBunch_m->setLocalBinCount(localBinCount[i], i);
            }
        }

        std::vector<size_t> localnum(bunchCount + 1);
        for (size_t i = 0; i < localnum.size() - 1; ++i) {
            localnum[i] = itsBunch_m->getLocalNumPerBunch(i) - locLostParticleNum[i];
            itsBunch_m->setLocalNumPerBunch(localnum[i], i);
        }

        /* We need to destroy the particles now
         * before we compute the means. We also
         * have to update the total number of particles
         * otherwise the statistics are wrong.
         */
        itsBunch_m->performDestroy(true);

        /* total number of particles of individual bunches
         * last index of vector contains total number over all
         * bunches, used as a check
         */
        std::vector<size_t> totalnum(bunchCount + 1);
        localnum[bunchCount] = itsBunch_m->getLocalNum();

        allreduce(localnum.data(), totalnum.data(), localnum.size(), std::plus<size_t>());
        itsBunch_m->setTotalNum(totalnum[bunchCount]);

        for (short i = 0; i < bunchCount; ++i) {
            itsBunch_m->setTotalNumPerBunch(totalnum[i], i);
        }

        size_t sum = std::accumulate(totalnum.begin(),
                                     totalnum.end() - 1, 0);

        if ( sum != totalnum[bunchCount] ) {
            throw OpalException("ParallelCyclotronTracker::deleteParticle()",
                                "Total number of particles " + std::to_string(totalnum[bunchCount]) +
                                " != " + std::to_string(sum) + " (sum over all bunches)");
        }

        size_t globLostParticleNum = 0;
        size_t locNumLost = std::accumulate(locLostParticleNum.begin(),
                                            locLostParticleNum.end(), 0);

        reduce(locNumLost, globLostParticleNum, 1, std::plus<size_t>());

        if ( globLostParticleNum > 0 ) {
            *gmsg << level3 << "At step " << step_m
                  << ", lost "  << globLostParticleNum << " particles" << endl;
        }

        if (totalnum[bunchCount] == 0) {
            IpplTimings::stopTimer(DelParticleTimer_m);
            return flagNeedUpdate;
        }

        Vector_t const meanR = calcMeanR();
        Vector_t const meanP = calcMeanP();

        // Bunch (local) azimuth at meanR w.r.t. y-axis
        double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * Physics::pi;

        // Bunch (local) elevation at meanR w.r.t. xy plane
        double const psi = 0.5 * Physics::pi - std::acos(meanP(2) / std::sqrt(dot(meanP, meanP)));

        // For statistics data, the bunch is transformed into a local coordinate system
        // with meanP in direction of y-axis -DW
        globalToLocal(itsBunch_m->R, phi, psi, meanR);
        globalToLocal(itsBunch_m->P, phi, psi, Vector_t(0.0)); // P should be rotated, but not shifted

        //itsBunch_m->R *= Vector_t(0.001); // mm --> m

        // Now destroy particles and update pertinent parameters in local frame
        // Note that update() will be called within boundp() -DW
        itsBunch_m->boundp();
        //itsBunch_m->update();

        itsBunch_m->calcBeamParameters();

        //itsBunch_m->R *= Vector_t(1000.0); // m --> mm

        localToGlobal(itsBunch_m->R, phi, psi, meanR);
        localToGlobal(itsBunch_m->P, phi, psi, Vector_t(0.0));

        // If particles were deleted, recalculate bingamma and reset BinID for remaining particles
        if ( isMultiBunch() )
            mbHandler_m->updateParticleBins(itsBunch_m);
    }

    IpplTimings::stopTimer(DelParticleTimer_m);
    return flagNeedUpdate;
}

void ParallelCyclotronTracker::initTrackOrbitFile() {

    std::string f = OpalData::getInstance()->getInputBasename() + std::string("-trackOrbit.dat");

    outfTrackOrbit_m.setf(std::ios::scientific, std::ios::floatfield);
    outfTrackOrbit_m.precision(8);

    if (myNode_m == 0) {

        if (OpalData::getInstance()->inRestartRun()) {
            outfTrackOrbit_m.open(f.c_str(), std::ios::app);
            outfTrackOrbit_m << "# Restart at integration step " << itsBunch_m->getLocalTrackStep() << std::endl;
        } else {
            outfTrackOrbit_m.open(f.c_str());
            outfTrackOrbit_m << "# The six-dimensional phase space data in the global Cartesian coordinates" << std::endl;
            outfTrackOrbit_m << "# Part. ID    x [m]       beta_x*gamma       y [m]      beta_y*gamma        z [m]      beta_z*gamma" << std::endl;
        }
    }
}

void ParallelCyclotronTracker::initDistInGlobalFrame() {

    if (!OpalData::getInstance()->inRestartRun()) {
        // Start a new run (no restart)

        double const initialReferenceTheta = referenceTheta * Units::deg2rad;

        // TODO: Replace with TracerParticle
        // Force the initial phase space values of the particle with ID = 0 to zero,
        // to set it as a reference particle.
        if (initialTotalNum_m > 2) {
            for (size_t i = 0; i < initialLocalNum_m; ++i) {
                if (itsBunch_m->ID[i] == 0) {
                    itsBunch_m->R[i] = Vector_t(0.0);
                    itsBunch_m->P[i] = Vector_t(0.0);
                }
            }
        }

        // Initialize global R
        //itsBunch_m->R *= Vector_t(1000.0); // m --> mm

        Vector_t const initMeanR = Vector_t(Units::mm2m * referenceR * cosRefTheta_m,
                                            Units::mm2m * referenceR * sinRefTheta_m,
                                            Units::mm2m * referenceZ);

        localToGlobal(itsBunch_m->R, initialReferenceTheta, initMeanR);

        // Initialize global P (Cartesian, but input P_ref is in Pr, Ptheta, Pz,
        // so translation has to be done before the rotation this once)
        // Cave: In the local frame, the positive y-axis is the direction of movement -DW
        for (size_t i = 0; i < initialLocalNum_m; ++i) {
            itsBunch_m->P[i](0) += referencePr;
            itsBunch_m->P[i](1) += referencePt;
            itsBunch_m->P[i](2) += referencePz;
        }

        // Out of the three coordinates of meanR (R, Theta, Z) only the angle
        // changes the momentum vector...
        localToGlobal(itsBunch_m->P, initialReferenceTheta);

        DistributionType distType = itsBunch_m->getDistType();
        if (distType == DistributionType::FROMFILE) {
            checkFileMomentum();
        }

        // Initialize the bin number of the first bunch to 0
        for (size_t i = 0; i < initialLocalNum_m; ++i) {
            itsBunch_m->Bin[i] = 0;
        }

        // Set time step per particle
        setTimeStep();

        // Backup initial distribution if multi bunch mode
        if ((initialTotalNum_m > 2) && isMultiBunch() && mbHandler_m->isForceMode()) {
            mbHandler_m->saveBunch(itsBunch_m);
        }

        // Else: Restart from the distribution in the h5 file
    } else {

        // Do a local frame restart (we have already checked that the old h5 file was saved in local
        // frame as well).
        if ((Options::psDumpFrame != DumpFrame::GLOBAL)) {

            *gmsg << "* Restart in the local frame" << endl;
            //itsBunch_m->R *= Vector_t(1000.0); // m --> mm

            Vector_t const initMeanR = Vector_t(Units::mm2m * referenceR * cosRefTheta_m,
                                                Units::mm2m * referenceR * sinRefTheta_m,
                                                Units::mm2m * referenceZ);

            // Do the transformations
            localToGlobal(itsBunch_m->R, referencePhi, referencePsi, initMeanR);
            localToGlobal(itsBunch_m->P, referencePhi, referencePsi);

            // Initialize the bin number of the first bunch to 0
            for (size_t i = 0; i < initialLocalNum_m; ++i) {
                itsBunch_m->Bin[i] = 0;
            }

            // Or do a global frame restart (no transformations necessary)
        } else {
            *gmsg << "* Restart in the global frame" << endl;

            pathLength_m = itsBunch_m->get_sPos();
            //itsBunch_m->R *= Vector_t(1000.0); // m --> mm
        }
    }

    // set the number of particles per bunch
    itsBunch_m->countTotalNumPerBunch();

    // ------------ Get some Values ---------------------------------------------------------- //
    Vector_t const meanR = calcMeanR();
    Vector_t const meanP = calcMeanP();

    // Bunch (local) azimuth at meanR w.r.t. y-axis
    double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * Physics::pi;

    // Bunch (local) elevation at meanR w.r.t. xy plane
    double const psi = 0.5 * Physics::pi - std::acos(meanP(2) / std::sqrt(dot(meanP, meanP)));

    double radius = std::sqrt(meanR[0] * meanR[0] + meanR[1] * meanR[1]);  // [m]

    if ( isMultiBunch() ) {
        mbHandler_m->setRadiusTurns(radius);
    }

    // Do boundp and repartition in the local frame at beginning of this run
    globalToLocal(itsBunch_m->R, phi, psi, meanR);
    globalToLocal(itsBunch_m->P, phi, psi); // P should be rotated, but not shifted

    //itsBunch_m->R *= Vector_t(0.001); // mm --> m

    itsBunch_m->boundp();

    checkNumPart(std::string("\n* Before repartition: "));
    repartition();
    checkNumPart(std::string("\n* After repartition:  "));

    itsBunch_m->calcBeamParameters();

    *gmsg << endl << "* *********************** Bunch information in local frame: ************************";
    *gmsg << *itsBunch_m << endl;

    //itsBunch_m->R *= Vector_t(1000.0); // m --> mm

    localToGlobal(itsBunch_m->R, phi, psi, meanR);
    localToGlobal(itsBunch_m->P, phi, psi);

    // Save initial distribution if not a restart
    if (!OpalData::getInstance()->inRestartRun()) {
        step_m -= 1;

        bunchDumpPhaseSpaceData();
        bunchDumpStatData();

        step_m += 1;
    }

    // Print out the Bunch information at beginning of the run. Because the bunch information
    // displays in units of m we have to change back and forth one more time.
    //itsBunch_m->R *= Vector_t(0.001); // mm --> m

    itsBunch_m->calcBeamParameters();

    // multi-bunch simulation only
    saveInjectValues();

    *gmsg << endl << "* *********************** Bunch information in global frame: ***********************";
    *gmsg << *itsBunch_m << endl;

    //itsBunch_m->R *= Vector_t(1000.0); // m --> mm
}

void ParallelCyclotronTracker::setTimeStep() {
    for (size_t i = 0; i < initialLocalNum_m; ++i) {
        itsBunch_m->dt[i] = itsBunch_m->getdT();
    }
}

void ParallelCyclotronTracker::checkFileMomentum() {

    double pTotalMean = 0.0;
    for (size_t i = 0; i < initialLocalNum_m; ++i) {
        pTotalMean += euclidean_norm(itsBunch_m->P[i]);
    }

    allreduce(pTotalMean, 1, std::plus<double>());

    pTotalMean /= initialTotalNum_m;

    if (std::abs(pTotalMean - referencePtot) / pTotalMean > 1e-2) {
        throw OpalException("ParallelCyclotronTracker::checkFileMomentum",
                            "The total momentum of the particle distribution\n"
                            "in the global reference frame: " +
                            std::to_string(pTotalMean) + ",\n"
                            "is different from the momentum given\n"
                            "in the \"BEAM\" command: " +
                            std::to_string(referencePtot) + ".\n"
                            "In Opal-cycl the initial distribution\n"
                            "is specified in the local reference frame.\n"
                            "When using a \"FROMFILE\" type distribution, the momentum \n"
                            "must be the same as the specified in the \"BEAM\" command,\n"
                            "which is in global reference frame.");
    }
}


//TODO: This can be completely rewritten with TracerParticle -DW
void ParallelCyclotronTracker::singleParticleDump() {
    IpplTimings::startTimer(DumpTimer_m);

    if (Ippl::getNodes() > 1 ) {

        double x;
        int id;
        dvector_t tmpr;
        ivector_t tmpi;

        int tag = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);

        // for all nodes, find the location of particle with ID = 0 & 1 in bunch containers
        int found[2] = {-1, -1};
        int counter = 0;

        for (size_t i = 0; i < itsBunch_m->getLocalNum(); ++i) {
            if (itsBunch_m->ID[i] == 0) {
                found[counter] = i;
                counter++;
            }
            if (itsBunch_m->ID[i] == 1) {
                found[counter] = i;
                counter++;
            }
        }

        if (myNode_m == 0) {
            int notReceived = Ippl::getNodes() - 1;
            int numberOfPart = 0;
            // receive other nodes
            while(notReceived > 0) {

                int node = COMM_ANY_NODE;
                Message *rmsg =  Ippl::Comm->receive_block(node, tag);

                if (rmsg == nullptr)
                    ERRORMSG("Could not receive from client nodes in main." << endl);

                --notReceived;

                rmsg->get(&numberOfPart);

                for (int i = 0; i < numberOfPart; ++i) {
                    rmsg->get(&id);
                    tmpi.push_back(id);
                    for (int ii = 0; ii < 6; ii++) {
                        rmsg->get(&x);
                        tmpr.push_back(x);
                    }
                }
                delete rmsg;
            }
            // own node
            for (int i = 0; i < counter; ++i) {

                tmpi.push_back(itsBunch_m->ID[found[i]]);

                for (int j = 0; j < 3; ++j) {
                    tmpr.push_back(itsBunch_m->R[found[i]](j));
                    tmpr.push_back(itsBunch_m->P[found[i]](j));
                }
            }
            // store
            dvector_t::iterator itParameter = tmpr.begin();

            for (auto tmpid : tmpi) {

                outfTrackOrbit_m << "ID" << tmpid;

                if (tmpid == 0) { // for stat file
                    itsBunch_m->RefPartR_m[0] = *itParameter;
                    itsBunch_m->RefPartR_m[1] = *(itParameter + 2);
                    itsBunch_m->RefPartR_m[2] = *(itParameter + 4);
                    itsBunch_m->RefPartP_m[0] = *(itParameter + 1);
                    itsBunch_m->RefPartP_m[1] = *(itParameter + 3);
                    itsBunch_m->RefPartP_m[2] = *(itParameter + 5);
                }
                for (int ii = 0; ii < 6; ii++) {
                    outfTrackOrbit_m << " " << *itParameter;
                    ++itParameter;
                }

                outfTrackOrbit_m << std::endl;
            }
        } else {
            // for other nodes
            Message *smsg = new Message();
            smsg->put(counter);

            for (int i = 0; i < counter; i++) {

                smsg->put(itsBunch_m->ID[found[i]]);

                for (int j = 0; j < 3; j++) {
                    smsg->put(itsBunch_m->R[found[i]](j));
                    smsg->put(itsBunch_m->P[found[i]](j));
                }
            }

            if (!Ippl::Comm->send(smsg, 0, tag)) {
                ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed " << endl);
            }
        }

        Ippl::Comm->barrier();

    } else {

        for (size_t i = 0; i < itsBunch_m->getLocalNum(); i++) {
            if (itsBunch_m->ID[i] == 0 || itsBunch_m->ID[i] == 1) {

                outfTrackOrbit_m << "ID" << itsBunch_m->ID[i] << " ";
                outfTrackOrbit_m << itsBunch_m->R[i](0) << " " << itsBunch_m->P[i](0) << " ";
                outfTrackOrbit_m << itsBunch_m->R[i](1) << " " << itsBunch_m->P[i](1) << " ";
                outfTrackOrbit_m << itsBunch_m->R[i](2) << " " << itsBunch_m->P[i](2) << std::endl;

                if (itsBunch_m->ID[i] == 0) { // for stat file
                    itsBunch_m->RefPartR_m = itsBunch_m->R[i];
                    itsBunch_m->RefPartP_m = itsBunch_m->P[i];
                }
            }
        }
    }

    IpplTimings::stopTimer(DumpTimer_m);
}

void ParallelCyclotronTracker::bunchDumpStatData(){

    IpplTimings::startTimer(DumpTimer_m);

    // dump stat file per bunch in case of multi-bunch mode
    if (isMultiBunch()) {
        double phi = 0.0, psi = 0.0;
        Vector_t meanR = calcMeanR();

        // Bunch (global) angle w.r.t. x-axis (cylinder coordinates)
        double theta = calculateAngle(meanR(0), meanR(1)) * Units::rad2deg;

        // fix azimuth_m --> make monotonically increasing
        dumpAngle(theta, prevAzimuth_m, azimuth_m);

        updateAzimuthAndRadius();

        if (Options::psDumpFrame != DumpFrame::GLOBAL) {
            Vector_t meanP = calcMeanP();

            // Bunch (local) azimuth at meanR w.r.t. y-axis
            phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * Physics::pi;

            // Bunch (local) elevation at meanR w.r.t. xy plane
            psi = 0.5 * Physics::pi - std::acos(meanP(2) / std::sqrt(dot(meanP, meanP)));

            // Rotate so Pmean is in positive y direction. No shift, so that normalized emittance and
            // unnormalized emittance as well as centroids are calculated correctly in
            // PartBunch::calcBeamParameters()
            globalToLocal(itsBunch_m->R, phi, psi, meanR);
            globalToLocal(itsBunch_m->P, phi, psi);
        }

        itsDataSink->writeMultiBunchStatistics(itsBunch_m, mbHandler_m.get());

        if (Options::psDumpFrame != DumpFrame::GLOBAL) {
            localToGlobal(itsBunch_m->R, phi, psi, meanR);
            localToGlobal(itsBunch_m->P, phi, psi);
        }

        IpplTimings::stopTimer(DumpTimer_m);
        return;
    }

    // --------------------------------- Get some Values ---------------------------------------- //
    double const temp_t = itsBunch_m->getT() * Units::s2ns;
    Vector_t meanR;
    Vector_t meanP;
    if (Options::psDumpFrame == DumpFrame::BUNCH_MEAN) {
        meanR = calcMeanR();
        meanP = calcMeanP();
    } else if (itsBunch_m->getLocalNum() > 0) {
        meanR = itsBunch_m->R[0];
        meanP = itsBunch_m->P[0];
    }
    double phi = 0;
    double psi = 0;

    // Bunch (global) angle w.r.t. x-axis (cylinder coordinates)
    double azimuth = calculateAngle(meanR(0), meanR(1)) * Units::rad2deg;

    // fix azimuth_m --> make monotonically increasing
    dumpAngle(azimuth, prevAzimuth_m, azimuth_m);

    // --------------  Calculate the external fields at the center of the bunch ----------------- //
    beamline_list::iterator DumpSindex = FieldDimensions.begin();

    extE_m = Vector_t(0.0, 0.0, 0.0);
    extB_m = Vector_t(0.0, 0.0, 0.0);

    (((*DumpSindex)->second).second)->apply(meanR, meanP, temp_t, extE_m, extB_m);

    // If we are saving in local frame, bunch and fields at the bunch center have to be rotated
    // TODO: Make decision if we maybe want to always save statistics data in local frame? -DW
    if (Options::psDumpFrame != DumpFrame::GLOBAL) {
        // -------------------- ----------- Do Transformations ---------------------------------- //
        // Bunch (local) azimuth at meanR w.r.t. y-axis
        phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * Physics::pi;

        // Bunch (local) elevation at meanR w.r.t. xy plane
        psi = 0.5 * Physics::pi - std::acos(meanP(2) / std::sqrt(dot(meanP, meanP)));

        // Rotate so Pmean is in positive y direction. No shift, so that normalized emittance and
        // unnormalized emittance as well as centroids are calculated correctly in
        // PartBunch::calcBeamParameters()
        globalToLocal(extB_m, phi, psi);
        globalToLocal(extE_m, phi, psi);
        globalToLocal(itsBunch_m->R, phi, psi);
        globalToLocal(itsBunch_m->P, phi, psi);
    }

    //itsBunch_m->R *= Vector_t(0.001); // mm -> m

    FDext_m[0] = extB_m * Units::kG2T;
    FDext_m[1] = extE_m;        // kV/mm? -DW

    // Save the stat file
    itsDataSink->dumpSDDS(itsBunch_m, FDext_m, azimuth_m);

    //itsBunch_m->R *= Vector_t(1000.0); // m -> mm

    // If we are in local mode, transform back after saving
    if (Options::psDumpFrame != DumpFrame::GLOBAL) {
        localToGlobal(itsBunch_m->R, phi, psi);
        localToGlobal(itsBunch_m->P, phi, psi);
    }

    IpplTimings::stopTimer(DumpTimer_m);
}


void ParallelCyclotronTracker::bunchDumpPhaseSpaceData() {
    // --------------------------------- Particle dumping --------------------------------------- //
    // Note: Don't dump when
    // 1. after one turn
    // in order to sychronize the dump step for multi-bunch and single bunch for comparison
    // with each other during post-process phase.
    // ------------------------------------------------------------------------------------------ //
    IpplTimings::startTimer(DumpTimer_m);

    // --------------------------------- Get some Values ---------------------------------------- //
    double const temp_t = itsBunch_m->getT() * Units::s2ns;

    Vector_t meanR;
    Vector_t meanP;

    // in case of multi-bunch mode take always bunch mean (although it takes all bunches)
    if (Options::psDumpFrame == DumpFrame::BUNCH_MEAN || isMultiBunch()) {
        meanR = calcMeanR();
        meanP = calcMeanP();
    } else if (itsBunch_m->getLocalNum() > 0) {
        meanR = itsBunch_m->R[0];
        meanP = itsBunch_m->P[0];
    }

    double const betagamma_temp = std::sqrt(dot(meanP, meanP));

    double const E = itsBunch_m->get_meanKineticEnergy();

    // Bunch (global) angle w.r.t. x-axis (cylinder coordinates)
    double const theta = std::atan2(meanR(1), meanR(0));

    // Bunch (local) azimuth at meanR w.r.t. y-axis
    double const phi = calculateAngle(meanP(0), meanP(1)) - 0.5 * Physics::pi;

    // Bunch (local) elevation at meanR w.r.t. xy plane
    double const psi = 0.5 * Physics::pi - std::acos(meanP(2) / std::sqrt(dot(meanP, meanP)));

    // ---------------- Re-calculate reference values in format of input values ----------------- //
    // Position:
    // New OPAL 2.0: Init in m (back to mm just for output) -DW
    referenceR = computeRadius(meanR); // includes m --> mm conversion
    referenceTheta = theta / Units::deg2rad;
    referenceZ = Units::m2mm * meanR(2);

    // Momentum in Theta-hat, R-hat, Z-hat coordinates at position meanR:
    referencePtot = betagamma_temp;
    referencePz = meanP(2);
    referencePr = meanP(0) * std::cos(theta) + meanP(1) * std::sin(theta);
    referencePt = std::sqrt(referencePtot * referencePtot - \
                       referencePz * referencePz - referencePr * referencePr);

    // -----  Calculate the external fields at the center of the bunch (Cave: Global Frame) ----- //
    beamline_list::iterator DumpSindex = FieldDimensions.begin();

    extE_m = Vector_t(0.0, 0.0, 0.0);
    extB_m = Vector_t(0.0, 0.0, 0.0);

    (((*DumpSindex)->second).second)->apply(meanR, meanP, temp_t, extE_m, extB_m);
    FDext_m[0] = extB_m * Units::kG2T;
    FDext_m[1] = extE_m;

    // -------------- If flag psDumpFrame is global, dump bunch in global frame ------------- //
    if (Options::psDumpFreq < std::numeric_limits<int>::max() ){
        bool dumpLocalFrame    = true;
        std::string dumpString = "local";
        if (Options::psDumpFrame == DumpFrame::GLOBAL) {
            dumpLocalFrame = false;
            dumpString     = "global";
        }

        if (dumpLocalFrame == true) {
            // ---------------- If flag psDumpFrame is local, dump bunch in local frame ---------------- //

            // The bunch is transformed into a local coordinate system with meanP in direction of y-axis
            globalToLocal(itsBunch_m->R, phi, psi, meanR);
            //globalToLocal(itsBunch_m->R, phi, psi, meanR * Vector_t(0.001));
            globalToLocal(itsBunch_m->P, phi, psi); // P should only be rotated

            globalToLocal(extB_m, phi, psi);
            globalToLocal(extE_m, phi, psi);
        }

        FDext_m[0] = extB_m * Units::kG2T;
        FDext_m[1] = extE_m;

        lastDumpedStep_m = itsDataSink->dumpH5(itsBunch_m, // Local and in m
                                               FDext_m, E,
                                               referencePr,
                                               referencePt,
                                               referencePz,
                                               referenceR,
                                               referenceTheta,
                                               referenceZ,
                                               phi / Units::deg2rad, // P_mean azimuth
                                               // at ref. R/Th/Z
                                               psi / Units::deg2rad, // P_mean elevation
                                               // at ref. R/Th/Z
                                               dumpLocalFrame);        // Flag localFrame

        if (dumpLocalFrame == true) {
            // Return to global frame
            localToGlobal(itsBunch_m->R, phi, psi, meanR);
            //localToGlobal(itsBunch_m->R, phi, psi, meanR * Vector_t(0.001));
            localToGlobal(itsBunch_m->P, phi, psi);
        }

        // Tell user in which mode we are dumping
        // New: no longer dumping for num part < 3, omit phase space dump number info
        if (lastDumpedStep_m == -1) {
          *gmsg << endl << "* Integration step " << step_m + 1
                << " (no phase space dump for <= 2 particles)" << endl;
        } else {
          *gmsg << endl << "* Phase space dump " << lastDumpedStep_m
                << " (" << dumpString << " frame) at integration step " << step_m + 1 << endl;
        }
    }

    // Print dump information on screen
    *gmsg << "* T = " << temp_t << " ns"
          << ", Live Particles: " << itsBunch_m->getTotalNum() << endl;
    *gmsg << "* E = " << E << " MeV"
          << ", beta * gamma = " << betagamma_temp << endl;
    *gmsg << "* Bunch position: R =  " << referenceR << " mm"
          << ", Theta = " << referenceTheta << " Deg"
          << ", Z = " << referenceZ << " mm" << endl;
    *gmsg << "* Local Azimuth = " << phi / Units::deg2rad << " Deg"
          << ", Local Elevation = " << psi / Units::deg2rad << " Deg" << endl;

    IpplTimings::stopTimer(DumpTimer_m);
}

bool ParallelCyclotronTracker::isTurnDone() {
    return (step_m > 10) && (((step_m + 1) %setup_m.stepsPerTurn) == 0);
}

void ParallelCyclotronTracker::update_m(double& t, const double& dt,
                                        const bool& finishedTurn)
{
    // Reference parameters
    t += dt;

    updateTime(dt);

    itsBunch_m->setLocalTrackStep((step_m + 1));
    if (!(step_m + 1 % 1000)) {
        *gmsg << "Step " << step_m + 1 << endl;
    }

    updatePathLength(dt);

    // Here is global frame, don't do itsBunch_m->boundp();

    if (itsBunch_m->getTotalNum()>0) {
        // Only dump last step if we have particles left.
        // Check separately for phase space (ps) and statistics (stat) data dump frequency
        if ( mode_m != TrackingMode::SEO && ( ((step_m + 1) % Options::psDumpFreq == 0) ||
                                      (Options::psDumpEachTurn && finishedTurn)))
        {
            // Write phase space data to h5 file
            bunchDumpPhaseSpaceData();
        }

        if ( mode_m != TrackingMode::SEO && ( ((step_m + 1) % Options::statDumpFreq == 0) ||
                                      (Options::psDumpEachTurn && finishedTurn)))
        {
            // Write statistics data to stat file
            bunchDumpStatData();
        }
    }

    if (Options::psDumpEachTurn && finishedTurn) {
        for (PluginElement* element : pluginElements_m) {
            element->save();
        }
    }
}


std::tuple<double, double, double> ParallelCyclotronTracker::initializeTracking_m() {
    // Read in some control parameters
    setup_m.scSolveFreq         = (spiral_flag) ? 1 : Options::scSolveFreq;
    setup_m.stepsPerTurn        = itsBunch_m->getStepsPerTurn();

    // Define 3 special azimuthal angles where dump particle's six parameters
    // at each turn into 3 ASCII files. only important in single particle tracking
    azimuth_angle_m.resize(3);
    azimuth_angle_m[0] = 0.0;
    azimuth_angle_m[1] = 22.5 * Units::deg2rad;
    azimuth_angle_m[2] = 45.0 * Units::deg2rad;

    double harm = getHarmonicNumber();
    double dt   = itsBunch_m->getdT() * Units::s2ns * harm;
    double t    = itsBunch_m->getT()  * Units::s2ns;

    double oldReferenceTheta = referenceTheta * Units::deg2rad;      // init here, reset each step
    setup_m.deltaTheta       = Physics::pi / (setup_m.stepsPerTurn); // half of the average angle per step

    //int stepToLastInj = itsBunch_m->getSteptoLastInj(); // TODO: Do we need this? -DW

    // Record how many bunches have already been injected. ONLY FOR MBM
    if (isMultiBunch()) {
        mbHandler_m->setNumBunch(itsBunch_m->getNumBunch());
     }

    initTrackOrbitFile();

    // Get data from h5 file for restart run and reset current step
    // to last step of previous simulation
    if (OpalData::getInstance()->inRestartRun()) {

        restartStep0_m = itsBunch_m->getLocalTrackStep();
        step_m = restartStep0_m;
        turnnumber_m = step_m / setup_m.stepsPerTurn + 1;

        *gmsg << "* Restart at integration step " << restartStep0_m
              << " at turn " << turnnumber_m - 1 << endl;

        initPathLength();
    }

    setup_m.stepsNextCheck = step_m + setup_m.stepsPerTurn; // Steps to next check for transition

    initDistInGlobalFrame();

    if (isMultiBunch()) {
        mbHandler_m->updateParticleBins(itsBunch_m);
    }

    // --- Output to user --- //
    *gmsg << "* Beginning of this run is at t = " << t << " [ns]" << endl;
    *gmsg << "* The time step is set to dt = " << dt << " [ns]" << endl;

    if ( isMultiBunch() ) {
        *gmsg << "* MBM: Time interval between neighbour bunches is set to "
              << setup_m.stepsPerTurn * dt << "[ns]" << endl;
        *gmsg << "* MBM: The particles energy bin reset frequency is set to "
              << Options::rebinFreq << endl;
    }

    *gmsg << "* Single particle trajectory dump frequency is set to " << Options::sptDumpFreq << endl;
    *gmsg << "* The frequency to solve space charge fields is set to " << setup_m.scSolveFreq << endl;
    *gmsg << "* The repartition frequency is set to " << Options::repartFreq << endl;

    switch ( mode_m ) {
        case TrackingMode::SEO: {
            *gmsg << endl;
            *gmsg << "* ------------------------- STATIC EQUILIBRIUM ORBIT MODE ----------------------------- *" << endl
                  << "* Instruction: When the total particle number is equal to 2, SEO mode is triggered      *" << endl
                  << "* automatically. This mode does NOT include any RF cavities. The initial distribution   *" << endl
                  << "* file must be specified. In the file the first line is for reference particle and the  *" << endl
                  << "* second line is for off-center particle. The tune is calculated by FFT routines based  *" << endl
                  << "* on these two particles.                                                               *" << endl
                  << "* ---------------- NOTE: SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE ------------------ *" << endl;

            if (Ippl::getNodes() != 1)
                throw OpalException("Error in ParallelCyclotronTracker::initializeTracking_m",
                                    "SEO MODE ONLY WORKS SERIALLY ON SINGLE NODE!");
            break;
        }
        case TrackingMode::SINGLE: {
            *gmsg << endl;
            *gmsg << "* ------------------------------ SINGLE PARTICLE MODE --------------------------------- *" << endl
                  << "* Instruction: When the total particle number is equal to 1, single particle mode is    *" << endl
                  << "* triggered automatically. The initial distribution file must be specified which should *" << endl
                  << "* contain only one line for the single particle                                         *" << endl
                  << "* ---------NOTE: SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON A SINGLE NODE ------------ *" << endl;

            if (Ippl::getNodes() != 1)
                throw OpalException("Error in ParallelCyclotronTracker::initializeTracking_m",
                                    "SINGLE PARTICLE MODE ONLY WORKS SERIALLY ON A SINGLE NODE!");

            // For single particle mode open output files
            openFiles(azimuth_angle_m.size() + 1, OpalData::getInstance()->getInputBasename());
            break;
        }
        case TrackingMode::BUNCH: {
            break;
        }
        case TrackingMode::UNDEFINED: {}
        default:  {
            throw OpalException("ParallelCyclotronTracker::initializeTracking_m()",
                                "No such tracking mode.");
        }
    }

    return std::make_tuple(t, dt, oldReferenceTheta);
}


void ParallelCyclotronTracker::finalizeTracking_m(dvector_t& Ttime,
                                                  dvector_t& Tdeltr,
                                                  dvector_t& Tdeltz, ivector_t& TturnNumber) {

    for (size_t ii = 0; ii < (itsBunch_m->getLocalNum()); ii++) {
        if (itsBunch_m->ID[ii] == 0) {
            double FinalMomentum2 = std::pow(itsBunch_m->P[ii](0), 2.0) + std::pow(itsBunch_m->P[ii](1), 2.0) + std::pow(itsBunch_m->P[ii](2), 2.0);
            double FinalEnergy = (std::sqrt(1.0 + FinalMomentum2) - 1.0) * itsBunch_m->getM() * Units::eV2MeV;
            *gmsg << "* Final energy of reference particle = " << FinalEnergy << " [MeV]" << endl;
            *gmsg << "* Total phase space dump number(includes the initial distribution) = " << lastDumpedStep_m + 1 << endl;
            *gmsg << "* One can restart simulation from the last dump step (--restart " << lastDumpedStep_m << ")" << endl;
        }
    }

    Ippl::Comm->barrier();

    switch ( mode_m ) {
        case TrackingMode::SEO: {
            // Calculate tunes after tracking.
            *gmsg << endl;
            *gmsg << "* **************** The result for tune calculation (NO space charge) ******************* *" << endl
                  << "* Number of tracked turns: " << TturnNumber.back() << endl;
            double nur, nuz;
            getTunes(Ttime, Tdeltr, Tdeltz, TturnNumber.back(), nur, nuz);
            break;
        }
        case TrackingMode::SINGLE:
            closeFiles();
            // fall through
        case TrackingMode::BUNCH: // we do nothing
        case TrackingMode::UNDEFINED:
        default: {
            // not for multibunch
            if (!isMultiBunch()) {
                *gmsg << "*" << endl;
                *gmsg << "* Finished during turn " << turnnumber_m << " (" << turnnumber_m - 1 << " turns completed)" << endl;
                *gmsg << "* Cave: Turn number is not correct for restart mode"<< endl;
            }
        }
    }

    Ippl::Comm->barrier();

    if (myNode_m == 0) {
        outfTrackOrbit_m.close();
    }

    *gmsg << endl << "* *********************** Bunch information in global frame: ***********************";

    if (itsBunch_m->getTotalNum() > 0) {
        // Print out the Bunch information at end of the run.
        itsBunch_m->calcBeamParameters();
        *gmsg << *itsBunch_m << endl;
    } else {
        *gmsg << endl << "* No Particles left in bunch!" << endl;
        *gmsg << "* **********************************************************************************" << endl;
    }
}

void ParallelCyclotronTracker::setTrackingMode() {
    if ( initialTotalNum_m == 1 ) {
        mode_m = TrackingMode::SINGLE;
    } else if ( initialTotalNum_m == 2 ) {
        mode_m = TrackingMode::SEO;
    } else if ( initialTotalNum_m > 2 ) {
        mode_m = TrackingMode::BUNCH;
    } else {
        mode_m = TrackingMode::UNDEFINED;
    }
}

void ParallelCyclotronTracker::seoMode_m(double& t, const double dt, bool& /*finishedTurn*/,
                                         dvector_t& Ttime, dvector_t& Tdeltr,
                                         dvector_t& Tdeltz, ivector_t& TturnNumber) {

    // 2 particles: Trigger SEO mode
    // (Switch off cavity and calculate betatron oscillation tuning)
    double r_tuning[2], z_tuning[2];

    IpplTimings::startTimer(IntegrationTimer_m);
    for (size_t i = 0; i < (itsBunch_m->getLocalNum()); i++) {

        if ((step_m % Options::sptDumpFreq == 0)) {
            outfTrackOrbit_m << "ID" << (itsBunch_m->ID[i]);
            outfTrackOrbit_m << " " << itsBunch_m->R[i](0)
                             << " " << itsBunch_m->P[i](0)
                             << " " << itsBunch_m->R[i](1)
                             << " " << itsBunch_m->P[i](1)
                             << " " << itsBunch_m->R[i](2)
                             << " " << itsBunch_m->P[i](2)
                             << std::endl;
        }

        double OldTheta = calculateAngle(itsBunch_m->R[i](0), itsBunch_m->R[i](1));
        r_tuning[i] = itsBunch_m->R[i](0) * std::cos(OldTheta) +
                      itsBunch_m->R[i](1) * std::sin(OldTheta);

        z_tuning[i] = itsBunch_m->R[i](2);

        // Integrate for one step in the lab Cartesian frame (absolute value).
        itsStepper_mp->advance(itsBunch_m, i, t, dt);

        if ( (i == 0) && isTurnDone() ) {
            turnnumber_m++;
        }

    } // end for: finish one step tracking for all particles for initialTotalNum_m == 2 mode
    IpplTimings::stopTimer(IntegrationTimer_m);

    // store dx and dz for future tune calculation if higher precision needed, reduce freqSample.
    if (step_m % Options::sptDumpFreq == 0) {
        Ttime.push_back(t * Units::ns2s);
        Tdeltz.push_back(z_tuning[1]);
        Tdeltr.push_back(r_tuning[1] - r_tuning[0]);
        TturnNumber.push_back(turnnumber_m);
    }
}


void ParallelCyclotronTracker::singleMode_m(double& t, const double dt,
                                            bool& finishedTurn, double& oldReferenceTheta) {
    // 1 particle: Trigger single particle mode

    // ********************************************************************************** //
    // * This was moved here b/c collision should be tested before the actual           * //
    // * timestep (bgf_main_collision_test() predicts the next step automatically)      * //

    // apply the plugin elements: probe, collimator, stripper, septum
    bool flagNeedUpdate = applyPluginElements(dt);

    // check if we lose particles at the boundary
    bgf_main_collision_test();
    // ********************************************************************************** //

    if (itsBunch_m->getLocalNum() == 0) return; // might happen if particle is in collimator

    IpplTimings::startTimer(IntegrationTimer_m);

    unsigned int i = 0; // we only have a single particle

    if ( step_m % Options::sptDumpFreq == 0 ) {
        outfTrackOrbit_m << "ID" <<itsBunch_m->ID[i]
                         << " " << itsBunch_m->R[i](0)
                         << " " << itsBunch_m->P[i](0)
                         << " " << itsBunch_m->R[i](1)
                         << " " << itsBunch_m->P[i](1)
                         << " " << itsBunch_m->R[i](2)
                         << " " << itsBunch_m->P[i](2)
                         << std::endl;
    }

    double temp_meanTheta = calculateAngle2(itsBunch_m->R[i](0),
                                            itsBunch_m->R[i](1)); // [-pi, pi]

    dumpThetaEachTurn_m(t, itsBunch_m->R[i], itsBunch_m->P[i],
                        temp_meanTheta, finishedTurn);

    dumpAzimuthAngles_m(t, itsBunch_m->R[i], itsBunch_m->P[i],
                        oldReferenceTheta, temp_meanTheta);

    oldReferenceTheta = temp_meanTheta;

    // used for gap crossing checking
    Vector_t Rold = itsBunch_m->R[i]; // [x,y,z]    (mm)
    Vector_t Pold = itsBunch_m->P[i]; // [px,py,pz] (beta*gamma)

    // integrate for one step in the lab Cartesian frame (absolute value).
    itsStepper_mp->advance(itsBunch_m, i, t, dt);

    // If gap crossing happens, do momenta kicking (not if gap crossing just happened)
    if (itsBunch_m->cavityGapCrossed[i] == true)
        itsBunch_m->cavityGapCrossed[i] = false;
    else
        gapCrossKick_m(i, t, dt, Rold, Pold);
    IpplTimings::stopTimer(IntegrationTimer_m);

    // Destroy particles if they are marked as Bin = -1 in the plugin elements
    // or out of the global aperture
    flagNeedUpdate |= (itsBunch_m->Bin[i] < 0);
    deleteParticle(flagNeedUpdate);
}


void ParallelCyclotronTracker::bunchMode_m(double& t, const double dt, bool& finishedTurn) {

    // Flag for transition from single-bunch to multi-bunches mode
    static bool flagTransition = false;

    // single particle dumping
    if (step_m % Options::sptDumpFreq == 0)
        singleParticleDump();

    // Find out if we need to do bunch injection
    injectBunch(flagTransition);

//     oldReferenceTheta = calculateAngle2(meanR(0), meanR(1));

    // Calculate SC field before each time step and keep constant during integration.
    // Space Charge effects are included only when total macropaticles number is NOT LESS THAN 1000.
    if (itsBunch_m->hasFieldSolver()) {

        if (step_m % setup_m.scSolveFreq == 0) {
            computeSpaceChargeFields_m();
        } else {
            // If we are not solving for the space charge fields at this time step
            // we will apply the fields from the previous step and have to rotate them
            // accordingly. For this we find the quaternion between the previous mean momentum (PreviousMeanP)
            // and the current mean momentum (meanP) and rotate the fields with this quaternion.

            Vector_t const meanP = calcMeanP();

            Quaternion_t quaternionToNewMeanP;

            getQuaternionTwoVectors(PreviousMeanP, meanP, quaternionToNewMeanP);

            // Reset PreviousMeanP. Cave: This HAS to be after the quaternion is calculated!
            PreviousMeanP = calcMeanP();

            // Rotate the fields
            globalToLocal(itsBunch_m->Ef, quaternionToNewMeanP);
            globalToLocal(itsBunch_m->Bf, quaternionToNewMeanP);
        }
    }

    // *** This was moved here b/c collision should be tested before the **********************
    // *** actual timestep (bgf_main_collision_test() predicts the next step automatically) -DW
    // Apply the plugin elements: probe, collimator, stripper, septum
    bool flagNeedUpdate = applyPluginElements(dt);

    // check if we lose particles at the boundary
    bgf_main_collision_test();

    IpplTimings::startTimer(IntegrationTimer_m);

    for (size_t i = 0; i < itsBunch_m->getLocalNum(); i++) {

        // used for gap crossing checking
        Vector_t Rold = itsBunch_m->R[i]; // [x,y,z]    (mm)
        Vector_t Pold = itsBunch_m->P[i]; // [px,py,pz] (beta*gamma)

        // Integrate for one step in the lab Cartesian frame (absolute value).
        itsStepper_mp->advance(itsBunch_m, i, t, dt);

        // If gap crossing happens, do momenta kicking (not if gap crossing just happened)
        if (itsBunch_m->cavityGapCrossed[i] == true) {
            itsBunch_m->cavityGapCrossed[i] = false;
        } else {
            gapCrossKick_m(i, t, dt, Rold, Pold);
        }
        flagNeedUpdate |= (itsBunch_m->Bin[i] < 0);
    }

    IpplTimings::stopTimer(IntegrationTimer_m);

    // Destroy particles if they are marked as Bin = -1 in the plugin elements
    // or out of global aperture
    deleteParticle(flagNeedUpdate);

    // Recalculate bingamma and reset the BinID for each particles according to its current gamma
    if (isMultiBunch() && step_m % Options::rebinFreq == 0) {
        mbHandler_m->updateParticleBins(itsBunch_m);
    }

    // Some status output for user after each turn
    if ( isTurnDone() ) {
        turnnumber_m++;
        finishedTurn = true;

        *gmsg << endl;
        *gmsg << "*** Finished turn " << turnnumber_m - 1
              << ", Total number of live particles: "
              << itsBunch_m->getTotalNum() << endl;
    }

    Ippl::Comm->barrier();
}


void ParallelCyclotronTracker::gapCrossKick_m(size_t i, double t,
                                              double dt,
                                              const Vector_t& Rold,
                                              const Vector_t& Pold) {

    for (beamline_list::iterator sindex = ++(FieldDimensions.begin());
        sindex != FieldDimensions.end(); ++sindex)
    {
        bool tag_crossing = false;
        double DistOld = 0.0; //mm
        RFCavity * rfcav;

        if (((*sindex)->first) == ElementType::RFCAVITY) {
            // here check gap cross in the list, if do , set tag_crossing to TRUE
            rfcav = static_cast<RFCavity *>(((*sindex)->second).second);
            tag_crossing = checkGapCross(Rold, itsBunch_m->R[i], rfcav, DistOld);
        }

        if ( tag_crossing ) {
            itsBunch_m->cavityGapCrossed[i] = true;

            double oldMomentum2  = dot(Pold, Pold);
            double oldBetgam = std::sqrt(oldMomentum2);
            double oldGamma = std::sqrt(1.0 + oldMomentum2);
            double oldBeta = oldBetgam / oldGamma;
            double dt1 = DistOld / (Physics::c * oldBeta * Units::m2mm / Units::s2ns); // ns, c in [mm/ns]
            double dt2 = dt - dt1;

            // retrack particle from the old postion to cavity gap point
            // restore the old coordinates and momenta
            itsBunch_m->R[i] = Rold;
            itsBunch_m->P[i] = Pold;

            if (dt / dt1 < 1.0e9) {
                itsStepper_mp->advance(itsBunch_m, i, t, dt1);
            }
            // Momentum kick
            RFkick(rfcav, t, dt1, i);

            /* Retrack particle from cavity gap point for
             * the left time to finish the entire timestep
             */
            if (dt / dt2 < 1.0e9) {
                itsStepper_mp->advance(itsBunch_m, i, t, dt2);
            }
        }
    }
}


void ParallelCyclotronTracker::dumpAzimuthAngles_m(const double& t,
                                                   const Vector_t& R,
                                                   const Vector_t& P,
                                                   const double& oldReferenceTheta,
                                                   const double& temp_meanTheta) {

    for (unsigned int i=0; i<=2; i++) {
        if ((oldReferenceTheta < azimuth_angle_m[i] - setup_m.deltaTheta) &&
            (  temp_meanTheta >= azimuth_angle_m[i] - setup_m.deltaTheta))
        {
            *(outfTheta_m[i]) << "#Turn number = " << turnnumber_m
                              << ", Time = " << t << " [ns]" << std::endl
                              << " " << std::hypot(R(0), R(1))
                              << " " << P(0) * std::cos(temp_meanTheta) + P(1) * std::sin(temp_meanTheta)
                              << " " << temp_meanTheta * Units::rad2deg
                              << " " << - P(0) * std::sin(temp_meanTheta) + P(1) * std::cos(temp_meanTheta)
                              << " " << R(2)
                              << " " << P(2) << std::endl;
        }
    }
}

void ParallelCyclotronTracker::dumpThetaEachTurn_m(const double& t,
                                                   const Vector_t& R,
                                                   const Vector_t& P,
                                                   const double& temp_meanTheta,
                                                   bool& finishedTurn) {
    if ( isTurnDone() ) {
        ++turnnumber_m;
        finishedTurn = true;

        *gmsg << "* SPT: Finished turn " << turnnumber_m - 1 << endl;

        *(outfTheta_m[3]) << "#Turn number = " << turnnumber_m
                          << ", Time = " << t << " [ns]" << std::endl
                          << " " << std::sqrt(R(0) * R(0) + R(1) * R(1))
                          << " " << P(0) * std::cos(temp_meanTheta) +
                                    P(1) * std::sin(temp_meanTheta)
                          << " " << temp_meanTheta / Units::deg2rad
                          << " " << - P(0) * std::sin(temp_meanTheta) +
                                      P(1) * std::cos(temp_meanTheta)
                          << " " << R(2)
                          << " " << P(2) << std::endl;
    }
}


void ParallelCyclotronTracker::computeSpaceChargeFields_m() {
    // Firstly reset E and B to zero before fill new space charge field data for each track step
    itsBunch_m->Bf = Vector_t(0.0);
    itsBunch_m->Ef = Vector_t(0.0);

    if (spiral_flag && itsBunch_m->getFieldSolverType() == FieldSolverType::SAAMG) {
        // --- Single bunch mode with spiral inflector --- //

        // If we are doing a spiral inflector simulation and are using the SAAMG solver
        // we don't rotate or translate the bunch and gamma is 1.0 (non-relativistic).

        //itsBunch_m->R *= Vector_t(0.001); // mm --> m

        itsBunch_m->setGlobalMeanR(Vector_t(0.0, 0.0, 0.0));
        itsBunch_m->setGlobalToLocalQuaternion(Quaternion_t(1.0, 0.0, 0.0, 0.0));

        itsBunch_m->computeSelfFields_cycl(1.0);

        //itsBunch_m->R *= Vector_t(1000.0); // m --> mm

    } else {

        Vector_t const meanR = calcMeanR(); // (m)

        PreviousMeanP = calcMeanP();

        // Since calcMeanP takes into account all particles of all bins (TODO: Check this! -DW)
        // Using the quaternion method with PreviousMeanP and yaxis should give the correct result

        Quaternion_t quaternionToYAxis;

        getQuaternionTwoVectors(PreviousMeanP, yaxis, quaternionToYAxis);

        globalToLocal(itsBunch_m->R, quaternionToYAxis, meanR);

        //itsBunch_m->R *= Vector_t(0.001); // mm --> m

        if ((step_m + 1) % Options::boundpDestroyFreq == 0) {
            itsBunch_m->boundp_destroyCycl();
        } else {
            itsBunch_m->boundp();
        }

        if (hasMultiBunch()) {
            // --- Multibunch mode --- //

            // Calculate gamma for each energy bin
            itsBunch_m->calcGammas_cycl();

            repartition();

            // Calculate space charge field for each energy bin
            for (int b = 0; b < itsBunch_m->getLastemittedBin(); b++) {

                itsBunch_m->setBinCharge(b, itsBunch_m->getChargePerParticle());
                //itsBunch_m->setGlobalMeanR(0.001 * meanR);
                itsBunch_m->setGlobalMeanR(meanR);
                itsBunch_m->setGlobalToLocalQuaternion(quaternionToYAxis);
                itsBunch_m->computeSelfFields_cycl(b);
            }

            itsBunch_m->Q = itsBunch_m->getChargePerParticle();

        } else {
            // --- Single bunch mode --- //

            double temp_meangamma = std::sqrt(1.0 + dot(PreviousMeanP, PreviousMeanP));

            repartition();

            //itsBunch_m->setGlobalMeanR(0.001 * meanR);
            itsBunch_m->setGlobalMeanR(meanR);
            itsBunch_m->setGlobalToLocalQuaternion(quaternionToYAxis);

            itsBunch_m->computeSelfFields_cycl(temp_meangamma);
        }

        //scale coordinates back
        //itsBunch_m->R *= Vector_t(1000.0); // m --> mm

        // Transform coordinates back to global
        localToGlobal(itsBunch_m->R, quaternionToYAxis, meanR);

        // Transform self field back to global frame (rotate only)
        localToGlobal(itsBunch_m->Ef, quaternionToYAxis);
        localToGlobal(itsBunch_m->Bf, quaternionToYAxis);
    }
}


bool ParallelCyclotronTracker::computeExternalFields_m(const size_t& i, const double& t,
                                                       Vector_t& Efield, Vector_t& Bfield) {

    beamline_list::iterator sindex = FieldDimensions.begin();

    // Flag whether a particle is out of field
    bool outOfBound = (((*sindex)->second).second)->apply(i, t, Efield, Bfield);

    Bfield *= Units::kG2T;
    Efield *= Units::kV2V / Units::mm2m;

    return outOfBound;
}

bool ParallelCyclotronTracker::computeExternalFields_m(const Vector_t& R, const Vector_t& P, const double& t,
                                                       Vector_t& Efield, Vector_t& Bfield) {

    beamline_list::iterator sindex = FieldDimensions.begin();
    // Flag whether a particle is out of field
    bool outOfBound = (((*sindex)->second).second)->apply(R, P, t, Efield, Bfield);

    Bfield *= Units::kG2T;
    Efield *= Units::kV2V / Units::mm2m;

    return outOfBound;
}



void ParallelCyclotronTracker::injectBunch(bool& flagTransition) {
    if (!isMultiBunch() || step_m != setup_m.stepsNextCheck) {
        return;
    }

    const short result = mbHandler_m->injectBunch(itsBunch_m, itsReference,
                                                  flagTransition);

    switch ( result ) {
        case 0: {
            // nothing happened
            break;
        }
        case 1: {
            // bunch got saved
            saveInjectValues();
            setup_m.stepsNextCheck += setup_m.stepsPerTurn;
            if ( flagTransition ) {
                *gmsg << "* MBM: Saving beam distribution at turn " << turnnumber_m << endl;
                *gmsg << "* MBM: After one revolution, Multi-Bunch Mode will be invoked" << endl;
            }
            break;
        }
        case 2: {
            // bunch got injected
            setup_m.stepsNextCheck += setup_m.stepsPerTurn;
            break;
        }
        default: {
            throw OpalException("ParallelCyclotronTracker::injectBunch()",
                                "Unknown return value " + std::to_string(result));
        }
    }
}


void ParallelCyclotronTracker::saveInjectValues() {
    if ( !isMultiBunch() ) {
        return;
    }

    Vector_t meanR = calcMeanR();

    // Bunch (global) angle w.r.t. x-axis (cylinder coordinates)
    double theta = calculateAngle(meanR(0), meanR(1)) * Units::rad2deg;

    // fix azimuth_m --> make monotonically increasing
    dumpAngle(theta, prevAzimuth_m, azimuth_m);

    const double radius = computeRadius(meanR);

    MultiBunchHandler::injection_t& inj = mbHandler_m->getInjectionValues();

    inj.time       = itsBunch_m->getT() * Units::s2ns;
    inj.pathlength = itsBunch_m->get_sPos();
    inj.azimuth    = azimuth_m;
    inj.radius     = radius;
}


void ParallelCyclotronTracker::updatePathLength(const double& dt) {
    /* the last element includes all particles,
     * all other elements are bunch number specific
     */
    std::vector<double> lpaths(1);

    if ( isMultiBunch() ) {
        lpaths.resize(mbHandler_m->getNumBunch() + 1);
    }

    computePathLengthUpdate(lpaths, dt);

    pathLength_m += lpaths.back();
    itsBunch_m->set_sPos(pathLength_m);

    if ( isMultiBunch() ) {
        mbHandler_m->updatePathLength(lpaths);
    }
}


void ParallelCyclotronTracker::updateTime(const double& dt) {
    // t is in seconds
    double t = itsBunch_m->getT();

    itsBunch_m->setT(t + dt * Units::ns2s);

    if ( isMultiBunch() ) {
        mbHandler_m->updateTime(dt);
    }
}


void ParallelCyclotronTracker::updateAzimuthAndRadius() {
    if (!isMultiBunch()) {
        return;
    }

    for (short b = 0; b < mbHandler_m->getNumBunch(); ++b) {
        Vector_t meanR = calcMeanR(b);
        MultiBunchHandler::beaminfo_t& binfo = mbHandler_m->getBunchInfo(b);

        binfo.radius  = computeRadius(meanR);
        double azimuth = calculateAngle(meanR(0), meanR(1)) * Units::rad2deg;
        dumpAngle(azimuth, binfo.prevAzimuth, binfo.azimuth, b);
    }
}


void ParallelCyclotronTracker::initPathLength() {
    if ( isMultiBunch() ) {
        // we need to reset the path length of each bunch
        itsDataSink->setMultiBunchInitialPathLengh(mbHandler_m.get());
    }
}
