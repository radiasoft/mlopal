//
// Class CCollimator
//   Interface for cyclotron collimator
//
// Copyright (c) 2018-2022, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "AbsBeamline/CCollimator.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Fields/Fieldmap.h"
#include "Solvers/ParticleMatterInteractionHandler.h"
#include "Structure/LossDataSink.h"

#include <cmath>
#include <fstream>

extern Inform *gmsg;

CCollimator::CCollimator():CCollimator("")
{}

CCollimator::CCollimator(const std::string& name):
    PluginElement(name) {
    setDimensions(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    setGeom(0.0);
}

CCollimator::CCollimator(const CCollimator& right):
    PluginElement(right),
    informed_m(right.informed_m) {
    setDimensions(right.xstart_m, right.xend_m,
                  right.ystart_m, right.yend_m,
                  right.zstart_m, right.zend_m,
                  right.width_m);
    setGeom(width_m);
}

CCollimator::~CCollimator() {}

void CCollimator::accept(BeamlineVisitor &visitor) const {
    visitor.visitCCollimator(*this);
}

bool CCollimator::doPreCheck(PartBunchBase<double, 3>* bunch) {
    Vector_t rmin, rmax;
    bunch->get_bounds(rmin, rmax);

    if (rmax(2) >= zstart_m && rmin(2) <= zend_m) {
        // interested in absolute minimum and maximum
        double xmin = std::min(std::abs(rmin(0)), std::abs(rmax(0)));
        double xmax = std::max(std::abs(rmin(0)), std::abs(rmax(0)));
        double ymin = std::min(std::abs(rmin(1)), std::abs(rmax(1)));
        double ymax = std::max(std::abs(rmin(1)), std::abs(rmax(1)));
        double rbunch_min = std::hypot(xmin, ymin);
        double rbunch_max = std::hypot(xmax, ymax);

        if ( rbunch_max > rmin_m && rbunch_min < rmax_m ){ // check similar to z
            return true;
        }
    }
    return false;
}

// rectangle collimators in cyclotron cylindrical coordinates
// when there is no particlematterinteraction, the particle hitting collimator is deleted directly
bool CCollimator::doCheck(PartBunchBase<double, 3>* bunch, const int turnnumber,
                          const double t, const double /*tstep*/) {

    bool flagNeedUpdate = false;
    size_t tempnum = bunch->getLocalNum();
    int pflag = 0;
    // now check each particle in bunch
    for (unsigned int i = 0; i < tempnum; ++i) {
        if (bunch->R[i](2) < zend_m && bunch->R[i](2) > zstart_m ) {
            // only now careful check in r
            pflag = checkPoint(bunch->R[i](0), bunch->R[i](1));
            /// bunch->Bin[i] != -1 makes sure the particle is not stored in more than one collimator
            if ((pflag != 0) && (bunch->Bin[i] != -1)) {
                lossDs_m->addParticle(OpalParticle(bunch->ID[i],
                                                   bunch->R[i], bunch->P[i],
                                                   t, bunch->Q[i], bunch->M[i]),
                                      std::make_pair(turnnumber, bunch->bunchNum[i]));

                bunch->Bin[i] = -1;
                flagNeedUpdate = true;
            }
        }
    }
    return flagNeedUpdate;
}

bool CCollimator::doFinaliseCheck(PartBunchBase<double, 3>* bunch, bool flagNeedUpdate) {

    reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());

    if (flagNeedUpdate && parmatint_m) {
        *gmsg << level2 << "============== START PARTICLE MATTER INTERACTION CALCULATION =============" << endl;
        do {
            parmatint_m->setFlagAllParticlesIn(false);

            bool collWithParticles = (parmatint_m->getParticlesInMat() > 0);

            unsigned int localNum = bunch->getLocalNum();
            unsigned int totalNum = 0;
            reduce(localNum, totalNum, OpAddAssign());

            bool allParticlesInMat = (totalNum == 0 && collWithParticles);
            if (allParticlesInMat) {
                parmatint_m->setFlagAllParticlesIn(true);
            }

            Vector_t rmin, rmax;
            bunch->get_bounds(rmin, rmax);
            std::pair<Vector_t, double> boundingSphere;
            boundingSphere.first = 0.5 * (rmax + rmin);
            boundingSphere.second = euclidean_norm(rmax - boundingSphere.first);

            parmatint_m->apply(bunch, boundingSphere);
            parmatint_m->print(*gmsg);

        } while (parmatint_m->stillActive());

        *gmsg << level2 << "============== END PARTICLE MATTER INTERACTION CALCULATION =============" << endl;
    }
    return flagNeedUpdate;
}

void CCollimator::doInitialise(PartBunchBase<double, 3>* /*bunch*/) {
    parmatint_m = getParticleMatterInteraction();
}

void CCollimator::goOnline(const double&) {
    print();
    online_m = true;
}

void CCollimator::doFinalise() {
    *gmsg << "* Finalize cyclotron collimator " << getName() << endl;
}

void CCollimator::print() {
    if (RefPartBunch_m == nullptr) {
        if (!informed_m) {
            std::string errormsg = Fieldmap::typeset_msg("BUNCH SIZE NOT SET", "warning");
            ERRORMSG(errormsg << endl);
            if (Ippl::myNode() == 0) {
                std::ofstream omsg("errormsg.txt", std::ios_base::app);
                omsg << errormsg << std::endl;
                omsg.close();
            }
            informed_m = true;
        }
        return;
    }
}

void CCollimator::setDimensions(double xstart, double xend,
                                double ystart, double yend,
                                double zstart, double zend,
                                double width) {
    setDimensions(xstart, xend, ystart, yend);
    zstart_m = zstart;
    zend_m   = zend;
    width_m  = width;
    // zstart and zend are independent from x, y
    if (zstart_m > zend_m) {
        std::swap(zstart_m, zend_m);
    }
    setGeom(width_m);
}

void CCollimator::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = 0.0;
    zEnd = getElementLength();
}

ElementType CCollimator::getType() const {
    return ElementType::CCOLLIMATOR;
}

void CCollimator::doSetGeom() {
    // calculate maximum r from these
    rmax_m = -std::numeric_limits<double>::max();
    for (int i=0; i<4; i++) {
        double r = std::hypot(geom_m[i].x, geom_m[i].y);
        rmax_m = std::max(rmax_m, r);
    }
    // some debug printout
    // for (int i=0; i<4; i++) {
    //     *gmsg << "point " << i << " ( " << geom_m[i].x << ", " << geom_m[i].y << ")" << endl;
    // }
    // *gmsg << "rmin " << rmin_m << " rmax " << rmax_m << endl;
}
