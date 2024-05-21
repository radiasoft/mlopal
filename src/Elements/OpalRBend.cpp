//
// Class OpalRBend
//   The RBEND element.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Elements/OpalRBend.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/RBendRep.h"
#include "Fields/BMultipoleField.h"
#include "Physics/Physics.h"
#include "Structure/OpalWake.h"
#include "Structure/ParticleMatterInteraction.h"
#include "Utilities/OpalException.h"

OpalRBend::OpalRBend():
    OpalBend("RBEND",
             "The \"RBEND\" element defines a rectangular bending magnet."),
    owk_m(0),
    parmatint_m(nullptr) {

    registerOwnership();

    setElement(new RBendRep("RBEND"));
}


OpalRBend::OpalRBend(const std::string& name, OpalRBend* parent):
    OpalBend(name, parent),
    owk_m(0),
    parmatint_m(nullptr) {
    setElement(new RBendRep(name));
}


OpalRBend::~OpalRBend() {
    delete owk_m;
    delete parmatint_m;
}


OpalRBend* OpalRBend::clone(const std::string& name) {
    return new OpalRBend(name, this);
}


void OpalRBend::update() {
    OpalElement::update();

    // Define geometry.
    RBendRep* bend =
        dynamic_cast<RBendRep*>(getElement());

    double length = Attributes::getReal(itsAttr[LENGTH]);
    double angle  = Attributes::getReal(itsAttr[ANGLE]);
    double e1     = Attributes::getReal(itsAttr[E1]);
    RBendGeometry& geometry = bend->getGeometry();
    geometry.setElementLength(length);
    geometry.setBendAngle(angle);

    // Define number of slices for map tracking
    bend->setNSlices(Attributes::getReal(itsAttr[NSLICES]));

    // Define pole face angles.
    bend->setEntryFaceRotation(Attributes::getReal(itsAttr[E1]));
    bend->setExitFaceRotation(Attributes::getReal(itsAttr[E2]));
    bend->setEntryFaceCurvature(Attributes::getReal(itsAttr[H1]));
    bend->setExitFaceCurvature(Attributes::getReal(itsAttr[H2]));

    // Define integration parameters.
    bend->setSlices(Attributes::getReal(itsAttr[SLICES]));
    bend->setStepsize(Attributes::getReal(itsAttr[STEPSIZE]));

    // Define field.
    double factor = OpalData::getInstance()->getP0() / Physics::c;
    BMultipoleField field;
    double k0 =
        itsAttr[K0] ? Attributes::getReal(itsAttr[K0]) :
        length ? 2 * sin(angle / 2) / length : angle;
    double k0s = itsAttr[K0S] ? Attributes::getReal(itsAttr[K0S]) : 0.0;
    //JMJ 4/10/2000: above line replaced
    //    length ? angle / length : angle;
    // to avoid closed orbit created by RBEND with default K0.
    field.setNormalComponent(1, factor * k0);
    field.setSkewComponent(1, factor * Attributes::getReal(itsAttr[K0S]));
    field.setNormalComponent(2, factor * Attributes::getReal(itsAttr[K1]));
    field.setSkewComponent(2, factor * Attributes::getReal(itsAttr[K1S]));
    field.setNormalComponent(3, factor * Attributes::getReal(itsAttr[K2]) / 2.0);
    field.setSkewComponent(3, factor * Attributes::getReal(itsAttr[K2S]) / 2.0);
    field.setNormalComponent(4, factor * Attributes::getReal(itsAttr[K3]) / 6.0);
    field.setSkewComponent(4, factor * Attributes::getReal(itsAttr[K3S]) / 6.0);
    bend->setField(field);

    // Set field amplitude or bend angle.
    if (itsAttr[ANGLE]) {
        if (bend->isPositioned() && angle < 0.0) {
            e1 = -e1;
            angle = -angle;

            Quaternion rotAboutZ(0, 0, 0, 1);
            CoordinateSystemTrafo g2l = bend->getCSTrafoGlobal2Local();
            bend->releasePosition();
            bend->setCSTrafoGlobal2Local(CoordinateSystemTrafo(g2l.getOrigin(),
                                                               rotAboutZ * g2l.getRotation()));
            bend->fixPosition();
        }
        bend->setBendAngle(angle);
    } else {
        bend->setFieldAmplitude(k0, k0s);
    }
        bend->setEntranceAngle(e1);

    if (itsAttr[ROTATION])
        throw OpalException("OpalRBend::update",
                            "ROTATION not supported any more; use PSI instead");


    if (itsAttr[FMAPFN]) {
        bend->setFieldMapFN(Attributes::getString(itsAttr[FMAPFN]));
    } else if (bend->getName() != "RBEND") {
        ERRORMSG(bend->getName() << ": No filename for a field map given. "
                 "Will assume the default map "
                 "\"1DPROFILE1-DEFAULT\"."
                 << endl);
        bend->setFieldMapFN("1DPROFILE1-DEFAULT");
    }

    // Energy in eV.
    if (itsAttr[DESIGNENERGY] && Attributes::getReal(itsAttr[DESIGNENERGY]) != 0.0) {
        bend->setDesignEnergy(Attributes::getReal(itsAttr[DESIGNENERGY]), false);
    } else if (bend->getName() != "RBEND") {
        throw OpalException("OpalRBend::update",
                            "RBend requires non-zero DESIGNENERGY");
    }

    double gap = Attributes::getReal(itsAttr[GAP]);
    bend->setFullGap(gap);

    if (itsAttr[APERT])
        throw OpalException("OpalRBend::update",
                            "APERTURE in RBEND not supported; use GAP and HAPERT instead");

    if (itsAttr[HAPERT]) {
        double hapert = Attributes::getReal(itsAttr[HAPERT]);
        bend->setAperture(ApertureType::RECTANGULAR, std::vector<double>({hapert, gap, 1.0}));
    } else {
        bend->setAperture(ApertureType::RECTANGULAR, std::vector<double>({0.5, gap, 1.0}));
    }

    if (itsAttr[WAKEF] && itsAttr[DESIGNENERGY] && owk_m == nullptr) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*bend);
        bend->setWake(owk_m->wf_m);
    }

    if (itsAttr[K1])
        bend->setK1(Attributes::getReal(itsAttr[K1]));
    else
        bend->setK1(0.0);

    if (itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == nullptr) {
        const std::string matterDescriptor = Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION]);
        ParticleMatterInteraction* orig = ParticleMatterInteraction::find(matterDescriptor);
        parmatint_m = orig->clone(matterDescriptor);
        parmatint_m->initParticleMatterInteractionHandler(*bend);
        bend->setParticleMatterInteraction(parmatint_m->handler_m);
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(bend);
}