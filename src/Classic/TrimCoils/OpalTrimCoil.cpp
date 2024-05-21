//
// Class OpalTrimCoil
//   A TRIMCOIL definition is used to define a trim coil which can be applied
//   to a Cyclotron.
//
// Copyright (c) 2018 - 2019, Matthias Frey and Jochem Snuverink,
//                            Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
// and the paper
// "Matching of turn pattern measurements for cyclotrons using multiobjective optimization"
// (https://doi.org/10.1103/PhysRevAccelBeams.22.064602)
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
#include "TrimCoils/OpalTrimCoil.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "TrimCoils/TrimCoilBFit.h"
#include "TrimCoils/TrimCoilPhaseFit.h"
#include "TrimCoils/TrimCoilMirrored.h"
#include "Utilities/OpalException.h"
#include "Utility/Inform.h"

extern Inform *gmsg;

// The attributes of class OpalTrimCoil.
namespace {
    enum {
        TYPE,       // The type of trim coil
        COEFNUM,    //
        COEFDENOM,  //
        COEFNUMPHI,
        COEFDENOMPHI,
        BMAX,       //
        PHIMIN,
        PHIMAX,
        RMIN,       //
        RMAX,       //
        SLPTC,
        SIZE
    };
}

OpalTrimCoil::OpalTrimCoil():
    Definition(SIZE, "TRIMCOIL",
               "The \"TRIMCOIL\" statement defines a trim coil."),
    trimcoil_m(nullptr) {
    itsAttr[TYPE]      = Attributes::makePredefinedString
                        ("TYPE", "Specifies the type of trim coil.", {"PSI-BFIELD", "PSI-PHASE", "PSI-BFIELD-MIRRORED"});

    itsAttr[COEFNUM]   = Attributes::makeRealArray
                         ("COEFNUM", "Radial profile: list of polynomial coefficients for the numerator (not for PSI-BFIELD-MIRRORED)");

    itsAttr[COEFDENOM] = Attributes::makeRealArray
                         ("COEFDENOM", "Radial profile: list of polynomial coefficients for the denominator (not for PSI-BFIELD-MIRRORED)");

    itsAttr[COEFNUMPHI]   = Attributes::makeRealArray
                         ("COEFNUMPHI", "Angular profile: list of polynomial coefficients for the numerator (not for PSI-BFIELD-MIRRORED)");

    itsAttr[COEFDENOMPHI] = Attributes::makeRealArray
                         ("COEFDENOMPHI", "Angular profile: list of polynomial coefficients for the denominator (not for PSI-BFIELD-MIRRORED)");

    itsAttr[BMAX]      = Attributes::makeReal
                         ("BMAX", "Maximum magnetic field in Tesla.");

    itsAttr[PHIMIN]    = Attributes::makeReal
                         ("PHIMIN", "Minimal azimuth [deg] (default 0)");

    itsAttr[PHIMAX]    = Attributes::makeReal
                         ("PHIMAX", "Maximal azimuth [deg] (default 360)");

    itsAttr[RMIN]      = Attributes::makeReal
                         ("RMIN", "Minimum radius [mm].");

    itsAttr[RMAX]      = Attributes::makeReal
                         ("RMAX", "Maximum radius [mm].");

    itsAttr[SLPTC]      = Attributes::makeReal
                         ("SLPTC", "Slopes of the rising edge [1/mm] (for PSI-BFIELD-MIRRORED)");


    registerOwnership(AttributeHandler::STATEMENT);

    OpalTrimCoil *defTrimCoil = clone("UNNAMED_TRIMCOIL");
    defTrimCoil->builtin = true;

    try {
        defTrimCoil->update();
        OpalData::getInstance()->define(defTrimCoil);
    } catch(...) {
        delete defTrimCoil;
    }
}


OpalTrimCoil::OpalTrimCoil(const std::string &name, OpalTrimCoil *parent):
    Definition(name, parent),
    trimcoil_m(nullptr)
{}


OpalTrimCoil::~OpalTrimCoil() {
}


bool OpalTrimCoil::canReplaceBy(Object *object) {
    // Can replace only by another trim coil.
    return dynamic_cast<OpalTrimCoil *>(object) != nullptr;
}


OpalTrimCoil *OpalTrimCoil::clone(const std::string &name) {
    return new OpalTrimCoil(name, this);
}


void OpalTrimCoil::execute() {
    update();
}


OpalTrimCoil *OpalTrimCoil::find(const std::string &name) {
    OpalTrimCoil *trimcoil = dynamic_cast<OpalTrimCoil *>(OpalData::getInstance()->find(name));

    if (trimcoil == nullptr) {
        throw OpalException("OpalTrimCoil::find()", "OpalTrimCoil \"" + name + "\" not found.");
    }
    return trimcoil;
}


void OpalTrimCoil::update() {
    // Set default name.
    if (getOpalName().empty()) setOpalName("UNNAMED_TRIMCOIL");
}


void OpalTrimCoil::initOpalTrimCoil() {
    if (trimcoil_m != nullptr) return;

    std::string type = Attributes::getString(itsAttr[TYPE]);

    double bmax   = Attributes::getReal(itsAttr[BMAX]);
    double phimin = Attributes::getReal(itsAttr[PHIMIN]);
    double phimax = Attributes::getReal(itsAttr[PHIMAX]);
    double rmin   = Attributes::getReal(itsAttr[RMIN]);
    double rmax   = Attributes::getReal(itsAttr[RMAX]);

    if (type == "PSI-BFIELD" || type == "PSI-PHASE") {
        std::vector<double> coefnum      = Attributes::getRealArray(itsAttr[COEFNUM]);
        std::vector<double> coefdenom    = Attributes::getRealArray(itsAttr[COEFDENOM]);
        std::vector<double> coefnumphi   = Attributes::getRealArray(itsAttr[COEFNUMPHI]);
        std::vector<double> coefdenomphi = Attributes::getRealArray(itsAttr[COEFDENOMPHI]);
        if (type == "PSI-BFIELD")
            trimcoil_m = std::unique_ptr<TrimCoilBFit>     (new TrimCoilBFit    (bmax, rmin, rmax, coefnum, coefdenom, coefnumphi, coefdenomphi));
        else // type == "PSI-PHASE"
            trimcoil_m = std::unique_ptr<TrimCoilPhaseFit> (new TrimCoilPhaseFit(bmax, rmin, rmax, coefnum, coefdenom, coefnumphi, coefdenomphi));

    } else if (type == "PSI-BFIELD-MIRRORED") {
        double slope = Attributes::getReal(itsAttr[SLPTC]);
        trimcoil_m = std::unique_ptr<TrimCoilMirrored>     (new TrimCoilMirrored(bmax, rmin, rmax, slope));
    }

    trimcoil_m->setAzimuth(phimin, phimax);

    *gmsg << level3 << *this << endl;
}

Inform& OpalTrimCoil::print(Inform &os) const {
    os << "* ******************************** T R I M C O I L ********************************\n"
       << "* TRIMCOIL       " << getOpalName() << '\n'
       << "* TYPE           " << Attributes::getString(itsAttr[TYPE]) << '\n';

    std::string type = Attributes::getString(itsAttr[TYPE]);
    if (type == "PSI-BFIELD" || type == "PSI-PHASE") {
        printPolynom(os,itsAttr[COEFNUM]);
        printPolynom(os,itsAttr[COEFDENOM]);
        printPolynom(os,itsAttr[COEFNUMPHI]);
        printPolynom(os,itsAttr[COEFDENOMPHI]);
    }

    os << "* BMAX           " << Attributes::getReal(itsAttr[BMAX]) << '\n'
       << "* RMIN           " << Attributes::getReal(itsAttr[RMIN]) << '\n'
       << "* RMAX           " << Attributes::getReal(itsAttr[RMAX]) << '\n';

    if (Attributes::getString(itsAttr[TYPE]) == "PSI-BFIELD-MIRRORED") {
        os << "* SLPTC          " << Attributes::getReal(itsAttr[SLPTC]) << '\n';
    }
    os << "* *********************************************************************************" << endl;
    return os;
}

void OpalTrimCoil::printPolynom(Inform& os, const Attribute& attr) const {
    std::stringstream ss;
    std::vector<double> coef = Attributes::getRealArray(attr);
    for (std::size_t i = 0; i < coef.size(); ++i) {
        ss << ((i > 0) ? "+ " : "") << coef[i]
           << ((i > 0) ? (" * x^" + std::to_string(i)) : "") << ' ';
    }
    os << "* POLYNOM " << attr.getName() << "   " << ss.str() << '\n';
}