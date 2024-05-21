//
// Class OpalCyclotron
//   The OpalCyclotron element.
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
#include "Elements/OpalCyclotron.h"

#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/CyclotronRep.h"
#include "Physics/Units.h"
#include "Structure/BoundaryGeometry.h"
#include "TrimCoils/OpalTrimCoil.h"
#include "TrimCoils/TrimCoil.h"
#include "Utilities/OpalException.h"

#include <numeric>


OpalCyclotron::OpalCyclotron():
    OpalElement(SIZE, "CYCLOTRON",
                "The \"CYCLOTRON\" defines an cyclotron"),
    obgeo_m(nullptr)  {
    itsAttr[CYHARMON] = Attributes::makeReal
        ("CYHARMON", "the harmonic number of the cyclotron");

    itsAttr[SYMMETRY] = Attributes::makeReal
        ("SYMMETRY", "defines how the field is stored");

    itsAttr[RINIT] = Attributes::makeReal
        ("RINIT", "Initial radius of the reference particle [mm]");

    itsAttr[PRINIT] = Attributes::makeReal
        ("PRINIT", "Initial radial momentum of the reference particle, pr = beta_r * gamma");

    itsAttr[PHIINIT] = Attributes::makeReal
        ("PHIINIT", "Initial azimuth of the reference particle [deg]");

    itsAttr[ZINIT] = Attributes::makeReal
        ("ZINIT", "Initial z-coordinate of the reference particle [mm]. Default = 0 mm", 0.0);

    itsAttr[PZINIT] = Attributes::makeReal
        ("PZINIT", "Initial vertical momentum of the reference particle pz = beta_z * gamma. Default = 0", 0.0);

    itsAttr[FMAPFN] = Attributes::makeString
        ("FMAPFN", "Filename for the B fieldmap");

    itsAttr[BSCALE] = Attributes::makeReal
         ("BSCALE", "Scale factor for the B-field", 1.0);

    itsAttr[RFFREQ] = Attributes::makeRealArray
        ("RFFREQ", "RF Frequency(ies) [MHz]");

    itsAttr[ESCALE] = Attributes::makeRealArray
        ("ESCALE", "Scale factor for the RF field(s)");

    itsAttr[SUPERPOSE] = Attributes::makeBoolArray
        ("SUPERPOSE", "If TRUE, all of the electric field maps are superposed, only used when TYPE = BANDRF");

    itsAttr[RFMAPFN] = Attributes::makeStringArray
        ("RFMAPFN", "Filename(s) for the RF fieldmap(s)");

    itsAttr[RFFCFN] = Attributes::makeStringArray
        ("RFFCFN", "Filename(s) for the RF Frequency Coefficients");

    itsAttr[RFVCFN] = Attributes::makeStringArray
        ("RFVCFN", "Filename(s) for the RF Voltage Coefficients");

    itsAttr[RFPHI] = Attributes::makeRealArray
        ("RFPHI", "Initial phase(s) of the electric field map(s) [deg]");

    itsAttr[MINZ] = Attributes::makeReal
        ("MINZ","Minimal vertical extent of the machine [mm]",-10000.0);

    itsAttr[MAXZ] = Attributes::makeReal
        ("MAXZ","Maximal vertical extent of the machine [mm]",10000.0);

    itsAttr[MINR] = Attributes::makeReal
        ("MINR","Minimal radial extent of the machine [mm]", 0.0);

    itsAttr[MAXR] = Attributes::makeReal
        ("MAXR","Maximal radial extent of the machine [mm]", 10000.0);

    itsAttr[GEOMETRY] = Attributes::makeString
        ("GEOMETRY", "Boundary Geometry for the Cyclotron");

    itsAttr[FMLOWE] = Attributes::makeReal
        ("FMLOWE", "Minimal Energy [GeV] the fieldmap can accept. Used in GAUSSMATCHED", -1.0);

    itsAttr[FMHIGHE] = Attributes::makeReal
        ("FMHIGHE","Maximal Energy [GeV] the fieldmap can accept. Used in GAUSSMATCHED", -1.0);

    itsAttr[SPIRAL] = Attributes::makeBool
        ("SPIRAL","Flag whether or not this is a spiral inflector simulation", false);

    itsAttr[TRIMCOILTHRESHOLD] = Attributes::makeReal
        ("TRIMCOILTHRESHOLD","Minimum magnetic field [T] for which trim coils are applied", 0.0);

    itsAttr[TRIMCOIL] = Attributes::makeStringArray
        ("TRIMCOIL", "List of trim coils");

    registerOwnership();

    setElement(new CyclotronRep("CYCLOTRON"));
}

OpalCyclotron::OpalCyclotron(const std::string& name, OpalCyclotron* parent):
    OpalElement(name, parent),
    obgeo_m(nullptr) {
    setElement(new CyclotronRep(name));
}


OpalCyclotron::~OpalCyclotron()
{}


OpalCyclotron* OpalCyclotron::clone(const std::string& name) {
    return new OpalCyclotron(name, this);
}


void OpalCyclotron::update() {

    CyclotronRep* cycl =
        dynamic_cast<CyclotronRep*>(getElement());

    std::string fmapfm = Attributes::getString(itsAttr[FMAPFN]);
    std::string type   = Attributes::getString(itsAttr[TYPE]);

    double harmnum  = Attributes::getReal(itsAttr[CYHARMON]);
    double symmetry = Attributes::getReal(itsAttr[SYMMETRY]);
    double rinit    = Attributes::getReal(itsAttr[RINIT]);
    double prinit   = Attributes::getReal(itsAttr[PRINIT]);
    double phiinit  = Attributes::getReal(itsAttr[PHIINIT]);
    double zinit    = Attributes::getReal(itsAttr[ZINIT]);
    double pzinit   = Attributes::getReal(itsAttr[PZINIT]);
    double bscale   = Attributes::getReal(itsAttr[BSCALE]);

    double minz = Attributes::getReal(itsAttr[MINZ]);
    double maxz = Attributes::getReal(itsAttr[MAXZ]);
    double minr = Attributes::getReal(itsAttr[MINR]);
    double maxr = Attributes::getReal(itsAttr[MAXR]);

    double fmLowE  = Attributes::getReal(itsAttr[FMLOWE]);
    double fmHighE = Attributes::getReal(itsAttr[FMHIGHE]);

    bool spiral_flag = Attributes::getBool(itsAttr[SPIRAL]);
    double trimCoilThreshold = Attributes::getReal(itsAttr[TRIMCOILTHRESHOLD]);

    cycl->setFieldMapFN(fmapfm);
    cycl->setSymmetry(symmetry);

    cycl->setRinit(rinit);
    cycl->setPRinit(prinit);
    cycl->setPHIinit(phiinit);
    cycl->setZinit(zinit);
    cycl->setPZinit(pzinit);

    cycl->setBScale(bscale);

    cycl->setCyclotronType(type);
    cycl->setCyclHarm(harmnum);

    cycl->setMinR(minr);
    cycl->setMaxR(maxr);
    cycl->setMinZ(minz);
    cycl->setMaxZ(maxz);

    cycl->setFMLowE(fmLowE * Units::GeV2MeV);
    cycl->setFMHighE(fmHighE * Units::GeV2MeV);

    cycl->setSpiralFlag(spiral_flag);
    cycl->setTrimCoilThreshold(trimCoilThreshold);

    cycl->setOutputFN(Attributes::getString(itsAttr[OUTFN]));

    std::vector<std::string> fm_str     = Attributes::getStringArray(itsAttr[RFMAPFN]);
    std::vector<std::string> rffcfn_str = Attributes::getStringArray(itsAttr[RFFCFN]);
    std::vector<std::string> rfvcfn_str = Attributes::getStringArray(itsAttr[RFVCFN]);
    std::vector<double> scale_str       = Attributes::getRealArray(itsAttr[ESCALE]);
    std::vector<double> phi_str         = Attributes::getRealArray(itsAttr[RFPHI]);
    std::vector<double> rff_str         = Attributes::getRealArray(itsAttr[RFFREQ]);
    std::vector<bool> superpose_str     = Attributes::getBoolArray(itsAttr[SUPERPOSE]);
    std::vector<std::string> trimcoil   = Attributes::getStringArray(itsAttr[TRIMCOIL]);

    if ( !trimcoil.empty() ) {

        std::vector<TrimCoil* > trimcoils;

        for (std::vector<std::string>::const_iterator tit = trimcoil.begin();
             tit != trimcoil.end(); ++tit)
        {
            OpalTrimCoil *tc = OpalTrimCoil::find(*tit);

            if ( tc ) {
                tc->initOpalTrimCoil();
                trimcoils.push_back(tc->trimcoil_m.get());
            }
        }
        cycl->setTrimCoils(trimcoils);
    }

    cycl->setRfPhi(phi_str);
    cycl->setEScale(scale_str);
    cycl->setRfFieldMapFN(fm_str);
    cycl->setRFFCoeffFN(rffcfn_str);
    cycl->setRFVCoeffFN(rfvcfn_str);
    cycl->setRfFrequ(rff_str);
    cycl->setSuperpose(superpose_str);

    if (itsAttr[GEOMETRY] && obgeo_m == nullptr) {
      obgeo_m = BoundaryGeometry::find(Attributes::getString(itsAttr[GEOMETRY]));
      if (obgeo_m) {
          cycl->setBoundaryGeometry(obgeo_m);
      }
    }

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(cycl);
}