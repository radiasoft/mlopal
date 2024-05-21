//
// Class OpalVerticalFFAMagnet
//   The class provides the user interface for the VERTICALFFA object
//
// Copyright (c) 2019 Chris Rogers
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
#include "Attributes/Attributes.h"

#include "AbsBeamline/EndFieldModel/Tanh.h" // classic
#include "AbsBeamline/VerticalFFAMagnet.h" // classic
#include "Elements/OpalVerticalFFAMagnet.h"
#include "Physics/Units.h"

std::string OpalVerticalFFAMagnet::docstring_m =
      std::string("The \"VerticalFFAMagnet\" element defines a vertical FFA ")+
      std::string("magnet, which has a field that increases in the vertical ")+
      std::string("direction while maintaining similar orbits.");
OpalVerticalFFAMagnet::OpalVerticalFFAMagnet() :
    OpalElement(SIZE, "VERTICALFFAMAGNET", docstring_m.c_str()) {
    itsAttr[B0] = Attributes::makeReal
          ("B0", "The nominal dipole field of the magnet at zero height [T].");
    itsAttr[FIELD_INDEX] = Attributes::makeReal("FIELD_INDEX",
      "Exponent term in the field index [m^(-1)].");
    itsAttr[WIDTH] = Attributes::makeReal("WIDTH",
      "The full width of the magnet. Particles moving more than WIDTH/2 horizontally, in either direction, are out of the aperture.");
    itsAttr[MAX_HORIZONTAL_POWER] = Attributes::makeReal("MAX_HORIZONTAL_POWER",
      "The maximum power in horizontal coordinate that will be considered in the field expansion.");
    itsAttr[END_LENGTH] = Attributes::makeReal("END_LENGTH",
      "The end length of the FFA fringe field [m].");
    itsAttr[CENTRE_LENGTH] = Attributes::makeReal("CENTRE_LENGTH",
      "The centre length of the FFA (i.e. length of the flat top) [m].");
    itsAttr[BB_LENGTH] = Attributes::makeReal("BB_LENGTH",
      "Determines the length of the bounding box. Magnet is situated symmetrically in the bounding box. [m]");
    itsAttr[HEIGHT_POS_EXTENT] = Attributes::makeReal("HEIGHT_POS_EXTENT",
      "Height of the magnet above z=0. Particles moving upwards more than HEIGHT_POS_EXTENT are out of the aperture [m].");
    itsAttr[HEIGHT_NEG_EXTENT] = Attributes::makeReal("HEIGHT_NEG_EXTENT",
      "Height of the magnet below z=0. Particles moving downwards more than HEIGHT_NEG_EXTENT are out of the aperture [m].");

    registerOwnership();

    VerticalFFAMagnet* magnet = new VerticalFFAMagnet("VerticalFFAMagnet");
    magnet->setEndField(new endfieldmodel::Tanh(1., 1., 1));
    setElement(magnet);
}


OpalVerticalFFAMagnet::OpalVerticalFFAMagnet(const std::string &name,
                                             OpalVerticalFFAMagnet *parent) :
    OpalElement(name, parent) {
    VerticalFFAMagnet* magnet = new VerticalFFAMagnet(name);
    magnet->setEndField(new endfieldmodel::Tanh(1., 1., 1));
    setElement(magnet);
}


OpalVerticalFFAMagnet::~OpalVerticalFFAMagnet() {
}


OpalVerticalFFAMagnet *OpalVerticalFFAMagnet::clone(const std::string &name) {
    return new OpalVerticalFFAMagnet(name, this);
}


void OpalVerticalFFAMagnet::update() {
    VerticalFFAMagnet *magnet =
              dynamic_cast<VerticalFFAMagnet*>(getElement());
    magnet->setB0(Attributes::getReal(itsAttr[B0]));
    int maxOrder = floor(Attributes::getReal(itsAttr[MAX_HORIZONTAL_POWER]));
    magnet->setMaxOrder(maxOrder);
    magnet->setFieldIndex(Attributes::getReal(itsAttr[FIELD_INDEX]));
    magnet->setNegativeVerticalExtent(Attributes::getReal(itsAttr[HEIGHT_NEG_EXTENT]));
    magnet->setPositiveVerticalExtent(Attributes::getReal(itsAttr[HEIGHT_POS_EXTENT]));
    magnet->setBBLength(Attributes::getReal(itsAttr[BB_LENGTH]));
    magnet->setWidth(Attributes::getReal(itsAttr[WIDTH]));

    // get centre length and end length in radians
    endfieldmodel::Tanh* endField = dynamic_cast<endfieldmodel::Tanh*>(magnet->getEndField());
    double end_length = Attributes::getReal(itsAttr[END_LENGTH]) * Units::m2mm;
    double centre_length = Attributes::getReal(itsAttr[CENTRE_LENGTH]) * Units::m2mm;
    endField->setLambda(end_length);
    // x0 is the distance between B=0.5*B0 and B=B0 i.e. half the centre length
    endField->setX0(centre_length/2.);
    endField->setTanhDiffIndices(maxOrder+2);
    magnet->initialise();
    setElement(magnet);
}