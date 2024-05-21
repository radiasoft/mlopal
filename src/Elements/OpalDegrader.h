//
// Class OpalDegrader
//   The DEGRADER element.
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
#ifndef OPAL_OpalDegrader_HH
#define OPAL_OpalDegrader_HH

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;

class OpalDegrader: public OpalElement {

public:

    /// The attributes of class OpalDegrader.
    enum {
        XSIZE = COMMON,  // not used
        YSIZE,           // not used
        SIZE
    };

    /// Exemplar constructor.
    OpalDegrader();

    virtual ~OpalDegrader();

    /// Make clone.
    virtual OpalDegrader* clone(const std::string& name);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalDegrader(const OpalDegrader&);
    void operator=(const OpalDegrader&);

    // Clone constructor.
    OpalDegrader(const std::string& name, OpalDegrader* parent);

    ParticleMatterInteraction* parmatint_m;
};

#endif // OPAL_OpalDegrader_HH
