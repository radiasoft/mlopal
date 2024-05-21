//
// Class ParticleMatterInteraction
//   The class for the OPAL PARTICLEMATTERINTERACTION command.
//
// Copyright (c) 2012 - 2021, Andreas Adelmann, Paul Scherrer Institut, Villigen PSI, Switzerland
//                            Christof Metzger-Kraus, Helmholtz-Zentrum Berlin
//                            Pedro Calvo, CIEMAT, Spain
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Optimizing the radioisotope production of the novel AMIT
// superconducting weak focusing cyclotron"
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
#ifndef OPAL_PARTICLEMATTERINTERACTION_HH
#define OPAL_PARTICLEMATTERINTERACTION_HH

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
#include "Solvers/ParticleMatterInteractionHandler.h"

class ElementBase;
class Inform;

class ParticleMatterInteraction: public Definition {

public:
    /// Exemplar constructor.
    ParticleMatterInteraction();

    virtual ~ParticleMatterInteraction();

    /// Test if replacement is allowed.
    virtual bool canReplaceBy(Object* object);

    /// Make clone.
    virtual ParticleMatterInteraction* clone(const std::string& name);

    /// Check the PARTICLEMATTERINTERACTION data.
    virtual void execute();

    /// Find named PARTICLEMATTERINTERACTION.
    static ParticleMatterInteraction* find(const std::string& name);

    /// Update the PARTICLEMATTERINTERACTION data.
    virtual void update();

    void print(std::ostream& os) const;

    void initParticleMatterInteractionHandler(ElementBase& element);

    void updateElement(ElementBase* element);

    ParticleMatterInteractionHandler* handler_m;

private:
    enum class InteractionType: unsigned short {
        SCATTERING,
        BEAMSTRIPPING
    };

    // Not implemented.
    ParticleMatterInteraction(const ParticleMatterInteraction&);
    void operator=(const ParticleMatterInteraction&);

    // Clone constructor.
    ParticleMatterInteraction(const std::string& name, ParticleMatterInteraction* parent);

    void getInteractionType();

    InteractionType type_m;
};

inline std::ostream& operator<<(std::ostream& os, const ParticleMatterInteraction& b) {
    b.print(os);
    return os;
}

#endif // OPAL_PARTICLEMATTERINTERACTION_HH
