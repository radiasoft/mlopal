//
// Class Multipole
//   The MULTIPOLE element defines a thick multipole.
//
// Copyright (c) 2012-2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "AbsBeamline/Multipole.h"

#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"

namespace{
    enum {
        DIPOLE,
        QUADRUPOLE,
        SEXTUPOLE,
        OCTUPOLE,
        DECAPOLE
    };
}

namespace {
    unsigned int factorial(unsigned int n) {
        static int fact[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320,
                             362880, 3628800, 39916800, 479001600};
        if (n > 12) {
            throw GeneralClassicException("Multipole::factorial", "Factorial can't be larger than 12");
        }
        return fact[n];
    }
}

Multipole::Multipole():
    Multipole("")
{}


Multipole::Multipole(const Multipole &right):
    Component(right),
    NormalComponents(right.NormalComponents),
    NormalComponentErrors(right.NormalComponentErrors),
    SkewComponents(right.SkewComponents),
    SkewComponentErrors(right.SkewComponentErrors),
    max_SkewComponent_m(right.max_SkewComponent_m),
    max_NormalComponent_m(right.max_NormalComponent_m),
    nSlices_m(right.nSlices_m) {
}


Multipole::Multipole(const std::string &name):
    Component(name),
    NormalComponents(1, 0.0),
    NormalComponentErrors(1, 0.0),
    SkewComponents(1, 0.0),
    SkewComponentErrors(1, 0.0),
    max_SkewComponent_m(1),
    max_NormalComponent_m(1),
    nSlices_m(1) {
}


Multipole::~Multipole()
{}


void Multipole::accept(BeamlineVisitor &visitor) const {
    visitor.visitMultipole(*this);
}


double Multipole::getNormalComponent(int n) const {
    if (n < max_NormalComponent_m) {
        return NormalComponents[n];
    }
    return 0.0;
}


double Multipole::getSkewComponent(int n) const {
    if (n < max_SkewComponent_m) {
        return SkewComponents[n];
    }
    return 0.0;
}


void Multipole::setNormalComponent(int n, double v, double vError) {
    //   getField().setNormalComponent(n, v);
    PAssert_GE(n, 1);

    if(n >  max_NormalComponent_m) {
        max_NormalComponent_m = n;
        NormalComponents.resize(max_NormalComponent_m, 0.0);
        NormalComponentErrors.resize(max_NormalComponent_m, 0.0);
    }
    switch(n-1) {
    case DIPOLE:
        NormalComponents[n - 1] = (v + vError) / 2;
        NormalComponentErrors[n - 1] = NormalComponents[n - 1];
        break;
    default:
        NormalComponents[n - 1]      = (v + vError) / factorial(n-1);
        NormalComponentErrors[n - 1] =      vError  / factorial(n-1);
    }
}

void Multipole::setSkewComponent(int n, double v, double vError) {
    //   getField().setSkewComponent(n, v);
    PAssert_GT(n, 1);

    if(n  > max_SkewComponent_m) {
        max_SkewComponent_m = n;
        SkewComponents.resize(max_SkewComponent_m, 0.0);
        SkewComponentErrors.resize(max_SkewComponent_m, 0.0);
    }
    switch(n-1) {
    case DIPOLE:
        SkewComponents[n - 1]      = (v + vError) / 2;
        SkewComponentErrors[n - 1] = SkewComponents[n - 1];
        break;
    default:
        SkewComponents[n - 1]      = (v + vError) / factorial(n-1);
        SkewComponentErrors[n - 1] =      vError  / factorial(n-1);
    }
}

//set the number of slices for map tracking
void Multipole::setNSlices(const std::size_t& nSlices) {
    nSlices_m = nSlices;
}

//get the number of slices for map tracking
std::size_t Multipole::getNSlices() const {
    return nSlices_m;
}

void Multipole::computeField(Vector_t R, Vector_t &/*E*/, Vector_t &B) {
    {
        std::vector<Vector_t> Rn(max_NormalComponent_m + 1);
        std::vector<double> fact(max_NormalComponent_m + 1);
        Rn[0] = Vector_t(1.0);
        fact[0] = 1;
        for (int i = 0; i < max_NormalComponent_m; ++ i) {
            switch(i) {
            case DIPOLE:
                B(1) += NormalComponents[i];
                break;

            case QUADRUPOLE:
                B(0) += NormalComponents[i] * R(1);
                B(1) += NormalComponents[i] * R(0);
                break;

            case SEXTUPOLE:
                B(0) += 2 * NormalComponents[i] * R(0) * R(1);
                B(1) += NormalComponents[i] * (Rn[2](0) - Rn[2](1));
                break;

            case OCTUPOLE:
                B(0) += NormalComponents[i] * (3 * Rn[2](0) * Rn[1](1) - Rn[3](1));
                B(1) += NormalComponents[i] * (Rn[3](0) - 3 * Rn[1](0) * Rn[2](1));
                break;

            case DECAPOLE:
                B(0) += 4 * NormalComponents[i] * (Rn[3](0) * Rn[1](1) - Rn[1](0) * Rn[3](1));
                B(1) += NormalComponents[i] * (Rn[4](0) - 6 * Rn[2](0) * Rn[2](1) + Rn[4](1));
                break;

            default:
                {
                    double powMinusOne = 1;
                    double Bx = 0.0, By = 0.0;
                    for (int j = 1; j <= (i + 1) / 2; ++ j) {
                        Bx += powMinusOne * NormalComponents[i] * (Rn[i - 2 * j + 1](0) * fact[i - 2 * j + 1] *
                                                                   Rn[2 * j - 1](1) * fact[2 * j - 1]);
                        By += powMinusOne * NormalComponents[i] * (Rn[i - 2 * j + 2](0) * fact[i - 2 * j + 2] *
                                                                   Rn[2 * j - 2](1) * fact[2 * j - 2]);
                        powMinusOne *= -1;
                    }

                    if ((i + 1) / 2 == i / 2) {
                        int j = (i + 2) / 2;
                        By += powMinusOne * NormalComponents[i] * (Rn[i - 2 * j + 2](0) * fact[i - 2 * j + 2] *
                                                                   Rn[2 * j - 2](1) * fact[2 * j - 2]);
                    }
                    B(0) += Bx;
                    B(1) += By;
                }
            }

            Rn[i + 1](0) = Rn[i](0) * R(0);
            Rn[i + 1](1) = Rn[i](1) * R(1);
            fact[i + 1] = fact[i] / (i + 1);
        }
    }

    {

        std::vector<Vector_t> Rn(max_SkewComponent_m + 1);
        std::vector<double> fact(max_SkewComponent_m + 1);
        Rn[0] = Vector_t(1.0);
        fact[0] = 1;
        for (int i = 0; i < max_SkewComponent_m; ++ i) {
            switch(i) {
            case DIPOLE:
                B(0) -= SkewComponents[i];
                break;

            case QUADRUPOLE:
                B(0) -= SkewComponents[i] * R(0);
                B(1) += SkewComponents[i] * R(1);
                break;

            case SEXTUPOLE:
                B(0) -= SkewComponents[i] * (Rn[2](0) - Rn[2](1));
                B(1) += 2 * SkewComponents[i] * R(0) * R(1);
                break;

            case OCTUPOLE:
                B(0) -= SkewComponents[i] * (Rn[3](0) - 3 * Rn[1](0) * Rn[2](1));
                B(1) += SkewComponents[i] * (3 * Rn[2](0) * Rn[1](1) - Rn[3](1));
                break;

            case DECAPOLE:
                B(0) -= SkewComponents[i] * (Rn[4](0) - 6 * Rn[2](0) * Rn[2](1) + Rn[4](1));
                B(1) += 4 * SkewComponents[i] * (Rn[3](0) * Rn[1](1) - Rn[1](0) * Rn[3](1));
                break;

            default:
                {
                    double powMinusOne = 1;
                    double Bx = 0, By = 0;
                    for (int j = 1; j <= (i + 1) / 2; ++ j) {
                        Bx -= powMinusOne * SkewComponents[i] * (Rn[i - 2 * j + 2](0) * fact[i - 2 * j + 2] *
                                                                 Rn[2 * j - 2](1) * fact[2 * j - 2]);
                        By += powMinusOne * SkewComponents[i] * (Rn[i - 2 * j + 1](0) * fact[i - 2 * j + 1] *
                                                                 Rn[2 * j - 1](1) * fact[2 * j - 1]);
                        powMinusOne *= -1;
                    }

                    if ((i + 1) / 2 == i / 2) {
                        int j = (i + 2) / 2;
                        Bx -= powMinusOne * SkewComponents[i] * (Rn[i - 2 * j + 2](0) * fact[i - 2 * j + 2] *
                                                                 Rn[2 * j - 2](1) * fact[2 * j - 2]);
                    }

                    B(0) += Bx;
                    B(1) += By;
                }
            }

            Rn[i + 1](0) = Rn[i](0) * R(0);
            Rn[i + 1](1) = Rn[i](1) * R(1);
            fact[i + 1] = fact[i] / (i + 1);
        }
    }
}


bool Multipole::apply(const size_t &i, const double &, Vector_t &E, Vector_t &B) {
    const Vector_t &R = RefPartBunch_m->R[i];
    if(R(2) < 0.0 || R(2) > getElementLength()) return false;
    if (!isInsideTransverse(R)) return getFlagDeleteOnTransverseExit();

    Vector_t Ef(0.0), Bf(0.0);
    computeField(R, Ef, Bf);

    for (unsigned int d = 0; d < 3; ++ d) {
        E[d] += Ef(d);
        B[d] += Bf(d);
    }

    return false;
}

bool Multipole::apply(const Vector_t &R, const Vector_t &, const double &, Vector_t &E, Vector_t &B) {
    if(R(2) < 0.0 || R(2) > getElementLength()) return false;
    if (!isInsideTransverse(R)) return getFlagDeleteOnTransverseExit();

    computeField(R, E, B);

    return false;
}

bool Multipole::applyToReferenceParticle(const Vector_t &R, const Vector_t &, const double &, Vector_t &E, Vector_t &B) {
    if(R(2) < 0.0 || R(2) > getElementLength()) return false;
    if (!isInsideTransverse(R)) return true;

    for (int i = 0; i < max_NormalComponent_m; ++ i) {
        NormalComponents[i] -= NormalComponentErrors[i];
    }
    for (int i = 0; i < max_SkewComponent_m; ++ i) {
        SkewComponents[i] -= SkewComponentErrors[i];
    }

    computeField(R, E, B);

    for (int i = 0; i < max_NormalComponent_m; ++ i) {
        NormalComponents[i] += NormalComponentErrors[i];
    }
    for (int i = 0; i < max_SkewComponent_m; ++ i) {
        SkewComponents[i] += SkewComponentErrors[i];
    }

    return false;
}

void Multipole::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    RefPartBunch_m = bunch;
    endField = startField + getElementLength();
    online_m = true;
}


void Multipole::finalise() {
    online_m = false;
}

bool Multipole::bends() const {
    return false;
}


void Multipole::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = 0.0;
    zEnd = getElementLength();
}


ElementType Multipole::getType() const {
    return ElementType::MULTIPOLE;
}


bool Multipole::isInside(const Vector_t &r) const {
    if (r(2) >= 0.0 && r(2) < getElementLength()) {
        return isInsideTransverse(r);
    }

    return false;
}

bool Multipole::isFocusing(unsigned int component) const {
    if (component >= NormalComponents.size()) throw GeneralClassicException("Multipole::isFocusing", "component too big");

    return NormalComponents[component] * std::pow(-1, component + 1) * RefPartBunch_m->getChargePerParticle() > 0.0;
}
