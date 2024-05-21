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
#ifndef CLASSIC_CCollimator_HH
#define CLASSIC_CCollimator_HH

#include "AbsBeamline/PluginElement.h"

class ParticleMatterInteractionHandler;

class CCollimator: public PluginElement {

public:
    /// Constructor with given name.
    explicit CCollimator(const std::string& name);

    CCollimator();
    CCollimator(const CCollimator& rhs);
    void operator=(const CCollimator&) = delete;
    virtual ~CCollimator();

    /// Apply visitor to CCollimator.
    virtual void accept(BeamlineVisitor&) const override;
    ///@{ Override implementation of PluginElement
    virtual void goOnline(const double& kineticEnergy) override;
    virtual ElementType getType() const override;
    virtual void getDimensions(double& zBegin, double& zEnd) const override;
    ///@}

    /// unused check method
    // bool checkCollimator(Vector_t r, Vector_t rmin, Vector_t rmax);

    /// Some debug print
    void print();

    /// Set dimensions and consistency checks
    void setDimensions(double xstart, double xend,
                       double ystart, double yend,
                       double zstart, double zend,
                       double width);
    /// unhide PluginElement::setDimensions(double xstart, double xend, double ystart, double yend)
    using PluginElement::setDimensions;

    ///@{ Member variable access
    double getZStart() ;
    double getZEnd() ;
    double getWidth() ;
    ///@}
private:
    /// Initialise particle matter interaction
    virtual void doInitialise(PartBunchBase<double, 3>* bunch) override;
    /// Record hits when bunch particles pass
    virtual bool doCheck(PartBunchBase<double, 3>* bunch, const int turnnumber, const double t, const double tstep) override;
    /// Calculate extend in r
    virtual void doSetGeom() override;
    /// Virtual hook for finalise
    virtual void doFinalise() override;
    /// Virtual hook for preCheck
    virtual bool doPreCheck(PartBunchBase<double, 3>*) override;
    /// Virtual hook for finaliseCheck
    virtual bool doFinaliseCheck(PartBunchBase<double, 3>* bunch, bool flagNeedUpdate) override;

    bool informed_m = false; ///< Flag if error information already printed

    ///@{ input geometry positions
    double zstart_m;
    double zend_m;
    double width_m;
    ///@}
    double rmax_m; ///< maximum extend in r

    ParticleMatterInteractionHandler* parmatint_m = nullptr;
};

inline
double CCollimator::getZStart() {
    return zstart_m;
}

inline
double CCollimator::getZEnd() {
    return zend_m;
}

inline
double CCollimator::getWidth() {
    return width_m;
}
#endif // CLASSIC_CCollimator_HH
