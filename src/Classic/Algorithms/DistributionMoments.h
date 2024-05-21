//
// Class DistributionMoments
//   Computes the statistics of particle distributions.
//
// Copyright (c) 2021, Christof Metzger-Kraus
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
#ifndef DISTRIBUTIONMOMENTS_H
#define DISTRIBUTIONMOMENTS_H

#include "FixedAlgebra/FMatrix.h"

#include "Vektor.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"

#include <vector>

class OpalParticle;
template<class T, unsigned Dim>
class PartBunchBase;

class DistributionMoments {
public:
    DistributionMoments();

    void compute(const std::vector<OpalParticle>::const_iterator &,
                 const std::vector<OpalParticle>::const_iterator &);
    void compute(PartBunchBase<double, 3> const&);
    void computeMeanKineticEnergy(PartBunchBase<double, 3> const&);
    void computeDebyeLength(PartBunchBase<double, 3> const&, double);
    void computePlasmaParameter(double);

    Vector_t getMeanPosition() const;
    Vector_t getStandardDeviationPosition() const;
    Vector_t getMeanMomentum() const;
    Vector_t getStandardDeviationMomentum() const;
    Vector_t getNormalizedEmittance() const;
    Vector_t getGeometricEmittance() const;
    Vector_t getStandardDeviationRP() const;
    Vector_t getHalo() const;
    Vector_t getMaxR() const;

    Vector_t get68Percentile() const;
    Vector_t get95Percentile() const;
    Vector_t get99Percentile() const;
    Vector_t get99_99Percentile() const;
    Vector_t getNormalizedEmittance68Percentile() const;
    Vector_t getNormalizedEmittance95Percentile() const;
    Vector_t getNormalizedEmittance99Percentile() const;
    Vector_t getNormalizedEmittance99_99Percentile() const;

    double getMeanTime() const;
    double getStdTime() const;
    double getMeanGamma() const;
    double getMeanKineticEnergy() const;
    double getTemperature() const;
    double getDebyeLength() const;
    double getPlasmaParameter() const;
    double getStdKineticEnergy() const;
    double getDx() const;
    double getDDx() const;
    double getDy() const;
    double getDDy() const;
    FMatrix<double, 6, 6> getMoments6x6() const;
    double getTotalCharge() const;
    double getTotalMass() const;
    double getTotalNumParticles() const;

private:
    bool isParticleExcluded(const OpalParticle &) const;
    template<class InputIt>
    void computeMeans(const InputIt &, const InputIt &);
    template<class InputIt>
    void computeStatistics(const InputIt &, const InputIt &);
    template<class InputIt>
    void computePercentiles(const InputIt &, const InputIt &);
    using iterator_t = std::vector<Vektor<double, 2>>::const_iterator;
    std::pair<double, iterator_t> determinePercentilesDetail(const iterator_t& begin, const iterator_t& end,
                                                             const std::vector<int>& globalAccumulatedHistogram,
                                                             const std::vector<int>& localAccumulatedHistogram,
                                                             unsigned int dimension,
                                                             int numRequiredParticles) const;
    double computeNormalizedEmittance(const iterator_t& begin, const iterator_t& end) const;
    void fillMembers(std::vector<double> const&);
    void reset();
    void resetPlasmaParameters();

    Vector_t meanR_m;
    Vector_t meanP_m;
    Vector_t stdR_m;
    Vector_t stdP_m;
    Vector_t stdRP_m;
    Vector_t normalizedEps_m;
    Vector_t geometricEps_m;
    Vector_t halo_m;
    Vector_t maxR_m;
    Vector_t minR_m;
    Vector_t sixtyEightPercentile_m;
    Vector_t normalizedEps68Percentile_m;
    Vector_t ninetyFivePercentile_m;
    Vector_t normalizedEps95Percentile_m;
    Vector_t ninetyNinePercentile_m;
    Vector_t normalizedEps99Percentile_m;
    Vector_t ninetyNine_NinetyNinePercentile_m;
    Vector_t normalizedEps99_99Percentile_m;

    double meanTime_m;
    double stdTime_m;
    double meanKineticEnergy_m;
    double temperature_m;
    double debyeLength_m;
    double plasmaParameter_m;
    double stdKineticEnergy_m;
    double meanGamma_m;
    double centroid_m[6];
    FMatrix<double, 6, 6> moments_m;

    double totalCharge_m;
    double totalMass_m;
    unsigned int totalNumParticles_m;

    static const double percentileOneSigmaNormalDist_m;
    static const double percentileTwoSigmasNormalDist_m;
    static const double percentileThreeSigmasNormalDist_m;
    static const double percentileFourSigmasNormalDist_m;
};

inline
Vector_t DistributionMoments::getMeanPosition() const
{
    return meanR_m;
}

inline
Vector_t DistributionMoments::getStandardDeviationPosition() const
{
    return stdR_m;
}

inline
Vector_t DistributionMoments::getMeanMomentum() const
{
    return meanP_m;
}

inline
Vector_t DistributionMoments::getStandardDeviationMomentum() const
{
    return stdP_m;
}

inline
Vector_t DistributionMoments::getNormalizedEmittance() const
{
    return normalizedEps_m;
}

inline
Vector_t DistributionMoments::getGeometricEmittance() const
{
    return geometricEps_m;
}

inline
Vector_t DistributionMoments::getStandardDeviationRP() const
{
    return stdRP_m;
}

inline
Vector_t DistributionMoments::getHalo() const
{
    return halo_m;
}

inline
double DistributionMoments::getMeanTime() const
{
    return meanTime_m;
}

inline
double DistributionMoments::getStdTime() const
{
    return stdTime_m;
}

inline
double DistributionMoments::getMeanGamma() const
{
    return meanGamma_m;
}

inline
double DistributionMoments::getMeanKineticEnergy() const
{
    return meanKineticEnergy_m;
}

// Compute and return the value of temperature in K
inline
double DistributionMoments::getTemperature() const
{
    return (temperature_m / 
           (Physics::kB * Units::eV2kg * 
            Physics::c * Physics::c));
}
inline
double DistributionMoments::getDebyeLength() const
{
    return debyeLength_m;
}
inline
double DistributionMoments::getPlasmaParameter() const
{
    return plasmaParameter_m;
}

inline
double DistributionMoments::getStdKineticEnergy() const
{
    return stdKineticEnergy_m;
}

inline
double DistributionMoments::getDx() const
{
    return moments_m(0, 5);
}

inline
double DistributionMoments::getDDx() const
{
    return moments_m(1, 5);
}

inline
double DistributionMoments::getDy() const
{
    return moments_m(2, 5);
}

inline
double DistributionMoments::getDDy() const
{
    return moments_m(3, 5);
}

inline 
FMatrix<double, 6, 6>  DistributionMoments::getMoments6x6() const
{
    return moments_m;
}

inline
double DistributionMoments::getTotalCharge() const
{
    return totalCharge_m;
}

inline
double DistributionMoments::getTotalMass() const
{
    return totalMass_m;
}

inline
double DistributionMoments::getTotalNumParticles() const
{
    return totalNumParticles_m;
}

inline
Vector_t DistributionMoments::get68Percentile() const
{
    return sixtyEightPercentile_m;
}

inline
Vector_t DistributionMoments::getNormalizedEmittance68Percentile() const
{
    return normalizedEps68Percentile_m;
}

inline
Vector_t DistributionMoments::get95Percentile() const
{
    return ninetyFivePercentile_m;
}

inline
Vector_t DistributionMoments::getNormalizedEmittance95Percentile() const
{
    return normalizedEps95Percentile_m;
}

inline
Vector_t DistributionMoments::get99Percentile() const
{
    return ninetyNinePercentile_m;
}

inline
Vector_t DistributionMoments::getNormalizedEmittance99Percentile() const
{
    return normalizedEps99Percentile_m;
}

inline
Vector_t DistributionMoments::get99_99Percentile() const
{
    return ninetyNine_NinetyNinePercentile_m;
}

inline
Vector_t DistributionMoments::getNormalizedEmittance99_99Percentile() const
{
    return normalizedEps99_99Percentile_m;
}

inline
Vector_t DistributionMoments::getMaxR() const
{
    Vector_t maxDistance;
    for (unsigned int i = 0; i < 3; ++ i) {
        maxDistance[i] = std::max(std::abs(maxR_m[i]), std::abs(minR_m[i]));
    }
    return maxDistance;
}

#endif
