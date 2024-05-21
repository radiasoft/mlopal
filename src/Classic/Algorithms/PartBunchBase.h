//
// Class PartBunchBase
//   Base class for representing particle bunches.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef PART_BUNCH_BASE_H
#define PART_BUNCH_BASE_H

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/DistributionMoments.h"
#include "Algorithms/OpalParticle.h"
#include "Algorithms/PBunchDefs.h"
#include "Algorithms/Quaternion.h"
#include "Algorithms/Vektor.h"
#include "Distribution/Distribution.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include "Particle/AbstractParticle.h"
#include "Particle/ParticleAttrib.h"
#include "Physics/ParticleProperties.h"
#include "Physics/Units.h"
#include "Structure/FieldSolver.h"
#include "Utilities/GeneralClassicException.h"
#include "Utility/IpplTimings.h"

#include <memory>
#include <utility>
#include <vector>

class Distribution;
class FieldSolver;
class PartBins;
class PartBinsCyc;
class PartData;

template <class T, unsigned Dim>
class PartBunchBase : std::enable_shared_from_this<PartBunchBase<T, Dim>>
{
public:
    typedef typename AbstractParticle<T, Dim>::ParticlePos_t   ParticlePos_t;
    typedef typename AbstractParticle<T, Dim>::ParticleIndex_t ParticleIndex_t;
    typedef typename AbstractParticle<T, Dim>::UpdateFlags     UpdateFlags_t;
    typedef typename AbstractParticle<T, Dim>::Position_t      Position_t;

    typedef std::pair<Vector_t, Vector_t> VectorPair_t;

    static const unsigned Dimension = Dim;

    enum UnitState_t { units = 0, unitless = 1 };

public:
    virtual ~PartBunchBase() { }

    PartBunchBase(AbstractParticle<T, Dim>* pb, const PartData* ref);

    PartBunchBase(const PartBunchBase& rhs) = delete; // implement if needed

    /*
     * Bunch common member functions
     */

    // This is required since we initialize the Layout and the RegionLayout with default constructor
    virtual void initialize(FieldLayout_t* fLayout) = 0;

    bool getIfBeamEmitting();

    int getLastEmittedEnergyBin();

    size_t getNumberOfEmissionSteps();

    int getNumberOfEnergyBins();

    void Rebin();

    void setEnergyBins(int numberOfEnergyBins);

    bool weHaveEnergyBins();

    //FIXME: unify methods, use convention that all particles have own dt
    void switchToUnitlessPositions(bool use_dt_per_particle = false);

    //FIXME: unify methods, use convention that all particles have own dt
    void switchOffUnitlessPositions(bool use_dt_per_particle = false);

    void setDistribution(Distribution* d,
                         std::vector<Distribution*> addedDistributions,
                         size_t& np);
    void setDistribution(Distribution* d,
                         size_t numberOfParticles,
                         double current, const Beamline& bl);

    bool isGridFixed() const;

    bool hasBinning() const;


    /*
       Energy bins related functions
     */

    void setTEmission(double t);

    double getTEmission();

    bool doEmission();

    bool weHaveBins() const;

    void setPBins(PartBins* pbin);

    void setPBins(PartBinsCyc* pbin);

    /** \brief Emit particles in the given bin
        i.e. copy the particles from the bin structure into the
        particle container
    */
    size_t emitParticles(double eZ);

    void updateNumTotal();

    void rebin();

    int getLastemittedBin();

    void setLocalBinCount(size_t num, int bin);

    /** \brief Compute the gammas of all bins */
    void calcGammas();

    void calcGammas_cycl();

    /** \brief Compute the (global) Debye length for the beam */
    void calcDebyeLength();

    /** \brief Get gamma of one bin */
    double getBinGamma(int bin);

    /** \brief Set the charge of one bin to the value of q and all other to zero */
    virtual void setBinCharge(int bin, double q);

    /** \brief Set the charge of all other the ones in bin to zero */
    virtual void setBinCharge(int bin);

    /** \brief returns the number of particles outside of a box defined by x */
    size_t calcNumPartsOutside(Vector_t x);

    void calcLineDensity(unsigned int nBins, std::vector<double>& lineDensity,
                         std::pair<double, double>& meshInfo);

    void setBeamFrequency(double v);

    /*
       Mesh and Field Layout related functions
     */

    virtual void boundp();

    /** delete particles which are too far away from the center of beam*/
    void boundp_destroyCycl();

    /** This is only temporary in order to get the collimator and pepperpot working */
    size_t boundp_destroyT();

    size_t destroyT();

    /*
       Read out coordinates
     */
    virtual double getPx(int i);
    virtual double getPy(int i);
    virtual double getPz(int i);

    virtual double getPx0(int i);
    virtual double getPy0(int i);

    virtual double getX(int i);
    virtual double getY(int i);
    virtual double getZ(int i);

    virtual double getX0(int i);
    virtual double getY0(int i);

    virtual void setZ(int i, double zcoo);

    void get_bounds(Vector_t& rmin, Vector_t& rmax) const;

    void getLocalBounds(Vector_t& rmin, Vector_t& rmax) const;

    std::pair<Vector_t, double> getBoundingSphere();

    std::pair<Vector_t, double> getLocalBoundingSphere();


    /*
       Compatibility function push_back
     */
    void push_back(OpalParticle const& p);

    void setParticle(FVector<double, 6> z, int ii);

    void setParticle(OpalParticle const& p, int ii);

    OpalParticle getParticle(int ii);

    class ConstIterator {
        friend class PartBunchBase<T, Dim>;

    public:
        ConstIterator():
            bunch_m(nullptr),
            index_m(0)
        {}
        ConstIterator(PartBunchBase const* bunch, unsigned int i):
            bunch_m(bunch),
            index_m(i)
        {}

        ~ConstIterator()
        {}

        bool operator == (ConstIterator const& rhs) const
        {
            return bunch_m == rhs.bunch_m && index_m == rhs.index_m;
        }

        bool operator != (ConstIterator const& rhs) const
        {
            return bunch_m != rhs.bunch_m || index_m != rhs.index_m;
        }

        OpalParticle operator*() const
        {
            if (index_m >= bunch_m->getLocalNum()) {
                throw GeneralClassicException("PartBunchBase::ConstIterator::operator*", "out of bounds");
            }
            return OpalParticle(bunch_m->ID[index_m],
                                bunch_m->R[index_m],
                                bunch_m->P[index_m],
                                bunch_m->getT(),
                                bunch_m->Q[index_m],
                                bunch_m->getM() * Units::eV2MeV);
        }

        ConstIterator operator++()
        {
            ++index_m;
            return *this;
        }

        ConstIterator operator++(int)
        {
            ConstIterator it = *this;
            ++index_m;

            return it;
        }

        int operator-(const ConstIterator& other) const
        {
            return index_m - other.index_m;
        }
    private:
        PartBunchBase const* bunch_m;
        unsigned int index_m;
    };

    ConstIterator begin() const {
        return ConstIterator(this, 0);
    }

    ConstIterator end() const {
        return ConstIterator(this, getLocalNum());
    }

    /// Return maximum amplitudes.
    //  The matrix [b]D[/b] is used to normalise the first two modes.
    //  The maximum normalised amplitudes for these modes are stored
    //  in [b]axmax[/b] and [b]aymax[/b].
    void maximumAmplitudes(const FMatrix<double, 6, 6>& D,
                           double& axmax, double& aymax);

    void   setdT(double dt);
    double getdT() const;

    void   setT(double t);
    void   incrementT();
    double getT() const;

    /**
     * get the spos of the primary beam
     *
     * @param none
     *
     */
    double get_sPos() const;

    void set_sPos(double s);

    double get_gamma() const;
    double get_meanKineticEnergy() const;
    double get_temperature() const;
    double get_debyeLength() const;
    double get_plasmaParameter() const;
    double get_rmsDensity() const;
    Vector_t get_origin() const;
    Vector_t get_maxExtent() const;
    Vector_t get_centroid() const;
    Vector_t get_rrms() const;
    Vector_t get_rprms() const;
    Vector_t get_rmean() const;
    Vector_t get_prms() const;
    Vector_t get_pmean() const;
    Vector_t get_pmean_Distribution() const;
    Vector_t get_emit() const;
    Vector_t get_norm_emit() const;
    Vector_t get_halo() const;
    Vector_t get_68Percentile() const;
    Vector_t get_95Percentile() const;
    Vector_t get_99Percentile() const;
    Vector_t get_99_99Percentile() const;
    Vector_t get_normalizedEps_68Percentile() const;
    Vector_t get_normalizedEps_95Percentile() const;
    Vector_t get_normalizedEps_99Percentile() const;
    Vector_t get_normalizedEps_99_99Percentile() const;
    virtual Vector_t get_hr() const;

    double get_Dx() const;
    double get_Dy() const;
    double get_DDx() const;
    double get_DDy() const;

    virtual void set_meshEnlargement(double dh);

    void gatherLoadBalanceStatistics();
    size_t getLoadBalance(int p) const;

    void get_PBounds(Vector_t &min, Vector_t &max) const;

    void calcBeamParameters();
    void calcBeamParametersInitial(); // Calculate initial beam parameters before emission.

    double getCouplingConstant() const;
    void setCouplingConstant(double c);

    // set the charge per simulation particle
    void setCharge(double q);
    // set the charge per simulation particle when total particle number equals 0
    void setChargeZeroPart(double q);

    // set the mass per simulation particle
    void setMass(double mass);
    void setMassZeroPart(double mass);

    /// get the total charge per simulation particle
    double getCharge() const;

    /// get the macro particle charge
    double getChargePerParticle() const;

    double getMassPerParticle() const;

    virtual void setSolver(FieldSolver *fs);

    bool hasFieldSolver();

    FieldSolverType getFieldSolverType() const;

    void setStepsPerTurn(int n);
    int getStepsPerTurn() const;

    /// step in multiple TRACK commands
    void setGlobalTrackStep(long long n);
    long long getGlobalTrackStep() const;

    /// step in a TRACK command
    void setLocalTrackStep(long long n);
    void incTrackSteps();
    long long getLocalTrackStep() const;

    void setNumBunch(short n);
    short getNumBunch() const;

    // used in ParallelCyclotronTracker for multi-bunch mode
    void setTotalNumPerBunch(size_t numpart, short n);
    size_t getTotalNumPerBunch(short n) const;

    void setLocalNumPerBunch(size_t numpart, short n);
    size_t getLocalNumPerBunch(short n) const;

    /* used in initializeTracking_m of ParallelCyclotronTracker
     * for multi-bunch mode
     */
    void countTotalNumPerBunch();

    void setGlobalMeanR(Vector_t globalMeanR);
    Vector_t getGlobalMeanR();
    void setGlobalToLocalQuaternion(Quaternion_t globalToLocalQuaternion);
    Quaternion_t getGlobalToLocalQuaternion();

    void setSteptoLastInj(int n);
    int getSteptoLastInj() const;

    /// calculate average angle of longitudinal direction of bins
    double calcMeanPhi();

    /// reset Bin[] for each particle according to the method given in paper PAST-AB(064402) by  G. Fubiani et al.
    bool resetPartBinID2(const double eta);

    bool resetPartBinBunch();

    ///@{ Access to reference data
    double getQ() const;
    double getM() const;
    double getP() const;
    double getE() const;
    ParticleOrigin getPOrigin() const;
    ParticleType getPType() const;
    double getInitialBeta() const;
    double getInitialGamma() const;
    ///@}
    ///@{ Set reference data
    void resetQ(double q);
    void resetM(double m);
    void setPOrigin(ParticleOrigin);
    void setPType(const std::string& type);
    ///@}
    double getdE() const;
    virtual double getGamma(int i);
    virtual double getBeta(int i);
    virtual void actT();

    const PartData* getReference() const;

    double getEmissionDeltaT();

    DistributionType getDistType() const;

    Quaternion_t getQKs3D();
    void         setQKs3D(Quaternion_t q);
    Vector_t     getKs3DRefr();
    void         setKs3DRefr(Vector_t r);
    Vector_t     getKs3DRefp();
    void         setKs3DRefp(Vector_t p);

    void iterateEmittedBin(int binNumber);

    void calcEMean();

    Inform& print(Inform& os);

    /*
     * (Pure) virtual member functions
     */

    virtual void runTests();

    virtual void do_binaryRepart();

    virtual void resetInterpolationCache(bool clearCache = false);

    //brief calculates back the max/min of the efield on the grid
    virtual VectorPair_t getEExtrema() = 0;

    virtual double getRho(int x, int y, int z) = 0;

    virtual void computeSelfFields() = 0;

    //brief used for self fields with binned distribution
    virtual void computeSelfFields(int bin) = 0;

    virtual void computeSelfFields_cycl(double gamma) = 0;
    virtual void computeSelfFields_cycl(int bin) = 0;

    virtual void swap(unsigned int i, unsigned int j);

    /*
       Mesh and Field Layout related functions
     */

    virtual void setBCAllPeriodic();
    virtual void setBCAllOpen();

    virtual void setBCForDCBeam();


//     virtual void setMesh(Mesh_t* mesh) = 0;
//     virtual Mesh_t &getMesh() = 0;

//     virtual void setFieldLayout(FieldLayout_t* fLayout) = 0;
    virtual FieldLayout_t& getFieldLayout() = 0;

    virtual void resizeMesh() { };

    /*
     * Wrapped member functions of IpplParticleBase
     */

    size_t getTotalNum() const;
    void setTotalNum(size_t n);
    void setLocalNum(size_t n);
    size_t getLocalNum() const;

    size_t getDestroyNum() const;
    size_t getGhostNum() const;

    ParticleLayout<T, Dim>& getLayout();
    const ParticleLayout<T, Dim>& getLayout() const;

    bool getUpdateFlag(UpdateFlags_t f) const;
    void setUpdateFlag(UpdateFlags_t f, bool val);


    ParticleBConds<Position_t, Dimension>& getBConds() {
        return pbase_m->getBConds();
    }

    void setBConds(const ParticleBConds<Position_t, Dimension>& bc) {
        pbase_m->setBConds(bc);
    }

    bool singleInitNode() const;

    void resetID();

    void update();
    void update(const ParticleAttrib<char>& canSwap);

    void createWithID(unsigned id);
    void create(size_t M);
    void globalCreate(size_t np);

    void destroy(size_t M, size_t I, bool doNow = false);
    void performDestroy(bool updateLocalNum = false);
    void ghostDestroy(size_t M, size_t I);

protected:
    size_t calcMoments();    // Calculates bunch moments using only emitted particles.

    /* Calculates bunch moments by summing over bins
     * (not accurate when any particles have been emitted).
     */
    void calcMomentsInitial();
    /// angle range [0~2PI) degree
    double calculateAngle(double x, double y);


private:
    virtual void updateDomainLength(Vektor<int, 3>& grid) = 0;

    virtual void updateFields(const Vector_t& hr, const Vector_t& origin);

    void setup(AbstractParticle<T, Dim>* pb);

public:
    /*
     * Bunch attributes
     */
    ParticlePos_t& R;
    ParticleIndex_t& ID;

    // Particle container attributes
    ParticleAttrib< Vector_t >     P;      // particle momentum //  ParticleSpatialLayout<double, 3>::ParticlePos_t P;
    ParticleAttrib< double >       Q;      // charge per simulation particle, unit: C.
    ParticleAttrib< double >       M;      // mass per simulation particle, for multi-species particle tracking, unit:GeV/c^2.
    ParticleAttrib< double >       Phi;    // the electric potential
    ParticleAttrib< Vector_t >     Ef;     // e field vector
    ParticleAttrib< Vector_t >     Eftmp;  // e field vector for gun simulations

    ParticleAttrib< Vector_t >     Bf;     // b field vector
    ParticleAttrib< int >          Bin;    // holds the bin in which the particle is in, if zero particle is marked for deletion
    ParticleAttrib< double >       dt;     // holds the dt timestep for particle
    ParticleAttrib< ParticleType > PType;  // particle names
    ParticleAttrib< ParticleOrigin > POrigin;  // we can distinguish dark current particles from primary particle
    ParticleAttrib< int >          TriID;  // holds the ID of triangle that the particle hit. Only for BoundaryGeometry case.
    ParticleAttrib< short >        cavityGapCrossed; // particle just crossed cavity gap (for ParallelCyclotronTracker)
    ParticleAttrib< short >        bunchNum; // bunch number to which particle belongs (multi-bunch mode)

    Vector_t RefPartR_m;
    Vector_t RefPartP_m;

    CoordinateSystemTrafo toLabTrafo_m;

    ParticleOrigin refPOrigin_m;
    ParticleType refPType_m;

    // The structure for particle binning
    PartBins* pbin_m;

    /// timer for IC, can not be in Distribution.h
    IpplTimings::TimerRef distrReload_m;
    IpplTimings::TimerRef distrCreate_m;

    // For AMTS integrator in OPAL-T
    double dtScInit_m, deltaTau_m;

    // get 2nd order momentum matrix
    FMatrix<double, 2 * Dim, 2 * Dim> getSigmaMatrix() const;

private:
    // save particles in case of one core
    std::unique_ptr<Inform> pmsg_m;
    std::unique_ptr<std::ofstream> f_stream;
    /// if the grid does not have to adapt
    bool fixed_grid;

protected:
    IpplTimings::TimerRef boundpTimer_m;
    IpplTimings::TimerRef boundpBoundsTimer_m;
    IpplTimings::TimerRef boundpUpdateTimer_m;
    IpplTimings::TimerRef statParamTimer_m;

    IpplTimings::TimerRef histoTimer_m;
    /// timer for selfField calculation
    IpplTimings::TimerRef selfFieldTimer_m;

    const PartData* reference;

    /*
       Member variables starts here
     */

    // unit state of PartBunch
    UnitState_t unit_state_;
    UnitState_t stateOfLastBoundP_;

    /// holds the centroid of the beam
    double centroid_m[2 * Dim];

    /// holds the timestep in seconds
    double dt_m;
    /// holds the actual time of the integration
    double t_m;
    /// the position along design trajectory
    double spos_m;

    /// Initialize the translation vector and rotation quaternion
    /// here. Cyclotron tracker will reset these values each timestep
    /// TTracker can just use 0 translation and 0 rotation (quat[1 0 0 0]).
    //Vector_t globalMeanR_m = Vector_t(0.0, 0.0, 0.0);
    //Quaternion_t globalToLocalQuaternion_m = Quaternion_t(1.0, 0.0, 0.0, 0.0);
    Vector_t globalMeanR_m;
    Quaternion_t globalToLocalQuaternion_m;

    /// maximal extend of particles
    Vector_t rmax_m;
    /// minimal extend of particles
    Vector_t rmin_m;

    //RMS number density of particles from grid
    double rmsDensity_m;

    /// meshspacing of cartesian mesh
    Vector_t hr_m;
    /// meshsize of cartesian mesh
    Vektor<int, 3> nr_m;

    /// stores the used field solver
    FieldSolver* fs_m;

    double couplingConstant_m;

    double qi_m;
    double massPerParticle_m;

    /// counter to store the distribution dump
    int distDump_m;

    /// Mesh enlargement
    double dh_m; /// relative enlargement of the mesh

    /// if larger than 0, emitt particles for tEmission_m [s]
    double tEmission_m;

    /// holds the gamma of the bin
    std::unique_ptr<double[]> bingamma_m;

    //FIXME: this should go into the Bin class!
    // holds number of emitted particles of the bin
    // jjyang: opal-cycl use *nBin_m of pbin_m
    std::unique_ptr<size_t[]> binemitted_m;

    /// steps per turn for OPAL-cycl
    int stepsPerTurn_m;

    /// step in a TRACK command
    long long localTrackStep_m;

    /// if multiple TRACK commands
    long long globalTrackStep_m;

    /// current bunch number
    short numBunch_m;

    /// number of particles per bunch
    std::vector<size_t> bunchTotalNum_m;
    std::vector<size_t> bunchLocalNum_m;

    /// this parameter records the current steps since last bunch injection
    /// it helps to inject new bunches correctly in the restart run of OPAL-cycl
    /// it is stored during phase space dump.
    int SteptoLastInj_m;

    /*
      Data structure for particle load balance information
    */

    std::unique_ptr<size_t[]> globalPartPerNode_m;

    Distribution *dist_m;
    DistributionMoments momentsComputer_m;

    // flag to tell if we are a DC-beam
    bool dcBeam_m;
    double periodLength_m;
    std::shared_ptr<AbstractParticle<T, Dim> > pbase_m;
};

template<class T, unsigned Dim>
typename PartBunchBase<T, Dim>::ConstIterator begin(PartBunchBase<T, Dim> const& bunch) {
    return bunch.begin();
}

template<class T, unsigned Dim>
typename PartBunchBase<T, Dim>::ConstIterator end(PartBunchBase<T, Dim> const& bunch) {
    return bunch.end();
}

#include "PartBunchBase.hpp"

#endif
