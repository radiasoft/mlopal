//
// Class ParallelCyclotronTracker
//   Tracker for OPAL-Cycl
//
// Copyright (c) 2007 - 2014, Jianjun Yang, Andreas Adelmann and Matthias Toggweiler,
//                            Paul Scherrer Institut, Villigen PSI, Switzerland
// Copyright (c) 2014,        Daniel Winklehner, MIT, Cambridge, MA, USA
// Copyright (c) 2012 - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects"
// and the paper
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects:
// Model, implementation, and application"
// (https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.13.064201)
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
#ifndef OPAL_ParallelCyclotronTracker_HH
#define OPAL_ParallelCyclotronTracker_HH

#include "AbsBeamline/ElementBase.h"
#include "Algorithms/MultiBunchHandler.h"
#include "Algorithms/Tracker.h"
#include "Steppers/Steppers.h"

#include <memory>
#include <tuple>
#include <vector>

class DataSink;
class PluginElement;
class LossDataSink;

template <class T, unsigned Dim>
class PartBunchBase;

struct CavityCrossData {
    RFCavity* cavity;
    double sinAzimuth;
    double cosAzimuth;
    double perpenDistance;
};

class ParallelCyclotronTracker: public Tracker {

public:
    typedef std::vector<double> dvector_t;
    typedef std::vector<int> ivector_t;
    typedef std::pair<double[8], Component*>      element_pair;
    typedef std::pair<ElementType, element_pair>        type_pair;
    typedef std::list<type_pair*>                 beamline_list;

    /// Constructor.
    //  The beam line to be tracked is "bl".
    //  The particle reference data are taken from "data".
    //  The particle bunch tracked is taken from [b]bunch[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    ParallelCyclotronTracker(const Beamline& bl, PartBunchBase<double, 3>* bunch, DataSink& ds,
                             const PartData& data, bool revBeam,
                             bool revTrack, int maxSTEPS,
                             Steppers::TimeIntegrator timeintegrator,
                             const int& numBunch,
                             const double& mbEta,
                             const double& mbPara,
                             const std::string& mbMode,
                             const std::string& mbBinning);

    virtual ~ParallelCyclotronTracker();

    /// Apply the algorithm to the top-level beamline.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void execute();

    /// Calculate the field map at an arbitrary point
    //  R - position
    //  P - momentum; note the field is returned in the lab frame
    //  t - time
    //  Efield - filled with values for the electric field
    //  Bfield - filled with values for the magnetic field
    //  Returns a boolean value, true if the particle was out of the nominal
    //  field map boundary, else false.
    bool computeExternalFields_m(const Vector_t& R,
                                 const Vector_t& P,
                                 const double& t,
                                 Vector_t& Efield,
                                 Vector_t& Bfield);

    /// Apply the algorithm to a beam line.
    //  overwrite the execute-methode from DefaultVisitor
    virtual void visitBeamline(const Beamline&);
    
    /// Apply the algorithm to a collimator.
    virtual void visitCCollimator(const CCollimator&);

    /// Apply the algorithm to a closed orbit corrector.
    virtual void visitCorrector(const Corrector&);

    /// Apply the algorithm to a cyclotron
    virtual void visitCyclotron(const Cyclotron& cycl);

    /// Apply the algorithm to a degrader
    virtual void visitDegrader(const Degrader&);

    /// Apply the algorithm to a drift space.
    virtual void visitDrift(const Drift&);

    /// Apply the algorithm to a flexible collimator
    virtual void visitFlexibleCollimator(const FlexibleCollimator&);

    /// Apply the algorithm to a marker.
    virtual void visitMarker(const Marker&);

    /// Apply the algorithm to a beam position monitor.
    virtual void visitMonitor(const Monitor&);

    /// Apply the algorithm to a multipole.
    virtual void visitMultipole(const Multipole&);

    /// Apply the algorithm to an arbitrary multipole.
    virtual void visitMultipoleT(const MultipoleT&);

    /// Apply the algorithm to an arbitrary straight multipole.
    virtual void visitMultipoleTStraight(const MultipoleTStraight&);

    /// Apply the algorithm to an arbitrary curved multipole of constant radius.
    virtual void visitMultipoleTCurvedConstRadius(const MultipoleTCurvedConstRadius&);

    /// Apply the algorithm to an arbitrary curved multipole of variable radius.
    virtual void visitMultipoleTCurvedVarRadius(const MultipoleTCurvedVarRadius&);

    /// Apply the algorithm to a offset (placement).
    virtual void visitOffset(const Offset&);

    /// Apply the algorithm to a probe.
    virtual void visitProbe(const Probe&);

    /// Apply the algorithm to a rectangular bend.
    virtual void visitRBend(const RBend&);

    /// Apply the algorithm to a RF cavity.
    virtual void visitRFCavity(const RFCavity&);

    /// Apply the algorithm to a ring
    virtual void visitRing(const Ring& ring);
    
    /// Apply the algorithm to a sector bend.
    virtual void visitSBend(const SBend&);

    /// Apply the algorithm to a sector bend with 3D field map.
    virtual void visitSBend3D(const SBend3D&);

    /// Apply the algorithm to a scaling FFA magnet.
    virtual void visitScalingFFAMagnet(const ScalingFFAMagnet& bend);

    /// Apply the algorithm to a septum.
    virtual void visitSeptum(const Septum&);

    /// Apply the algorithm to a solenoid.
    virtual void visitSolenoid(const Solenoid&);

    /// Apply the algorithm to a particle stripper.
    virtual void visitStripper(const Stripper&);

    /// Apply the algorithm to a vacuum space.
    virtual void visitVacuum(const Vacuum&);

    /// Apply the algorithm to a variabel RF cavity.
    virtual void visitVariableRFCavity(const VariableRFCavity& cav);

    /// Apply the algorithm to a variable RF cavity with Fringe Field.
    virtual void visitVariableRFCavityFringeField(const VariableRFCavityFringeField& cav);

    /// Apply the algorithm to a vertical FFA magnet.
    virtual void visitVerticalFFAMagnet(const VerticalFFAMagnet& bend);

    /// set last dumped step
    inline void setLastDumpedStep(const int para) { lastDumpedStep_m = para; }

    ///@{ Method for restart
    inline void setPr(double x) {referencePr = x;}
    inline void setPt(double x) {referencePt = x;}
    inline void setPz(double x) {referencePz = x;}
    inline void setR(double x) {referenceR = x;}
    inline void setTheta(double x) {referenceTheta = x;}
    inline void setZ(double x) {referenceZ = x;}
    inline void setBeGa(double x) {bega = x;}
    inline void setPhi(double x) {referencePhi = x;}
    inline void setPsi(double x) {referencePsi = x;}
    inline void setPreviousH5Local(bool x) {previousH5Local = x;}
    ///@}
    void bgf_main_collision_test();
    void initializeBoundaryGeometry();

private:
    enum class TrackingMode: unsigned short {
        UNDEFINED,
        SINGLE,
        SEO,
        BUNCH
    };

    // Not implemented.
    ParallelCyclotronTracker();
    ParallelCyclotronTracker(const ParallelCyclotronTracker&);
    void operator=(const ParallelCyclotronTracker&);

    beamline_list FieldDimensions;
    std::list<Component*> myElements;
    Beamline* itsBeamline;
    std::vector<PluginElement*> pluginElements_m;
    std::vector<CavityCrossData> cavCrossDatas_m;

    DataSink* itsDataSink;

    BoundaryGeometry* bgf_m;

    Cyclotron* cycl_m;

    /// The maximal number of steps the system is integrated
    int maxSteps_m;

    /// The positive axes unit vectors
    static Vector_t const xaxis;
    static Vector_t const yaxis;
    static Vector_t const zaxis;

    /// The reference variables
    double bega;
    double referenceR;
    double referenceTheta;
    double referenceZ = 0.0;

    double referencePr;
    double referencePt;
    double referencePz = 0.0;
    double referencePtot;

    double referencePsi;
    double referencePhi;

    bool spiral_flag = false;

    Vector_t PreviousMeanP;

    bool previousH5Local;

    double sinRefTheta_m;
    double cosRefTheta_m;

    // only used for multi-bunch mode
    std::unique_ptr<MultiBunchHandler> mbHandler_m;

    int lastDumpedStep_m;

    // for each bunch we have a path length
    double pathLength_m;
    // Multi time step tracker
    void MtsTracker();

    void GenericTracker();
    bool getFieldsAtPoint(const double& t, const size_t& Pindex, Vector_t& Efield, Vector_t& Bfield);

    /*
      Local Variables both used by the integration methods
    */

    long long step_m;
    long long restartStep0_m;

    int turnnumber_m;

    // only for dumping
    double azimuth_m;
    double prevAzimuth_m;

    /* only for dumping to stat-file --> make azimuth monotonically increasing
     * @param theta computed azimuth [deg]
     * @param prevAzimuth previous angle [deg]
     * @param azimuth to be updated [deg]
     * @param bunchNr in restart mode only --> to compute initial azimuth
     */
    void dumpAngle(const double& theta,
                   double& prevAzimuth,
                   double& azimuth,
                   const short& bunchNr = 0);

    double computeRadius(const Vector_t& meanR) const;

    void computePathLengthUpdate(std::vector<double>& dl, const double& dt);

    // external field arrays for dumping
    Vector_t FDext_m[2], extE_m, extB_m;

    const int myNode_m;
    const size_t initialLocalNum_m;
    const size_t initialTotalNum_m;

    /// output coordinates at different azimuthal angles and one after every turn
    std::vector<std::unique_ptr<std::ofstream> > outfTheta_m;
    /// the different azimuthal angles for the outfTheta_m output files
    std::vector<double> azimuth_angle_m;
    ///@ open / close output coordinate files
    void openFiles(size_t numFiles, std::string fn);
    void closeFiles();
    ///@}
    /// output file for six dimensional phase space
    std::ofstream outfTrackOrbit_m;

    void buildupFieldList(double BcParameter[], ElementType elementType, Component* elptr);

    // angle range [0~2PI) degree
    double calculateAngle(double x, double y);
    // angle range [-PI~PI) degree
    double calculateAngle2(double x, double y);

    bool checkGapCross(Vector_t Rold, Vector_t Rnew, RFCavity* rfcavity, double& DistOld);
    bool RFkick(RFCavity* rfcavity, const double t, const double dt, const int Pindex);

    bool getTunes(dvector_t& t, dvector_t& r, dvector_t& z, int lastTurn, double& nur, double& nuz);

    IpplTimings::TimerRef IntegrationTimer_m;
    IpplTimings::TimerRef DumpTimer_m ;
    IpplTimings::TimerRef TransformTimer_m;
    IpplTimings::TimerRef BinRepartTimer_m;
    IpplTimings::TimerRef PluginElemTimer_m;
    IpplTimings::TimerRef DelParticleTimer_m;

    /*
     * @param bunchNr if <= -1 --> take all particles in simulation, if bunchNr > -1,
     * take only particles of *this* bunch
     */
    Vector_t calcMeanR(short bunchNr = -1) const;

    Vector_t calcMeanP() const;

    void repartition(); // Do repartition between nodes if step_m is multiple of Options::repartFreq

    // Transform the x- and y-parts of a particle attribute (position, momentum, fields) from the
    // global reference frame to the local reference frame.

    // phi is the angle of the bunch measured counter-clockwise from the positive x-axis.
    void globalToLocal(ParticleAttrib<Vector_t>& vectorArray,
                       double phi, Vector_t const translationToGlobal = 0);

    // Transform the x- and y-parts of a particle attribute (position, momentum, fields) from the
    // local reference frame to the global reference frame.
    void localToGlobal(ParticleAttrib<Vector_t>& vectorArray,
                       double phi, Vector_t const translationToGlobal = 0);

    // Overloaded version of globalToLocal using a quaternion for 3D rotation
    inline void globalToLocal(ParticleAttrib<Vector_t>& vectorArray,
                              Quaternion_t const quaternion,
                              Vector_t const meanR = Vector_t(0.0));

    // Overloaded version of localToGlobal using a quaternion for 3D rotation
    inline void localToGlobal(ParticleAttrib<Vector_t>& vectorArray,
                              Quaternion_t const quaternion,
                              Vector_t const meanR = Vector_t(0.0));

    // Overloaded version of globalToLocal using phi and theta for pseudo 3D rotation
    inline void globalToLocal(ParticleAttrib<Vector_t>& particleVectors,
                              double const phi, double const psi,
                              Vector_t const meanR = Vector_t(0.0));

    // Overloaded version of localToGlobal using phi and theta for pseudo 3D rotation
    inline void localToGlobal(ParticleAttrib<Vector_t>& particleVectors,
                              double const phi, double const psi,
                              Vector_t const meanR = Vector_t(0.0));

    // Overloaded version of globalToLocal using phi and theta for pseudo 3D rotation, single vector
    inline void globalToLocal(Vector_t& myVector,
                              double const phi, double const psi,
                              Vector_t const meanR = Vector_t(0.0));

    // Overloaded version of localToGlobal using phi and theta for pseudo 3D rotation, single vector
    inline void localToGlobal(Vector_t& myVector,
                              double const phi, double const psi,
                              Vector_t const meanR = Vector_t(0.0));

    // Rotate the particles by an angle and axis defined in the quaternion (4-vector w,x,y,z)
    inline void rotateWithQuaternion(ParticleAttrib<Vector_t>& vectorArray, Quaternion_t const quaternion);

    // Returns the quaternion of the rotation from vector u to vector v
    inline void getQuaternionTwoVectors(Vector_t u, Vector_t v, Quaternion_t& quaternion);

    // Normalization of a quaternion
    inline void normalizeQuaternion(Quaternion_t& quaternion);

    // Normalization of a quaternion
    inline void normalizeVector(Vector_t& vector);

    // Some rotations
    inline void rotateAroundZ(ParticleAttrib<Vector_t>& particleVectors, double const phi);
    inline void rotateAroundX(ParticleAttrib<Vector_t>& particleVectors, double const psi);
    inline void rotateAroundZ(Vector_t& myVector, double const phi);
    inline void rotateAroundX(Vector_t& myVector, double const psi);

    // Push particles for time h.
    // Apply effects of RF Gap Crossings.
    // Unit assumptions: [itsBunch_m->R] = m, [itsBunch_m->P] = 1, [h] = s, [c] = m/s, [itsBunch_m->getT()] = s
    bool push(double h);
    // Kick particles for time h
    // The fields itsBunch_m->Bf, itsBunch_m->Ef are used to calculate the forces
    bool kick(double h);

    // Apply the trilogy half push - kick - half push,
    // considering only external fields
    void borisExternalFields(double h);

    // apply the plugin elements: probe, collimator, stripper, septum
    bool applyPluginElements(const double dt);

    // destroy particles if they are marked as Bin=-1 in the plugin elements or out of global aperture
    bool deleteParticle(bool = false);

    void initTrackOrbitFile();

    void singleParticleDump();

    void bunchDumpStatData();

    void bunchDumpPhaseSpaceData();

    void evaluateSpaceChargeField();

    void initDistInGlobalFrame();

    void setTimeStep();

    void checkFileMomentum();

    void checkNumPart(std::string s);

    // we store a pointer explicitly to the Ring
    Ring* opalRing_m;

    // to save geometry losses
    std::unique_ptr<LossDataSink> lossDs_m;

    // If Ring is defined take the harmonic number from Ring; else use
    // cyclotron
    double getHarmonicNumber() const;

    typedef std::function<bool(const double&,
                               const size_t&,
                               Vector_t&,
                               Vector_t&)> function_t;

    std::unique_ptr< Stepper<function_t> > itsStepper_mp;

    struct settings {
        int scSolveFreq;
        int stepsPerTurn;
        double deltaTheta;
        int stepsNextCheck;
    } setup_m;

    TrackingMode mode_m;

    Steppers::TimeIntegrator stepper_m;

    /// Check if turn done
    bool isTurnDone();

    /// Update time and path length, write to output files
    void update_m(double& t, const double& dt, const bool& finishedTurn);

    /*!
     * @returns the time t [ns], time step dt [ns] and the azimuth angle [rad]
     */
    std::tuple<double, double, double> initializeTracking_m();

    void finalizeTracking_m(dvector_t& Ttime,
                            dvector_t& Tdeltr,
                            dvector_t& Tdeltz,
                            ivector_t& TturnNumber);

    void setTrackingMode();

    void seoMode_m(double& t, const double dt, bool& finishedTurn,
                   dvector_t& Ttime, dvector_t& Tdeltr,
                   dvector_t& Tdeltz, ivector_t& TturnNumber);

    void singleMode_m(double& t, const double dt, bool& finishedTurn, double& oldReferenceTheta);

    void bunchMode_m(double& t, const double dt, bool& finishedTurn);

    void gapCrossKick_m(size_t i, double t, double dt,
                        const Vector_t& Rold, const Vector_t& Pold);


    inline void dumpAzimuthAngles_m(const double& t,
                                    const Vector_t& R,
                                    const Vector_t& P,
                                    const double& oldReferenceTheta,
                                    const double& temp_meanTheta);

    inline void dumpThetaEachTurn_m(const double& t,
                                    const Vector_t& R,
                                    const Vector_t& P,
                                    const double& temp_meanTheta,
                                    bool& finishedTurn);

    void computeSpaceChargeFields_m();

    bool computeExternalFields_m(const size_t& i,
                                 const double& t,
                                 Vector_t& Efield,
                                 Vector_t& Bfield);

    void injectBunch(bool& flagTransition);

    void saveInjectValues();

    bool isMultiBunch() const;

    bool hasMultiBunch() const;

    /*
     * @param dt time step in ns
     */
    void updatePathLength(const double& dt);

    /*
     * @param dt time step in ns
     */
    void updateTime(const double& dt);

    void updateAzimuthAndRadius();

    /** multi-bunch mode: set the path length of each bunch in case of restart mode
     *
     * At creation of DataSink the lines are rewinded properly --> the last entry of
     * the path length is therefore the initial path length at restart.
     * @pre In order to work properly in restart mode, the lines in the multi-bunch
     * *.smb files need to be rewinded first.
     */
    void initPathLength();
};

/**
 * @param x
 * @param y
 *
 * @return angle range [0~2PI) degree
 */
inline
double ParallelCyclotronTracker::calculateAngle(double x, double y) {
    double thetaXY = std::atan2(y, x);

    if (thetaXY < 0) return thetaXY + Physics::two_pi;
    return thetaXY;
}

/**
 * @param x
 * @param y
 *
 * @return angle range [-PI~PI) degree
 */
inline
double ParallelCyclotronTracker::calculateAngle2(double x, double y) {
    return std::atan2(y,x);
}


inline
bool ParallelCyclotronTracker::isMultiBunch() const {
    return ( (mbHandler_m != nullptr) && itsBunch_m->weHaveBins() );
}


inline
bool ParallelCyclotronTracker::hasMultiBunch() const {
    return ( isMultiBunch() && mbHandler_m->getNumBunch() > 1 );
}

#endif // OPAL_ParallelCyclotronTracker_HH
