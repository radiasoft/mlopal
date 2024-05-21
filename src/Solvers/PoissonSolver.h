#ifndef POISSON_SOLVER_H_
#define POISSON_SOLVER_H_

//////////////////////////////////////////////////////////////
#include "Algorithms/PBunchDefs.h"

#include "Field/Field.h"

#ifdef ENABLE_AMR
    #include "Utilities/OpalException.h"
#endif
//////////////////////////////////////////////////////////////
template <class T, unsigned Dim>
class PartBunchBase;

class PoissonSolver {

protected:
    typedef Field<int, 3, Mesh_t, Center_t>                  IField_t;
    typedef Field<std::complex<double>, 3, Mesh_t, Center_t> CxField_t;

public:

    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential in open space
    virtual void computePotential(Field_t &rho, Vector_t hr) = 0;
    
#ifdef ENABLE_AMR
    /**
     * AMR solver calls
     * 
     * @param rho right-hand side charge density on grid [C / m]
     * @param phi electrostatic potential (unknown) [V]
     * @param efield electric field [V / m]
     * @param baseLevel for solve
     * @param finestLevel for solve
     * @param prevAsGuess use of previous solution as initial guess
     */
    virtual void solve(AmrScalarFieldContainer_t &/*rho*/,
                       AmrScalarFieldContainer_t &/*phi*/,
                       AmrVectorFieldContainer_t &/*efield*/,
                       unsigned short /*baseLevel*/,
                       unsigned short /*finestLevel*/,
                       bool /*prevAsGuess*/ = true)
    {
        throw OpalException("PoissonSolver::solve()", "Not supported for non-AMR code.");
    };
    
    /**
     * Tell solver to regrid
     */
    virtual void hasToRegrid() {
        throw OpalException("PoissonSolver::hasToRegrid()", "Not supported for non-AMR code.");
    }
#endif
                                  
    virtual void computePotential(Field_t &rho, Vector_t hr, double zshift) = 0;

    virtual double getXRangeMin(unsigned short level = 0) = 0;
    virtual double getXRangeMax(unsigned short level = 0) = 0;
    virtual double getYRangeMin(unsigned short level = 0) = 0;
    virtual double getYRangeMax(unsigned short level = 0) = 0;
    virtual double getZRangeMin(unsigned short level = 0) = 0;
    virtual double getZRangeMax(unsigned short level = 0) = 0;
    virtual void test(PartBunchBase<double, 3> *bunch) = 0 ;
    virtual ~PoissonSolver(){};

    virtual void resizeMesh(Vector_t& /*origin*/, Vector_t& /*hr*/,
                            const Vector_t& /*rmin*/, const Vector_t& /*rmax*/,
                            double /*dh*/)
    { };
    virtual double getinteractionRadius() const { return 0.0; }

    virtual void calculatePairForces(PartBunchBase<double, 3>* /*bunch*/, 
                                     double /*gammaz*/) 
    { };

};

inline Inform &operator<<(Inform &os, const PoissonSolver &/*fs*/) {
    return os << "";
}
#endif
