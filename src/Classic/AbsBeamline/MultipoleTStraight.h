/*
 *  Copyright (c) 2017, Titus Dascalu
 *  Copyright (c) 2018, Martin Duy Tat
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */


#ifndef CLASSIC_MULTIPOLET_STRAIGHT_H
#define CLASSIC_MULTIPOLET_STRAIGHT_H

/** ---------------------------------------------------------------------
  *
  * MultipoleTiStraight defines a straight combined function magnet (up 
  * to arbitrary multipole component) with fringe fields
  *
  * ---------------------------------------------------------------------
  *
  * Class category: AbsBeamline \n
  * $Author: Titus Dascalu, Martin Duy Tat, Chris Rogers
  *
  * ---------------------------------------------------------------------
  *
  * The field is obtained from the scalar potential \n
  *     @f[ V = f_0(x,s) z + f_1 (x,s) \frac{z^3}{3!} + f_2 (x,s) \frac{z^5}{5!}
  *      + ...  @f] \n
  *     (x,z,s) -> Frenet-Serret local coordinates along the magnet \n
  *     z -> vertical component \n
  *     assume mid-plane symmetry \n 
  *     set field on mid-plane -> @f$ B_z = f_0(x,s) = T(x) \cdot S(s) @f$ \n
  *     T(x) -> transverse profile; this is a polynomial describing
  *             the field expansion on the mid-plane inside the magnet
  *             (not in the fringe field);
  *             1st term is the dipole strength, 2nd term is the 
  *             quadrupole gradient * x, etc. \n
  *          -> when setting the magnet, one gives the multipole
  *             coefficients of this polynomial (i.e. dipole strength,  
  *             quadrupole gradient, etc.) \n
  * \n
  * ------------- example ----------------------------------------------- \n
  *     Setting a combined function magnet with dipole, quadrupole and 
  *     sextupole components: \n
  *     @f$ T(x) = B_0 + B_1 \cdot x + B_2 \cdot x^2 @f$\n
  *     user gives @f$ B_0, B_1, B_2 @f$ \n
  * ------------- example end ------------------------------------------- \n
  * \n
  *     S(s) -> fringe field \n
  *     recursion -> @f$ f_n (x,s) = (-1)^n \cdot \sum_{i=0}^{n} C_n^i 
  *     \cdot T^{(2i)} \cdot S^{(2n-2i)} @f$ \n
  *     for curved magnets the above recursion is more complicated \n
  *     @f$ C_n^i @f$ -> binomial coeff; 
  *     @f$ T^{(n)} @f$ -> n-th derivative
  *
  * ---------------------------------------------------------------------
  */

#include "BeamlineGeometry/StraightGeometry.h"
#include "AbsBeamline/MultipoleTBase.h"
#include <vector>

class MultipoleTStraight: public MultipoleTBase {
public: 
    /** Constructor
     *  \param name -> User-defined name
     */
    explicit MultipoleTStraight(const std::string &name);
    /** Copy constructor */
    MultipoleTStraight(const MultipoleTStraight &right);
    /** Destructor */ 
    ~MultipoleTStraight();
    /** Inheritable copy constructor */
    virtual ElementBase* clone() const override;
    /** Accept a beamline visitor */
    void accept(BeamlineVisitor &visitor) const override;
    /** Set the number of terms used in calculation of field components \n
     *  Maximum power of z in Bz is 2 * maxOrder_m
     *  \param maxOrder -> Number of terms in expansion in z
     */
    virtual void setMaxOrder(const std::size_t &maxOrder) override;
    /** Return the cell geometry */
    StraightGeometry& getGeometry() override;
    /** Return the cell geometry */
    const StraightGeometry& getGeometry() const override;
    /** Initialise the MultipoleT
     *  \param bunch -> Bunch the global bunch object
     *  \param startField -> Not used
     *  \param endField -> Not used
     */
    virtual void initialise(PartBunchBase<double, 3>* bunch,
                            double &startField,
                            double &endField) override;
private:
    MultipoleTStraight operator=(const MultipoleTStraight &rhs);
    /** Geometry */
    StraightGeometry straightGeometry_m;
    /** Transform to Frenet-Serret coordinates for sector magnets */
    virtual void transformCoords(Vector_t &R) override;
    /** Transform B-field from Frenet-Serret coordinates to lab coordinates */
    virtual void transformBField(Vector_t &B, const Vector_t &R) override;
    /** Radius of curvature \n
     *  Straight magnet, infinite radius, infinity (1.0e300) is returned
     *  \param s -> Coordinate s
     */
    virtual double getRadius(const double &s) override;
    /** Returns the scale factor @f$ h_s = 1@f$
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    virtual double getScaleFactor(const double &x, const double &s) override;
    /** Get x-component of the B-field \n
     *  This function has been overloaded because calculating \n
     *  the B-field directly is quicker and more accurate
     */
    virtual double getBx (const Vector_t &R) override;
    /** Get s-component of the B-field \n
     *  This function has been overloaded because calculating \n
     *  the B-field directly is quicker and more accurate
     */
    virtual double getBs (const Vector_t &R) override;
    /** Calculate fn(x, s) by expanding the differential operator
     *  (from Laplacian and scalar potential) in terms of polynomials
     *  \param n -> nth derivative
     *  \param x -> Coordinate x
     *  \param s -> Coordinate s
     */
    virtual double getFn(const std::size_t &n,
                         const double &x,
                         const double &s) override;
};

inline
    void MultipoleTStraight::accept(BeamlineVisitor &visitor) const {
        visitor.visitMultipoleTStraight(*this);
}
inline
    void MultipoleTStraight::transformBField(Vector_t &/*B*/, const Vector_t &/*R*/) {
}
inline
    double MultipoleTStraight::getRadius(const double &/*s*/) {
        return 1e300;
}
inline
    double MultipoleTStraight::getScaleFactor(const double &/*x*/,
                                              const double &/*s*/) {
    return 1.0;
}
inline
    StraightGeometry& MultipoleTStraight::getGeometry() {
        return straightGeometry_m;
}
inline
    const StraightGeometry& MultipoleTStraight::getGeometry() const {
        return straightGeometry_m;
}
inline
    void MultipoleTStraight::initialise(PartBunchBase<double, 3>* bunch,
                                        double &/*startField*/,
                                        double &/*endField*/) {
        RefPartBunch_m = bunch;
        straightGeometry_m.setElementLength(2 * getBoundingBoxLength());
}

#endif
