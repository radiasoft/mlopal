/* 
 *  Copyright (c) 2012, Chris Rogers
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

#ifndef OPAL_OPALSBEND3D_HH
#define OPAL_OPALSBEND3D_HH

#include "Elements/OpalBend.h"

/** OpalSBend3D provides user interface information for the SBEND3D object
 *
 *  Defines three parameters - field map name, units for field, length for field
 */
class OpalSBend3D: public OpalElement {
  public:
    /** enum maps string to integer value for UI definitions */
    enum {
        FMAPFN = COMMON,
        FIELD_UNITS,
        LENGTH_UNITS,
        SIZE // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalSBend3D();

    /** Destructor does nothing */
    virtual ~OpalSBend3D();

    /** Inherited copy constructor */
    virtual OpalSBend3D *clone(const std::string &name);

    /** Update the SBend3D with new parameters from UI parser */
    virtual void update();

  private:
    // Not implemented.
    OpalSBend3D(const OpalSBend3D &);
    void operator=(const OpalSBend3D &);

    // Clone constructor.
    OpalSBend3D(const std::string &name, OpalSBend3D *parent);
};

#endif // OPAL_OpalSBend_HH
