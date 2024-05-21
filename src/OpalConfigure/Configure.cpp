//
// Namespace Configure
//   The OPAL configurator.
//   This class must be modified to configure the commands to be contained
//   in an executable OPAL program. For each command an exemplar object
//   is constructed and linked to the main directory. This exemplar is then
//   available to the OPAL parser for cloning.
//   This class could be part of the class OpalData.  It is separated from
//   that class and opale into a special module in order to reduce
//   dependencies between modules.
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
#include "OpalConfigure/Configure.h"
#include "AbstractObjects/OpalData.h"

#include "Distribution/Distribution.h"

// Basic action commands.
#include "BasicActions/Call.h"
#include "BasicActions/DumpFields.h"
#include "BasicActions/DumpEMFields.h"
#include "BasicActions/Echo.h"
#include "BasicActions/Help.h"
#include "BasicActions/Option.h"
#include "BasicActions/Select.h"
#include "BasicActions/Stop.h"
#include "BasicActions/Quit.h"
#include "BasicActions/System.h"
#include "BasicActions/PSystem.h"
#include "BasicActions/Title.h"
#include "BasicActions/Value.h"

// Macro command.
#include "OpalParser/MacroCmd.h"

// Commands introducing a special mode.
#include "Track/TrackCmd.h"

// Table-related commands.
#include "Structure/Beam.h"
#include "Structure/FieldSolver.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/OpalWake.h"
#include "Structure/ParticleMatterInteraction.h"
#include "Utilities/OpalFilter.h"
#include "TrimCoils/OpalTrimCoil.h"
#include "Tables/List.h"

// Value definitions commands.
#include "ValueDefinitions/BoolConstant.h"
#include "ValueDefinitions/RealConstant.h"
#include "ValueDefinitions/RealVariable.h"
#include "ValueDefinitions/RealVector.h"
#include "ValueDefinitions/StringConstant.h"

// Element commands.
#include "Elements/OpalAsymmetricEnge.h"
#include "Elements/OpalCavity.h"
#include "Elements/OpalCCollimator.h"
#include "Elements/OpalCyclotron.h"
#include "Elements/OpalDegrader.h"
#include "Elements/OpalDrift.h"
#include "Elements/OpalECollimator.h"
#include "Elements/OpalEnge.h"
#include "Elements/OpalFlexibleCollimator.h"
#include "Elements/OpalHKicker.h"
#include "Elements/OpalKicker.h"
#include "Elements/OpalMarker.h"
#include "Elements/OpalMonitor.h"
#include "Elements/OpalMultipole.h"
#include "Elements/OpalMultipoleT.h"
#include "Elements/OpalMultipoleTStraight.h"
#include "Elements/OpalMultipoleTCurvedConstRadius.h"
#include "Elements/OpalMultipoleTCurvedVarRadius.h"
#include "Elements/OpalOctupole.h"
#include "Elements/OpalOffset/OpalLocalCartesianOffset.h"
#include "Elements/OpalOffset/OpalLocalCylindricalOffset.h"
#include "Elements/OpalOffset/OpalGlobalCartesianOffset.h"
#include "Elements/OpalOffset/OpalGlobalCylindricalOffset.h"
#include "Elements/OpalPepperPot.h"
#include "Elements/OpalPolynomialTimeDependence.h"
#include "Elements/OpalProbe.h"
#include "Elements/OpalQuadrupole.h"
#include "Elements/OpalRBend.h"
#include "Elements/OpalRBend3D.h"
#include "Elements/OpalRCollimator.h"
#include "Elements/OpalRingDefinition.h"
#include "Elements/OpalSBend.h"
#include "Elements/OpalSBend3D.h"
#include "Elements/OpalScalingFFAMagnet.h"
#include "Elements/OpalSeptum.h"
#include "Elements/OpalSextupole.h"
#include "Elements/OpalSlit.h"
#include "Elements/OpalSolenoid.h"
#include "Elements/OpalSource.h"
#include "Elements/OpalStripper.h"
#include "Elements/OpalTravelingWave.h"
#include "Elements/OpalVacuum.h"
#include "Elements/OpalVariableRFCavity.h"
#include "Elements/OpalVariableRFCavityFringeField.h"
#include "Elements/OpalVerticalFFAMagnet.h"
#include "Elements/OpalVKicker.h"

#ifdef ENABLE_OPAL_FEL
#include "Elements/OpalUndulator.h"
#endif

// Structure-related commands.
#include "Lines/Line.h"
#include "Lines/Sequence.h"

// Optimize command
#include "Optimize/OptimizeCmd.h"
#include "Optimize/DVar.h"
#include "Optimize/Objective.h"
#include "Optimize/Constraint.h"

// Sample command
#include "Sample/SampleCmd.h"
#include "Sample/OpalSample.h"

#include "changes.h"

// Modify these methods to add new commands.
// ------------------------------------------------------------------------

namespace {

    void makeActions() {
        OpalData *opal = OpalData::getInstance();
        opal->create(new Call());
        opal->create(new DumpFields());
        opal->create(new DumpEMFields());
        opal->create(new Echo());
        opal->create(new Help());
        opal->create(new List());
        opal->create(new Option());
        opal->create(new OptimizeCmd());
        opal->create(new SampleCmd());
        opal->create(new Select());
        opal->create(new Stop());
        opal->create(new Quit());
        opal->create(new PSystem());
        opal->create(new System());
        opal->create(new Title());
        opal->create(new TrackCmd());
        opal->create(new Value());
    }

    void makeDefinitions() {
        OpalData *opal = OpalData::getInstance();
        // Must create the value definitions first.
        opal->create(new BoolConstant());
        opal->create(new RealConstant());
        opal->create(new RealVariable());
        opal->create(new RealVector());
        opal->create(new StringConstant());

        opal->create(new Beam());
        opal->create(new FieldSolver());
        opal->create(new BoundaryGeometry());
        opal->create(new OpalWake());
        opal->create(new ParticleMatterInteraction());

        opal->create(new OpalFilter());
        opal->create(new OpalTrimCoil());

        opal->create(new Distribution());

        opal->create(new MacroCmd());

        opal->create(new DVar());
        opal->create(new Objective());
        opal->create(new Constraint());

        opal->create(new OpalSample());
    }

    void makeElements() {
        OpalData *opal = OpalData::getInstance();
        opal->create(new OpalAsymmetricEnge());
        opal->create(new OpalCavity());
        opal->create(new OpalCCollimator());
        opal->create(new OpalCyclotron());
        opal->create(new OpalDegrader());
        opal->create(new OpalDrift());
        opal->create(new OpalECollimator());
        opal->create(new OpalFlexibleCollimator());
        opal->create(new OpalHKicker());
        opal->create(new OpalKicker());
        opal->create(new OpalMarker());
        opal->create(new OpalMonitor());
        opal->create(new OpalMultipole());
        opal->create(new OpalMultipoleT());
        opal->create(new OpalMultipoleTStraight());
        opal->create(new OpalMultipoleTCurvedConstRadius());
        opal->create(new OpalMultipoleTCurvedVarRadius());
        opal->create(new OpalOctupole());
        opal->create(new OpalOffset::OpalLocalCartesianOffset());
//        opal->create(new OpalOffset::OpalLocalCylindricalOffset());
//        opal->create(new OpalOffset::OpalGlobalCartesianOffset());
//        opal->create(new OpalOffset::OpalGlobalCylindricalOffset());
        opal->create(new OpalPepperPot());
        opal->create(new OpalPolynomialTimeDependence());
        opal->create(new OpalProbe());
        opal->create(new OpalQuadrupole());
        opal->create(new OpalRBend());
        opal->create(new OpalRBend3D());
        opal->create(new OpalRCollimator());
        opal->create(new OpalRingDefinition());
        opal->create(new OpalSBend());
        opal->create(new OpalSBend3D());
        opal->create(new OpalScalingFFAMagnet());
        opal->create(new OpalSeptum());
        opal->create(new OpalSextupole());
        opal->create(new OpalSlit());
        opal->create(new OpalSolenoid());
        opal->create(new OpalSource());
        opal->create(new OpalStripper());
        opal->create(new OpalTravelingWave());
#ifdef ENABLE_OPAL_FEL        
        opal->create(new OpalUndulator());
#endif
        opal->create(new OpalVacuum());
        opal->create(new OpalVariableRFCavity());
        opal->create(new OpalVariableRFCavityFringeField());
        opal->create(new OpalVerticalFFAMagnet());
        opal->create(new OpalVKicker());

        opal->create(new Line());
        opal->create(new Sequence());
    }
};

namespace Configure {
    void configure() {
        makeDefinitions();
        makeElements();
        makeActions();
        Versions::fillChanges();
    }
};
