//
// Class Sampler
//   This class creates, dispatches and dumps new individuals.
//
// Copyright (c) 2018, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
//                     Yves Ineichen, ETH Zürich
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
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
#ifndef __OPAL_SAMPLER_H__
#define __OPAL_SAMPLER_H__

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <queue>
#include <utility>
#include <fstream>
#include <list>
#include <memory>
#include <queue>

#include "Comm/types.h"
#include "Util/Types.h"
#include "Util/CmdArguments.h"

#include "Optimizer/Optimizer.h"
#include "Sample/SampleIndividual.h"
#include "Sample/SamplingMethod.h"

#include <boost/smart_ptr.hpp>

class Sampler : public Optimizer {

public:

    /** This constructor should never be called.
     *  It's provided due to the inheritance of SamplePilot from
     *  Pilot
     *
     */
    Sampler(Expressions::Named_t objectives,
            Expressions::Named_t constraints,
            DVarContainer_t dvars,
            size_t dim, Comm::Bundle_t comms,
            CmdArguments_t args,
            std::vector<double> hypervolRef,
            int nrWorkerGroups);


    /**
     *  Retrieves all (for the sampler) relevant arguments specified on the
     *  command line, initializes the variator and sets up statistics and
     *  debug traces.
     *
     *  @param[in] sampleMethods per design variable (dvar)
     *  @param[in] dvars of sampling
     *  @param[in] comms available to the sampler
     *  @param[in] args the user passed on the command line
     */
    Sampler(const std::map< std::string,
                            std::shared_ptr<SamplingMethod>
                >& sampleMethods,
            Expressions::Named_t objectives,
            DVarContainer_t dvars,
            Comm::Bundle_t comms,
            CmdArguments_t args);

    /// Initialization and start algorithm
    virtual void initialize();

    /// type used in solution state exchange with other optimizers
    typedef std::vector< SampleIndividual > SolutionState_t;

protected:

    // implementing poller hooks
    bool onMessage(MPI_Status status, size_t length);
    void postPoll();

    void setupPoll() {}
    void prePoll() {}
    void onStop() {}

    // helper sending evaluation requests to the pilot
    void dispatch_forward_solves();

private:

    std::map<std::string,
             std::shared_ptr<SamplingMethod>
        > sampleMethods_m;

    // global index (for job id)
    int gid;

    int my_local_pid_;

    typedef SampleIndividual  Individual_t;

    /// communicator bundle for the optimizer
    Comm::Bundle_t comms_;

    /// mapping from unique job ID to individual
    std::map<size_t, boost::shared_ptr<Individual_t> > jobmapping_m;

    std::queue<boost::shared_ptr<Individual_t> > individuals_m;

    /// bounds on each specified gene
    bounds_t dVarBounds_m;

    /// design variables
    DVarContainer_t dvars_m;

    /// objectives
    Expressions::Named_t objectives_m;

    int nSamples_m;


    /// command line arguments specified by the user
    CmdArguments_t args_;

    /// current generation
    int act_sample_m;

    int done_sample_m;

    enum State {
        SUBMIT,
        STOP,
        TERMINATE
    };

    State curState_m;

    /// Dumps id, design variables and bound
    std::size_t jsonDumpFreq_m;
    std::string jsonFname_m;
    void writeJsonHeader();
    std::list<Individual_t> individualsToDump_m;

    void dumpIndividualsToJSON();
    void addIndividualToJSON(const boost::shared_ptr<Individual_t>& ind);

    void runStateMachine();

    void createNewIndividual();
};

#endif