//
// File Types
//
// Copyright (c) 2010 - 2013, Yves Ineichen, ETH Zürich
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Toward massively parallel multi-objective optimization with application to
// particle accelerators" (https://doi.org/10.3929/ethz-a-009792359)
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
#ifndef __TYPES_H__
#define __TYPES_H__

#include <map>
#include <utility>
#include <vector>

#include "boost/tuple/tuple.hpp"

#include "Expression/Expression.h"

//FIXME: add namespaces

/// roles a processor can attain
enum Role_t {WORKER, OPTIMIZER, POLLER, UNASSIGNED};


// Variables

// information the optimizer can request: either wants a derivative or a
// function evaluation.
enum InfoType_t {EVALUATE, DERIVATE};
enum Position_t {TIME, POSITION};

typedef std::pair<std::string, double> namedVariable_t;
typedef std::map<std::string, double>  namedVariableCollection_t;

typedef namedVariableCollection_t Param_t;
typedef namedVariableCollection_t variableDictionary_t;

/// type of an expression value is either a single double in case of
/// objectives and for constraints we include the value of LHS and RHS.
//typedef boost::variant< double, boost::tuple<double, double, double> >
    //reqVarValue_t;

//FIXME: do we need InfoType_t ?
/** a requested variable has the following form:
 *
 *  - type of the request (derivation, evaluation, ...)
 *  - vector value (objectives just a single value; constraints return three
 *    values: actual value, lhs value, rhs value)
 *  - boolean denoting the status of the evaluation and if the returned result
 *    is valid
 */
typedef struct reqVarInfo {
    InfoType_t          type;
    std::vector<double> value;
    bool                is_valid;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/) {
        ar & type;
        ar & value;
        ar & is_valid;
    }
} reqVarInfo_t;

typedef std::pair<std::string, reqVarInfo_t> namedReqVar_t;
typedef std::map<std::string, reqVarInfo_t> reqVarContainer_t;



/// type of design variables
typedef boost::tuple<std::string, double, double> DVar_t;
enum DVar_tIdx {
    VAR_NAME,
    LOWER_BOUND,
    UPPER_BOUND
};

typedef std::pair<std::string, DVar_t> namedDVar_t;
typedef std::map<std::string, DVar_t> DVarContainer_t;

#endif
