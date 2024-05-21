//
// Class SDDSWriter
//   This class is the base class for all SDDS writers.
//
// Copyright (c) 2019, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
//                     Christof Metzger-Kraus, Open Sourcerer
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
#include "Structure/SDDSWriter.h"

#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "OPALconfig.h"
#include "Util/SDDSParser.h"
#include "Utilities/Util.h"

#include "Utility/IpplInfo.h"

#include <queue>


SDDSWriter::SDDSWriter(const std::string& fname, bool restart)
    : fname_m(fname)
    , mode_m(std::ios::out)
    , indent_m("        ")
{
    namespace fs = boost::filesystem;

    if (fs::exists(fname_m) && restart) {
        mode_m = std::ios::app;
        INFOMSG("* Appending data to existing data file: '" << fname_m << "'" << endl);
    } else {
        INFOMSG("* Creating new file for data: '" << fname_m << "'" << endl);
    }
}


void SDDSWriter::rewindLines(size_t numberOfLines) {
    if (numberOfLines == 0 || Ippl::myNode() != 0) {
        return;
    }

    std::string line;
    std::queue<std::string> allLines;
    std::fstream fs;

    fs.open (fname_m.c_str(), std::fstream::in);

    if (!fs.is_open()) return;

    while (getline(fs, line)) {
        allLines.push(line);
    }
    fs.close();


    fs.open (fname_m.c_str(), std::fstream::out);

    if (!fs.is_open()) return;

    while (allLines.size() > numberOfLines) {
        fs << allLines.front() << "\n";
        allLines.pop();
    }
    fs.close();
}


void SDDSWriter::replaceVersionString() {

    if (Ippl::myNode() != 0)
        return;

    std::string versionFile;
    SDDS::SDDSParser parser(fname_m);
    parser.run();
    parser.getParameterValue("revision", versionFile);

    std::string line;
    std::queue<std::string> allLines;
    std::fstream fs;

    fs.open (fname_m.c_str(), std::fstream::in);

    if (!fs.is_open()) return;

    while (getline(fs, line)) {
        allLines.push(line);
    }
    fs.close();


    fs.open (fname_m.c_str(), std::fstream::out);

    if (!fs.is_open()) return;

    while (!allLines.empty()) {
        line = allLines.front();

        if (line != versionFile) {
            fs << line << "\n";
        } else {
            fs << OPAL_PROJECT_NAME << " "
               << OPAL_PROJECT_VERSION << " git rev. #"
               << Util::getGitRevision() << "\n";
        }

        allLines.pop();
    }

    fs.close();
}


double SDDSWriter::getLastValue(const std::string& column) {
    SDDS::SDDSParser parser(fname_m);
    parser.run();
    double val = 0.0;
    parser.getValue(-1, column, val);
    return val;
}


void SDDSWriter::open() {
    if ( Ippl::myNode() != 0 || os_m.is_open() )
        return;

    os_m.open(fname_m.c_str(), mode_m);
    os_m.precision(precision_m);
    os_m.setf(std::ios::scientific, std::ios::floatfield);
}


void SDDSWriter::close() {
    if ( Ippl::myNode() == 0 && os_m.is_open() ) {
        os_m.close();
    }
}


void SDDSWriter::writeHeader() {
    if ( Ippl::myNode() != 0 || mode_m == std::ios::app )
        return;

    this->writeDescription();

    this->writeParameters();

    this->writeColumns();

    this->writeInfo();

    mode_m = std::ios::app;
}


void SDDSWriter::writeDescription() {
    os_m << "SDDS1" << std::endl;
    os_m << "&description\n"
         << indent_m << "text=\"" << desc_m.first << "\",\n"
         << indent_m << "contents=\"" << desc_m.second << "\"\n"
         << "&end\n";
}


void SDDSWriter::writeParameters() {
    while ( !params_m.empty() ) {
        param_t param = params_m.front();

        os_m << "&parameter\n"
             << indent_m << "name=" << std::get<0>(param) << ",\n"
             << indent_m << "type=" << std::get<1>(param) << ",\n"
             << indent_m << "description=\"" << std::get<2>(param) << "\"\n"
             << "&end\n";

        params_m.pop();
    }
}


void SDDSWriter::writeColumns() {
    columns_m.writeHeader(os_m, indent_m);
}


void SDDSWriter::writeInfo() {
    os_m << "&data\n"
         << indent_m << "mode=" << info_m.first << ",\n"
         << indent_m << "no_row_counts=" << info_m.second << "\n"
         << "&end";

    while (!paramValues_m.empty()) {
        os_m << "\n" << paramValues_m.front();

        paramValues_m.pop();
    }

    os_m << std::endl;
}

void SDDSWriter::addDefaultParameters() {
    std::stringstream revision;
    revision << OPAL_PROJECT_NAME << " "
             << OPAL_PROJECT_VERSION << " "
             << "git rev. #" << Util::getGitRevision();

    std::string flavor;
    if (OpalData::getInstance()->isInOPALTMode()) {
        flavor = "opal-t";
    } else if (OpalData::getInstance()->isInOPALCyclMode()) {
        flavor = "opal-cycl";
    } else {
        flavor = "opal-map";
    }

    addParameter("processors", "long", "Number of Cores used", Ippl::getNodes());

    addParameter("revision", "string", "git revision of opal", revision.str());

    addParameter("flavor", "string", "OPAL flavor that wrote file", flavor);
}