//
// Class PeakFinder
//   Find peaks of radial profile.
//   It computes a histogram based on the radial distribution of the particle
//   bunch. After that all peaks of the histogram are searched.
//   The radii are written in ASCII format to a file.
//   This class is used for the cyclotron probe element.
//
// Copyright (c) 2017 - 2021, Matthias Frey, Jochem Snuverink, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef PEAKFINDER_H
#define PEAKFINDER_H

#include "Algorithms/Vektor.h"

#include <fstream>
#include <string>
#include <vector>
#include <list>

class PeakFinder {

public:
    using container_t = std::vector<double>;

public:

    PeakFinder() = delete;

    PeakFinder(std::string elem, double min, double max, double binwidth, bool singlemode);

    /*!
     * Append the particle coordinates to the container
     * @param R is a particle position (x, y, z)
     */
    void addParticle(const Vector_t& R);

    /*!
     * Evaluate the centroid of a turn.
     */
    void evaluate(const int& turn);

    void save();

private:

    // compute global histogram, involves some inter-node communication
    void createHistogram_m();

    /***************
     * Output file *
     ***************/
    /// Open output file
    void open_m();
    /// Open output file in append mode
    void append_m();
    /// Close output file
    void close_m();
    /// Write to output file
    void saveASCII_m();

    void computeCentroid_m();

private:
    container_t radius_m;
    /// global histogram values
    container_t globHist_m;

    /// filename with extension (.peaks)
    std::string fileName_m;

    /// histogram filename with extension (.hist)
    std::string hist_m;

    /// used to write out the data
    std::ofstream os_m;

    /// used to write out the histrogram
    std::ofstream hos_m;

    /// Element/probe name, for name output file
    std::string outputName_m;

    // Histogram details
    /// Number of bins
    unsigned int nBins_m;
    /// Bin width in mm
    double binWidth_m;
    ///@{ histogram size
    double min_m, max_m;

    int turn_m;
    double peakRadius_m;
    int registered_m;
    std::list<double> peaks_m;
    bool singlemode_m;
    bool first_m;
    bool finished_m;
};

#endif
