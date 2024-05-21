//
// Struct SumErrSqRadialPeak
//   A simple expression computing the sum of all peak errors (given as
//   first and second argument) for a range of peaks (third argument and fourth argument)
//   according to
//
//   \f[
//     result = \frac{1}{n} * \sqrt{\sum_{i=start}^end (measurement_i - value_i)^2}
//   \f]
//
// Copyright (c) 2017 - 2018, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
// and the paper
// "Matching of turn pattern measurements for cyclotrons using multiobjective optimization"
// (https://doi.org/10.1103/PhysRevAccelBeams.22.064602)
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
#ifndef __SUMERRSQRADIALPEAK_H__
#define __SUMERRSQRADIALPEAK_H__

#include <string>

#include "boost/type_traits/remove_cv.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"
#include "boost/smart_ptr.hpp"

#include "Util/Types.h"
#include "Util/PeakReader.h"
#include "Expression/Parser/function.hpp"

struct SumErrSqRadialPeak {

    static const std::string name;

    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 4) {
            throw OptPilotException("SumErrSqRadialPeak::operator()",
                                    "sumErrSqRadialPeak expects 4 arguments, " + std::to_string(args.size()) + " given");
        }

        meas_filename_ = boost::get<std::string>(args[0]);
        sim_filename_  = boost::get<std::string>(args[1]);
        begin_         = boost::get<double>(args[2]);
        end_           = boost::get<double>(args[3]);

        bool is_valid = true;

        boost::scoped_ptr<PeakReader> meas_peaks(new PeakReader(meas_filename_));
        boost::scoped_ptr<PeakReader> sim_peaks(new PeakReader(sim_filename_));
        try {
            sim_peaks->parseFile();
            meas_peaks->parseFile();

            if ( end_ < begin_ || end_ < 0 || begin_ < 0 )
                throw OptPilotException("SumErrSqRadialPeak::operator()",
                                        "Error check turn number range");

        } catch (OptPilotException &ex) {
            std::cout << "Caught exception: " << ex.what() << std::endl;
            is_valid = false;
        }

        double sum = 0;
        int nPeaks = end_ - begin_ + 1;

        for (int turn = begin_; turn < end_ + 1; ++turn) {
            double sim_value = 0.0, meas_value = 0.0;
            try {
                sim_peaks->getPeak(turn, sim_value);
                meas_peaks->getPeak(turn, meas_value);
            } catch(OptPilotException &e) {
                std::cout << "Exception while getting value "
                          << "from peak file: " << e.what()
                          << std::endl;
                is_valid = false;
            }
            double val = meas_value - sim_value;
            sum += val * val;
        }

        return boost::make_tuple(std::sqrt(sum) / (double)nPeaks, is_valid);
    }

private:
    std::string meas_filename_;
    std::string sim_filename_;
    int begin_;
    int end_;

    // define a mapping to arguments in argument vector
    boost::tuple<std::string, std::string, int, int> argument_types;
    // :FIXME: remove unused enum
#if 0
    enum {
          meas_filename
        , sim_filename
        , begin
        , end
    } argument_type_id;
#endif
};

#endif
