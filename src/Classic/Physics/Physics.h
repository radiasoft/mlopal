//
// Class Physics
//   A namespace defining various mathematical and physical constants.
//
// Copyright (c) 2015-2021, Christof Metzger-Kraus, Helmholtz-Zentrum Berlin
//                          Pedro Calvo, CIEMAT, Spain
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Optimizing the radioisotope production of the novel AMIT
// superconducting weak focusing cyclotron"
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
#ifndef CLASSIC_Physics_HH
#define CLASSIC_Physics_HH


namespace Physics {

    /// The value of \f[ \pi \f]
    constexpr double pi         = 3.14159265358979323846;

    /// The value of \f[2 \pi \f]
    constexpr double two_pi     = 2 * pi;

    /// The value of \f[ \frac{1}{2} \pi \f]
    constexpr double u_two_pi   = 1.0 / two_pi;

    /// The value of \f[ e \f]
    constexpr double e          = 2.7182818284590452354;

    /// The logarithm of $e$ to the base 10
    constexpr double log10e     = 0.43429448190325182765;

    /// The velocity of light in m/s
    constexpr double c          = 299792458.0;

    /// The permeability of vacuum in Vs/Am
    constexpr double mu_0       = 1.25663706212e-06;

    /// The permittivity of vacuum in As/Vm
    constexpr double epsilon_0  = 8.8541878128e-12;

    /// The reduced Planck constant in GeVs
    constexpr double h_bar      = 6.582119569e-25;

    /// The Avogadro's number
    constexpr double Avo        = 6.02214076e23;

    /// Boltzman's constant in eV/K.
    constexpr double kB         = 8.617333262e-5;

    /// Rydberg's energy (Rydberg's constant times hc) in GeV
    constexpr double E_ryd      = 13.605693122994e-9;

    /// Bohr radius in m
    constexpr double a0         = 5.29177210903e-11;

    /// The elementary charge in As
    constexpr double q_e        = 1.602176634e-19;

    /// The fine structure constant, no dimension
    constexpr double alpha      = 7.2973525693e-03;

    /// The atomic mass unit energy equivalent in GeV
    constexpr double amu        = 0.93149410242;

    /// The electron rest mass in GeV
    constexpr double m_e        = 0.51099895000e-03;

    /// The classical electron radius in m
    constexpr double r_e        = 2.8179403262e-15;

    /// The reduced Compton wave length for electrons in m
    constexpr double lamda_e    = 3.8615926796e-13;

    /// The magnetic momentum anomaly for electrons, no dimension
    constexpr double a_e        = 1.15965218128e-03;

    /// The proton rest mass in GeV
    constexpr double m_p        = 0.93827208816e+00;

    /// The classical proton radius in m
    constexpr double r_p        = 1.53469857e-18;

    /// The reduced Compton wave length for protons in m
    constexpr double lamda_p    = 2.10308910336e-16;

    /// The magnetic momentum anomaly for protons, no dimension
    constexpr double a_p        = 1.792847386e+00;

    /// The charge of proton
    constexpr double z_p        = 1;

    /// The carbon-12 (fully-stripped) rest mass in GeV
    constexpr double m_c        = 11.9967074146787 * amu;

    /// The H- rest mass in GeV
    constexpr double m_hm       = 1.00837 * amu;

    /// The uranium-238 (fully-stripped) rest mass in GeV
    constexpr double m_u        = 237.999501 * amu;

    /// The muon rest mass in GeV
    constexpr double m_mu       = 0.1056583755;

    /// The deuteron rest mass in GeV
    constexpr double m_d        = 2.013553212745 * amu;

    /// The xenon-129 (fully-stripped) rest mass in GeV
    constexpr double m_xe       = 128.87494026 * amu;

    /// The alpha particle rest mass in GeV
    constexpr double m_alpha    = 4.001506179127 * amu;

    /// The hydrogen atom rest mass in GeV
    constexpr double m_h        = 1.00782503224 * amu;

    /// The H2+ rest mass in GeV
    constexpr double m_h2p       = 2.01510 * amu;

    /// The H3+ rest mass in GeV
    constexpr double m_h3p       = 3.02293 * amu;

    constexpr double PMASS      = 1.67262192369e-27; // kg

    constexpr double EMASS      = 9.1093837015e-31; // kg

    constexpr double PCHARGE    = 1.602176634e-19; // C

    // alfven current
    constexpr double Ia         = 17.045148e+03;
    // e/m - electron charge to mass quotient
    constexpr double e0m        = 1.75882001076e+11;
    // e/mc
    constexpr double e0mc       = e0m / c;
};

#endif // CLASSIC_Physics_HH
