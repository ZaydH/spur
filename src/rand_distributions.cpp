/**
 * rand_distributions.cpp
 *
 * Purpose: Static class for getting a random variable distributed according to the binomial
 * distribution.
 *
 * @author Zayd Hammoudeh <zayd@ucsc.edu>
 * @version 0.00.00
 *
 * Copyright (C) 2018 Zayd Hammoudeh.
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms of the MIT license.  See the
 * LICENSE file for details.
 */

#include <cstdlib>

#include "rand_distributions.h"

// Satisfy the linker by initializing the static variables
SolverConfiguration * Random::master_config_;
std::random_device Random::rd_;     // only used once to initialise (seed) engine
std::mt19937 Random::rng_;    // random-number engine used (Mersenne-Twister in this case)
std::uniform_int_distribution<int> Random::uni_(INT_MIN, INT_MAX);

bool Random::Mpz::use_approx_binom_ = true;
gmp_randstate_t Random::Mpz::rand_state_;
