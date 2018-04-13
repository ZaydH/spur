/**
 * stack.cpp
 *
 * Purpose: Defines methods for the StackLevel() class.
 *
 * @author Zayd Hammoudeh <zayd@ucsc.edu>
 * @version 0.00.00
 *
 * Copyright (C) 2018 Zayd Hammoudeh.
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms of the MIT license.  See the
 * LICENSE file for details.
 *
 * Original Author: Marc Thurley.
 */

#include <string>
#include "stack.h"

SolverConfiguration* StackLevel::config_ = nullptr;
DataAndStatistics* StackLevel::statistics_ = nullptr;
VariableIndex StackLevel::component_split_depth_ = 0;


void StackLevel::includeSolutionVar(const mpz_class &solutions, mpz_class solution_var[],
                                    const std::string &name_var) {
  if (branch_found_unsat_[active_branch_]) {
    assert(solution_var[active_branch_] == 0);
    return;
  }
  if (solutions == 0)
    branch_found_unsat_[active_branch_] = true;
  if (solution_var[active_branch_] == 0) {
    solution_var[active_branch_] = solutions;
    if (config_->verbose && solutions > 0)
      std::cout << "\t" << name_var << ": Solution count MPZ set to " << solutions << std::endl;
  } else {
    solution_var[active_branch_] *= solutions;
    if (config_->verbose)
      std::cout << "\t" << name_var << ": Solution count MPZ multiplied by "
                << solutions << " for a product of " << solution_var[active_branch_] << std::endl;
  }
}

SampleAssignment StackLevel::random_cache_sample() const {
  // Only get the sample after both branches explored.
  assert(active_branch_);
  // Must have at least one valid sample to have a random sample.
  assert(branch_model_count_[0] > 0 || branch_model_count_[1] > 0);

  // If either side is UNSAT, make the easy choice
  if (branch_model_count_[0] == 0)
    return cache_sample_[1];
  if (branch_model_count_[1] == 0)
    return cache_sample_[0];

  // Randomly select either the positive or negative side.
  mpz_class rand_mpz, total_count = branch_model_count_[0] + branch_model_count_[1];
//    SamplesManager::get_rand_mpz(total_count, uniform_mpz);
  Random::Mpz::uniform(total_count, rand_mpz);
  if (rand_mpz < branch_model_count_[0])
    return cache_sample_[0];
  else
    return cache_sample_[1];
}

void StackLevel::includeSolution(unsigned solutions) {
  if (branch_found_unsat_[active_branch_]) {
    assert(branch_model_count_[active_branch_] == 0);
    return;
  }
  if (solutions == 0)
    branch_found_unsat_[active_branch_] = true;
  if (branch_model_count_[active_branch_] == 0) {
    branch_model_count_[active_branch_] = solutions;
    if (config_->verbose) {
      std::cout << "\tSolution count set to " << solutions << std::endl;
    }
  } else {
    branch_model_count_[active_branch_] *= solutions;
    if (config_->verbose) {
      std::cout << "\tSolution count multiplied by " << solutions << " for a product of "
                << branch_model_count_[active_branch_] << std::endl;
    }
  }
}
