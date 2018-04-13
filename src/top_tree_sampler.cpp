/**
 * top_tree_sampler.h
 *
 * Purpose: Defines the static attributes of the TopTreeSampler() class.
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

#include <gmp.h>

#include <fstream>
#include <string>

#include "primitive_types.h"
#include "top_tree_sampler.h"
#include "solver_config.h"
#include "statistics.h"

const VariableIndex TopTreeSampler::TOP_COMPONENT_SPLIT_DEPTH_ = 0;

SolverConfiguration * TopTreeSampler::config_;
SolverConfiguration * TopTreeSampler::root_config_;
DataAndStatistics * TopTreeSampler::statistics_;
DataAndStatistics * TopTreeSampler::root_statistics_;

std::ofstream TopTreeSampler::all_top_tree_fout_;
std::string TopTreeSampler::all_top_tree_filename_;

mpz_class TopTreeSampler::total_models_ = 0;

std::vector<VariableIndex> TopTreeSampler::literal_stack_start_loc_;
std::vector<LiteralID> TopTreeSampler::debug_lit_stack_;

TopTreeSampler::RemainingFormulas TopTreeSampler::remaining_formulas;
TopTreeSampler::BaseAssignments TopTreeSampler::output_partials;

SampleSize TopTreeSampler::tot_num_rough_samples_ = 0;
float TopTreeSampler::sample_count_scalar_ = -1;
