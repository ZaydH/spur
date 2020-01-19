/**
 * solver_config.h
 *
 * Purpose: Defines the SolverConfiguration() class stores the state of the sampler and counter.
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

#ifndef SOLVER_CONFIG_H_
#define SOLVER_CONFIG_H_

#include <cstdlib>
#include <cassert>
#include <string>

#include "primitive_types.h"
#include "sampler_tools.h"

/**
 * Singleton object used for storing the configuration
 * settings for the current execution.
 */
struct SolverConfiguration {
  // Support for these features is removed in the sampler.
//  /**
//   * Support to modify the default state of this feature is
//   * disabled in the sampler.
//   *
//   * Forces variable backtracking in the case of a conflict
//   * being founded.
//   */
//  bool perform_non_chron_back_track = true;
  /**
    * Support to modify the default state of this feature is
    * disabled in the sampler.
   */
  bool perform_component_caching = true;
  /**
    * Support to modify the default state of this feature is
    * disabled in the sampler.
   */
  bool perform_failed_lit_test = true;
  /**
   * Support to modify the default state of this feature is
   * disabled in the sampler.
   */
  bool perform_pre_processing = true;
  /**
   * Rather than performing standard model counting, sample satisfying models
   * uniformly at random.
   */
  bool perform_random_sampling_ = true;
  /**
   * Performs sampling of the top of the tree.
   */
  bool perform_top_tree_sampling = false;
  /**
   * Disable writing the samples file.
   */
  bool disable_samples_write_ = false;
  /**
   * When random sampling is enabled, this flag indicates whether the models should
   * currently be stored or other processing is occurring.
   *
   * @return True if sample models should be stored.
   */
  bool store_sampled_models() {
    return perform_random_sampling_ && !perform_top_tree_sampling;
  }
  /**
   * Location to write the output sample and sampler information.
   */
  std::string samples_output_file = "";
  /**
   * Debug mode has special settings to make debugging the program easier.
   */
  bool debug_mode = false;
  /**
   * The fixed seed. If -1 then time is used instead.
   */
  int fixed_seed = -1;

  /**
   * Stores the maximum execution time of the sampler.
   *
   * Default is approximately 27 hours.
   */
  uint64_t time_bound_seconds = UINT64_MAX;
  /**
   * Controls whether verbose state tracking is enabled.  When
   * verbose is enabled, a basic trace printing is also enabled.
   *
   * In some versions of the console, this verbose printing
   * may be multicolored for easier visual tracking.
   */
  bool verbose = false;
  /**
   * Prevents printing to the console.
   *
   * Cannot be true if verbose is true.
   */
  bool quiet = false;
  /**
   * Cache Sampling State Accessor
   *
   * Returns whether the cache currently supports sampling.
   *
   * @return true if samples are being stored in the cache.
   */
  bool perform_sample_caching() {
    assert(perform_random_sampling_ || !perform_sample_caching_);
//    return perform_random_sampling_ && perform_sample_caching_;
    return perform_sample_caching_;
  }
  /**
   * Sample Caching Enabler
   *
   * This enables storing of samples in the cache.  It can only be
   * enabled.  It can never be disabled by design since once
   * the caching stage has been enabled, it can never be disabled.
   */
  void EnableSampleCaching() {
    assert(perform_random_sampling_);
    perform_sample_caching_ = true;
  }
  /**
   * Forces two pass sample construction.  This applies only to the
   * case where the number of samples requested is 1.
   */
  bool perform_two_pass_sampling_ = false;
  /**
   * Number of samples that will be cached at once.
   */
  unsigned num_samples_to_cache_ = 1;
  /**
   * DEBUG_ONLY MODE - This is an unsupported feature used for
   * debug and development purposes only. Almost no one should
   * use it beyond the core development team.
   */
  bool skip_partial_assignment_fill = false;
  /**
   * Number of satisfying assignments to randomly sample.  This is a command
   * line parameter.
   */
  SampleSize num_samples_ = 0;
  /**
   * When trimming the top of the tree, this represents the maximum branch variable depth
   * before sample building stops.
   */
  TreeNodeIndex max_top_tree_depth_ = 25;
  /**
   * Maximum number of samples that can come from a single top tree node.
   */
  SampleSize max_top_tree_leaf_sample_count = 50;
  /**
   * Location to write the top tree samples.
   */
  std::string top_tree_samples_output_file_ = "." FILE_PATH_SEPARATOR
                                              "__final_top_tree_samples.txt";

 private:
  /**
   * Stores whether storing of samples in the cache is enabled.
   */
  bool perform_sample_caching_ = false;
};

#endif /* SOLVER_CONFIG_H_ */
