/**
 * statistics.h
 *
 * Purpose: Defines the DataAndStatistics() class.
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

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <gmpxx.h>

#include <string>
#include <cstdint>
#include <vector>
#include <cfloat>

#include "structures.h"
#include "component_types/cacheable_component.h"

#include "primitive_types.h"
#include "solver_config.h"


class DataAndStatistics {
 public:
  DataAndStatistics() {
    time(&start_time_);
    for (unsigned i = 0; i < static_cast<unsigned>(TopTreeNodeType::NUM_TREE_NODE_TYPES); i++) {
      num_models_by_tree_node_type_[i] = 0;
      num_tree_nodes_by_type_[i] = 0;
    }
  }
  /**
   * Path to the file containing the CNF formula used in the current run.
   */
  std::string input_file_;
  //double time_elapsed_ = 0.0;
  time_t start_time_;
  uint64_t maximum_cache_size_bytes_ = 0;

  SolverExitState exit_state_ = SolverExitState::NO_STATE;
  // different variable counts
  // number of variables  and clauses before preprocessing
  VariableIndex num_original_variables_ = 0;
  ClauseIndex num_original_clauses_ = 0;
  ClauseIndex num_original_binary_clauses_ = 0;
  ClauseIndex num_original_unit_clauses_ = 0;

  /**
   * number of variables remaining
   */
  VariableIndex num_variables_ = 0;
  /**
   * Number of variables that actually occurs in clauses
   */
  VariableIndex num_used_variables_ = 0;
  /**
   * Number of variables that do not appear in any clause.
   */
  VariableIndex num_free_variables_ = 0;
  /**
   * Number of variables in the independent support.
   */
  VariableIndex num_indep_support_variables_ = 0;
  /**
   * Number of long clauses (i.e., greater than 3 unique literals) after pre-processing.
   */
  ClauseIndex num_long_clauses_ = 0;
  /**
   * Number of binary clauses after pre-processing.
   */
  ClauseIndex num_binary_clauses_ = 0;

  ClauseIndex num_long_conflict_clauses_ = 0;
  ClauseIndex num_binary_conflict_clauses_ = 0;

  ClauseIndex times_conflict_clauses_cleaned_ = 0;

  ClauseIndex num_unit_clauses_ = 0;
  /**
   * number of all decisions made - MT
   *
   * Number of times a branch variable is set.  It
   * does not include BCP variable setting.
   */
  unsigned long num_decisions_ = 0;
//  /// number of all implications derived
//  unsigned long num_implications_ = 0;
  //
  /**
   * Number of all failed literal detected.
   */
  unsigned long num_failed_literals_detected_ = 0;
  unsigned long num_failed_literal_tests_ = 0;
  /**
   * Number of all conflicts encountered.
   */
  unsigned long num_conflicts_ = 0;

  // number of clauses overall learned
  unsigned long num_clauses_learned_ = 0;
  /**
   * Maximum number of embedded components splits.
   */
  VariableIndex max_component_split_depth_ = 0;
  /**
   * Maximum depth of branch variable setting.  This excludes variables
   * implied via UP or implied BCP.
   */
  VariableIndex max_branch_var_depth_ = 0;
  /**
   * Number of top tree nodes observed.
   */
  TreeNodeIndex num_top_tree_nodes_ = 0;
  /**
 * Number of nodes of each of the top tree node types
 */
  TreeNodeIndex num_tree_nodes_by_type_[static_cast<long>(TopTreeNodeType::NUM_TREE_NODE_TYPES)];
  /**
   * Number of models for each of the top tree node types.
   */
  mpz_class num_models_by_tree_node_type_[static_cast<long>(TopTreeNodeType::NUM_TREE_NODE_TYPES)];
  /**
   * Updates the number of nodes and total model count based on the top tree node type.
   *
   * @param node_type Type of top tree node
   * @param num_models Number of new models for the specified node type.
   */
  void UpdateNodeTypeStatistics(TopTreeNodeType node_type, const mpz_class & num_models) {
    auto node_type_id = static_cast<long>(node_type);
    num_tree_nodes_by_type_[node_type_id]++;
    num_models_by_tree_node_type_[node_type_id] += num_models;
  }
  /**
   * Accesor for the number of top tree nodes of a specified type.
   *
   * @param node_type Type of top tree node.
   *
   * @return Total number of top tree nodes of the specified type.
   */
  TreeNodeIndex num_tree_nodes(TopTreeNodeType node_type) {
    auto node_type_id = static_cast<long>(node_type);
    return num_tree_nodes_by_type_[node_type_id];
  }
  /**
   * Accessor for the total number of models associated with a specified type of top tree node.
   *
   * @param node_type Type of top tree node.
   *
   * @return Total number of models associated with the specific top tree node type.
   */
  mpz_class num_tree_node_models(TopTreeNodeType node_type) {
    auto node_type_id = static_cast<long>(node_type);
    return num_models_by_tree_node_type_[node_type_id];
  }
  /**
   * Calculates the percentage of all models that are of a specific type.
   *
   * @param node_type Type of top tree node (e.g., cylinder, max depth, component split, etc.)
   * @return Percent of all models that fall under the specified node type
   */
  double percent_tree_node_models(TopTreeNodeType node_type) {
    if (final_solution_count_ == 0)
      return -1;
    mpf_class percent_models = num_tree_node_models(node_type);
    percent_models /= static_cast<mpf_class>(final_solution_count_);
    return percent_models.get_d();
  }

  /* cache statistics */
  uint64_t num_cache_hits_ = 0;
  uint64_t num_cache_look_ups_ = 0;
  uint64_t sum_cache_hit_sizes_ = 0;

  uint64_t num_cached_components_ = 0;
  uint64_t sum_size_cached_components_ = 0;

  /**
   * Number of bytes occupied by all cached components.
   */
  uint64_t sum_bytes_cached_components_ = 0;
  // the same number, summing over all components ever stored
  uint64_t overall_bytes_components_stored_ = 0;

  // the above numbers, but without any overhead,
  // counting only the pure data size of the components - without model counts
  uint64_t sum_bytes_pure_cached_component_data_ = 0;
  // the same number, summing over all components ever stored
  uint64_t overall_bytes_pure_stored_component_data_ = 0;


  uint64_t sys_overhead_sum_bytes_cached_components_ = 0;
    // the same number, summing over all components ever stored
  uint64_t sys_overhead_overall_bytes_components_stored_ = 0;

  uint64_t cache_infrastructure_bytes_memory_usage_ = 0;

  uint64_t overall_num_cache_stores_ = 0;

  /*end statistics */

  void reset_statistics() {
    // Reset cache statistics
    num_cache_hits_ = 0;
    num_cache_look_ups_ = 0;

    // Reset literal search operations
    num_failed_literals_detected_ = 0;
    num_failed_literal_tests_ = 0;
    num_conflicts_ = 0;
    num_decisions_ = 0;
//    num_implications_ = 0;
    num_unit_clauses_ = 0;

    // Reset execution parameters
    sampler_time_elapsed_ = 0.0;
    exit_state_ = SolverExitState::NO_STATE;

    // Reset all sampler only variables
    numb_second_pass_vars_.clear();

    sampler_time_elapsed_ = DBL_MAX;
    sampler_pass_1_time_ = DBL_MAX;
    sampler_pass_2_time_ = DBL_MAX;
  }

  /**
   * Checks whether the cache is full (i.e., exceeds 100% capacity).
   *
   * @return true if the cache is full and false otherwise.
   */
  bool cache_full() const {
    return cache_bytes_memory_usage() >= maximum_cache_size_bytes_;
  }

  /**
   *
   * @return Sum of the infrastructure bytes and the sized of the cached components.
   */
  uint64_t cache_bytes_memory_usage() const {
    return cache_infrastructure_bytes_memory_usage_
           + sum_bytes_cached_components_;
  }

  uint64_t overall_cache_bytes_memory_stored() {
      return cache_infrastructure_bytes_memory_usage_
             + overall_bytes_components_stored_;
    }

  void incorporate_cache_store(CacheableComponent &ccomp) {
    sum_bytes_cached_components_ += ccomp.SizeInBytes();
    sum_size_cached_components_ += ccomp.num_variables();
    num_cached_components_++;
    overall_bytes_components_stored_ += ccomp.SizeInBytes();
    overall_num_cache_stores_ += ccomp.num_variables();
    sys_overhead_sum_bytes_cached_components_ += ccomp.sys_overhead_SizeInBytes();
    sys_overhead_overall_bytes_components_stored_ += ccomp.sys_overhead_SizeInBytes();


    sum_bytes_pure_cached_component_data_ += ccomp.data_only_byte_size();
    overall_bytes_pure_stored_component_data_ += ccomp.data_only_byte_size();
  }
  void incorporate_cache_erase(CacheableComponent &ccomp) {
    sum_bytes_cached_components_ -= ccomp.SizeInBytes();
    sum_size_cached_components_ -= ccomp.num_variables();
    num_cached_components_--;
    sum_bytes_pure_cached_component_data_ -= ccomp.data_only_byte_size();

    sys_overhead_sum_bytes_cached_components_ -= ccomp.sys_overhead_SizeInBytes();
  }

  void incorporate_cache_hit(CacheableComponent &ccomp) {
    num_cache_hits_++;
    sum_cache_hit_sizes_ += ccomp.num_variables();
  }
  unsigned long cache_MB_memory_usage() {
    return cache_bytes_memory_usage() / 1000000;
  }
  mpz_class final_solution_count_ = 0;

  double implicitBCP_miss_rate() {
    if (num_failed_literal_tests_ == 0) return 0.0;
    return (num_failed_literal_tests_ - num_failed_literals_detected_)
           / static_cast<double>(num_failed_literal_tests_);
  }
  unsigned long num_clauses() const {
    return num_long_clauses_ + num_binary_clauses_ + num_unit_clauses_;
  }
  unsigned long num_conflict_clauses() {
    return num_long_conflict_clauses_ + num_binary_conflict_clauses_;
  }

  unsigned long clause_deletion_interval() {
    return 10000 + 10 * times_conflict_clauses_cleaned_;
  }
  /**
   * Updates the solver field "final_solution_count_" with the solution count multiplied by
   * 2^{#unused variables}.  Hence:
   *
   * final_solution_count_ = count * 2^(num_variables_ - num_used_variables_)
   *
   * @param count model count.
   */
  void set_final_solution_count(const mpz_class &count) {
    // set final_solution_count_ = count * 2^(num_variables_ - num_used_variables_)
    mpz_mul_2exp(final_solution_count_.get_mpz_t(), count.get_mpz_t(),
                 num_variables_ - num_used_variables_);
  }
  /**
   * Accessor for the final solution count
   *
   * @return Final solution count.
   */
  const mpz_class &final_solution_count() const {
    return final_solution_count_;
  }

//  void incorporateConflictClauseData(const std::vector<LiteralID> &clause) {
//    if (clause.size() == 1)
//      num_unit_clauses_++;
//    else if (clause.size() == 2)
//      num_binary_conflict_clauses_++;
//    num_long_conflict_clauses_++;
//  }
  void incorporateClauseData(const std::vector<LiteralID> &clause) {
    if (clause.size() == 1)
      num_unit_clauses_++;
    else if (clause.size() == 2)
      num_binary_clauses_++;
    else
      num_long_clauses_++;
  }

//  /**
//   * Prints to the console the file solution count.
//   */
//  void print_final_solution_count();
//  /**
//   * Prints statistic information regarding the CNF to a file including, the final
//   * solution count, number of decisions, original formula variable and clause
//   * counts, and elpased time.
//   *
//   * If the formula is not satisfiable, that will be specially printed.
//   *
//   * Format of the output is basic HTML table format.
//   *
//   * @param file_name Path of output file to write
//   */
//  void writeToFile(const std::string & file_name);

  void printShort();
  /**
   * Prints to the console information about the formula including:
   * number of long, binary, & unit clauses, number of used variables,
   * and total number of variables.
   */
  void printShortFormulaInfo() {
    std::cout << "variables (all/used/free): \t"
              << num_variables_ << "/" << num_used_variables_ << "/"
              << num_variables_ - num_used_variables_ << "\n"
              << "independent support size:  \t" << num_indep_support_variables_ << "\n";

    std::cout << "clauses (all/long/binary/unit): "
              << num_clauses() << "/" << num_long_clauses_
              << "/" << num_binary_clauses_ << "/" << num_unit_clauses_ << std::endl;
  }
  /**
   * Number of decisions (i.e., branch variable selections) made.  This
   * is used for determining when to remove components from the cache.
   *
   * @return Number of branch variable selections
   */
  unsigned long getTime() const {
    return num_decisions_;
  }

  long double getAvgComponentSize() {
    return sum_size_cached_components_ / static_cast<long double>(num_cached_components_);
  }
  /**
   * Access the number of components in the cache.
   *
   * @return
   */
  unsigned long cached_component_count() {
    return num_cached_components_;
  }
  /**
   * Number of cache look-ups where the component requested was in the cache.
   *
   * @return Cache look-up count.
   */
  unsigned long cache_hits() {
    return num_cache_hits_;
  }
  /**
   * The rate of cache look-ups where the requested component was not in the cache.
   *
   * @return 0 if a cache look-up never occurred. Otherwise, it is the number of
   * caches misses divided by the total number of cache look-ups.
   */
  double cache_miss_rate() {
    if (num_cache_look_ups_ == 0) return 0.0;
    return (num_cache_look_ups_ - num_cache_hits_) / static_cast<double>(num_cache_look_ups_);
  }
  /**
   * The average size of a reused (i.e., hit) components in the cache.
   *
   * @return 0 if a cache hit has never occurred.  Otherwise, it is the
   * average size of cached components that have ever been HIT.
   */
  long double getAvgCacheHitSize() {
    if (num_cache_hits_ == 0) return 0.0;
    return sum_cache_hit_sizes_ / static_cast<long double>(num_cache_hits_);
  }

  //-----------------------------------------//
  //        Sampler Fields and Methods       //
  //-----------------------------------------//
  /**
   * Amount of time required to complete the full sampling.
   */
  double sampler_time_elapsed_ = DBL_MAX;
  /**
   * Amount of time required to complete phase 1 and build all partial
   * assignments (or just the time required to collect a single
   * sample when two pass is disabled).
   */
  double sampler_pass_1_time_ = DBL_MAX;
  /**
   * Time required to convert all partial assignments into
   * complete assignments.  This represents "pass 2" of the sampler.
   */
  double sampler_pass_2_time_ = DBL_MAX;
  /**
   * When performing two pass testing, this records the number of variables
   * that remain to be set at the start of the second pass for each
   * sample.
   */
  std::vector<VariableIndex> numb_second_pass_vars_;
};

#endif /* STATISTICS_H_ */
