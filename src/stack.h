/**
 * stack.h
 *
 * Purpose: Defines the StackLevel() class that encapsulates a single level in the decision stack.
 * The decision stack itself is stored in the DecisionStack() class.
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

#ifndef STACK_H_
#define STACK_H_


#include <gmpxx.h>

#include <cassert>
#include <vector>
#include <string>

#include "model_sampler.h"
#include "cached_assignment.h"
#include "solver_config.h"


/**
 * Represents a single level in the decision stack.
 */
class StackLevel {
  /**
   * MT - Active component, once initialized, it should not change
   */
  const VariableIndex super_component_ = 0;
  /**
   * False if exploring the negative site of the branch
   * and true if exploring the positive side.
   */
  bool active_branch_ = false;
  /**
   * Offset in the literal stack where to store set literals
   */
  const VariableIndex literal_stack_ofs_ = 0;
  /**
   * Stores whether the current stack level and active branch
   * corresponds with a component split.
   */
  bool is_component_split_ = false;
  /**
   * This is used for pushing extra components onto the stack.
   */
  bool first_component_ = false;
  /**
   * Stores the sample used to build cached samples.  To ensure
   * randomness for the two active branches, this function keeps
   * two possible cache samples - one for the positive and the other
   * for the negative branch.
   */
  SampleAssignment cache_sample_[2];
  /**
   * Stores the model count for the first and second explored
   * literal values.  These do NOT necessarily correspond to
   * false and true respectively.
   */
  mpz_class branch_model_count_[2] = {0, 0};
  /**
   * Stores whether it was determined if either side of the dependent
   * subtress was unsatisfiable.
   */
  bool branch_found_unsat_[2] = {false, false};
  /**
   * Generic helper function to update the value of one local
   * mpz_class object based off the solution count found.
   *
   * @param solutions Number of new solutions to add
   * @param solution_var Variable storing the solution count
   * @param name_var Name of the variable being updated.  This is for debug printing.
   */
  void includeSolutionVar(const mpz_class &solutions, mpz_class solution_var[],
                          const std::string &name_var);
  /** MT - The start offset in the component stack for the remaining
   * components in this decision level all remaining components can
   * hence be found in the range:
   *
   * [remaining_components_ofs_, "nextLevel".remaining_components_begin_)
   *
   * <b>Note</b>: remaining_component_ofs is inclusive.
   */
  const VariableIndex remaining_components_ofs_ = 0;
  /**
   * MT - Boundary of the stack marking which components still need to
   * be processed all components to be processed can be found in
   *
   * [remaining_components_ofs_, unprocessed_components_end_)
   *
   * <b>Note</b>: remaining_components_ofs_ is inclusive while
   * unprocessed_components_end_ is exclusive.
   *
   * also, all processed, can be found in the range:
   *
   * [unprocessed_components_end_, component_stack.size())
   */
  VariableIndex unprocessed_components_end_ = 0;
  /**
   * Current depth of component splitting.
   */
  static VariableIndex component_split_depth_;

  /*---------------------------------------*
   *       Begin Sampler Variables         *
   *---------------------------------------*/

  /**
   * When deciding through the stack, some variables may become "free" since
   * they no longer appear in any remaining clauses in the residual formula.
   *
   * When a leaf node is reached, the solution count must be multiplied by this
   * multiplier count.
   */
  mpz_class stack_solution_count_multiplier_[2] = {1, 1};
  /**
   * All variables freed at the current level of the stack
   * and in all preceding levels.
   */
  std::vector<std::vector<VariableIndex>> stack_emancipated_vars_;
  /**
   * Store the identification number for the cached entries.
   */
  std::vector<CacheEntryID> cache_comp_ids_[2];
  /**
   * Solver configuration built at the creation of the solver.
   */
  static SolverConfiguration* config_;
  /**
   * Solver statistics used for collecting information about the run.
   */
  static DataAndStatistics* statistics_;
  /**
   * Stored information concerning the cached assignment.
   */
  CachedAssignment cached_assn_;

 public:
  bool hasUnprocessedComponents() const {
    assert(unprocessed_components_end_ >= remaining_components_ofs_);
    return unprocessed_components_end_ > remaining_components_ofs_;
  }
  /**
   * Mark the current unprocessed component as processed.
   *
   * It does verify that at least one unprocessed component exists at this decision level.
   *
   * This function does not return anything.  Rather, the other code operates directly
   * on the component stack and accesses the top of it via the @see StackLevel:top() method.
   */
  void nextUnprocessedComponent() {
    assert(unprocessed_components_end_ > remaining_components_ofs_);
    unprocessed_components_end_--;
  }
  /**
   * Mark all components in this decision level as processed.
   */
  void resetRemainingComps() {
    unprocessed_components_end_ = remaining_components_ofs_;
  }
  const VariableIndex super_component() const {
    return super_component_;
  }
  /**
   * MT - the start offset in the component stack for the remaining components in this
   * decision level all remaining components can hence be found in
   *
   * [remaining_components_ofs_, "nextLevel".remaining_components_begin_)
   *
   * @return Level in the component stack of the first component for this decision level.
   */
  const VariableIndex remaining_components_ofs() const {
    return remaining_components_ofs_;
  }
  /**
   * One greater than the location in the component stack of unprocess components
   * for the current decision level.
   *
   * @return One greater than the location in the component stack for unprocessed
   * components for teh current point in the decision stack.
   */
  const VariableIndex unprocessed_components_end() const {
    return unprocessed_components_end_;
  }
  /**
   * Update the location of the unprocessed component reference for this decision level.
   *
   * @param end New location of the set of unprocessed components for this decision level.
   */
  void set_unprocessed_components_end(VariableIndex end) {
    unprocessed_components_end_ = end;
    assert(remaining_components_ofs_ <= unprocessed_components_end_);
  }

  StackLevel() {
    stack_emancipated_vars_.resize(2);
  }

  StackLevel(VariableIndex super_comp, VariableIndex lit_stack_ofs,
             ClauseOfs comp_stack_ofs) :
      super_component_(super_comp),
      literal_stack_ofs_(lit_stack_ofs),
      remaining_components_ofs_(comp_stack_ofs),
      unprocessed_components_end_(comp_stack_ofs) {
    assert(super_comp < comp_stack_ofs);
    stack_emancipated_vars_.resize(2);
  }
  /**
   * Access for the total number of remaining components across all decision levels not
   * just the current one.
   *
   * @return Total number of remaining components
   */
  VariableIndex currentRemainingComponent() {
    assert(remaining_components_ofs_ <= unprocessed_components_end_ - 1);
    return unprocessed_components_end_ - 1;
  }
  /**
   * Checks if the second (i.e., "true") branch is now active.
   *
   * @return true if the second ("true") branch is active.
   */
  bool isSecondBranch() {
    return active_branch_;
  }
  /**
   * Set the active branch (i.e., the one being searched to
   * the positive branch.
   */
  void changeBranch() {
    active_branch_ = true;
  }
  /**
   * Checks if there is another processible component left.  Lack of processible components
   * could be due to either the branch being unsatisfiable or there being no processible
   * components left
   *
   * @return true if a processible component remains.
   */
  bool anotherCompProcessible() {
    return (!branch_found_unsat()) && hasUnprocessedComponents();
  }

  VariableIndex literal_stack_ofs() const {
    return literal_stack_ofs_;
  }
  /**
   * Updates the solution count.  If no solution have been found,
   * then the increase is additive.  If at least one solution has
   * already been found then the increase is multiplicative.
   *
   * @param solutions Number of solutions to increase by.
   */
  void includeSolution(const mpz_class &solutions) {
    includeSolutionVar(solutions, branch_model_count_, "branch_model_count");
  }
  /**
   * Updates the multiplier for when doing solution sampling.
   *
   * @param solutions Number of solutions to increase by.
   */
  void includeSolutionSampleMultiplier(const mpz_class &solutions) {
    includeSolutionVar(solutions, stack_solution_count_multiplier_,
                       "sampler_solution_multiplier");
  }
  /**
   * Updates the solution count.  If no solution have been found,
   * then the increase is additive.  If at least one solution has
   * already been found then the increase is multiplicative.
   *
   * @param solutions Number of solutions to increase by.
   */
  void includeSolution(unsigned solutions);
  /**
   * When variables become free in a residual formula (i.e., they no longer appear
   * in any clauses), the solution count multiplier is stored.  This needs to be
   * pushed down the stack to ensure that the sample weights are correct.
   *
   * This function pushes the current scalar multiplier down the stack.
   *
   * @param top Current top of the decision stack.
   * @param decision_level Number of branch variables already assigned.  This does not include any
   * implicit partial assignment.
   */
  void configureNewLevel(const StackLevel &top, const DecisionLevel decision_level) {
    // Do not push down any counts at the top of the decision stack when doing reservoir sampling
    // This is required in some cases in top tree sampling since the solution count at each leaf
    // is required while traversing the tree.
    if (decision_level == 0 && !(config_->perform_top_tree_sampling && !top.isComponentSplit()))
      return;
    // Push down the solution multiplier and freed variable list
    // Make sure to only use the active branch.
    for (auto &multiplier : stack_solution_count_multiplier_)
      multiplier = top.getSamplerSolutionMultiplier();

    for (auto &freed_variables : stack_emancipated_vars_)
      freed_variables = top.emancipated_vars();

    for (auto &cache_comp_ids : cache_comp_ids_)
      cache_comp_ids = top.cached_comp_ids();
  }
  /**
   * Freed Variables Adder
   *
   * When performing variable assignments, it is common that variables will become free
   * (i.e., no longer appear in any clauses).  Those variables will not appear in
   * the literal stack since they are not set.  However, they must be assigned as
   * part of the sample generator.  This function stores those freed variables.
   *
   * @param freed_vars Newly freed variables
   */
  void addFreeVariables(const std::vector<VariableIndex> &freed_vars) {
    int idx = (active_branch_) ? 1 : 0;
    stack_emancipated_vars_[idx].insert(stack_emancipated_vars_[idx].end(),
                                       freed_vars.begin(), freed_vars.end());
  }
  /**
   * Store the component IDs for the cached components.
   *
   * @param ids List of cached entry identification numbers
   */
  void addCachedCompIds(const std::vector<CacheEntryID> &ids) {
    int idx = (active_branch_) ? 1 : 0;
//    assert(cached_comp_ids_.empty());
    if (ids.empty())
      return;
    cache_comp_ids_[idx].insert(cache_comp_ids_[idx].end(), ids.begin(), ids.end());
  }
  /**
   * Checks if the currently active branch variable has been
   * determined to be unsatisfiable.
   *
   * @return true if the currently active branch is unsatisfiable.
   */
  inline bool branch_found_unsat() const {
    return branch_found_unsat_[active_branch_];
  }
  /**
   * Sets the currently active branch to being unsatisfiable.
   */
  inline void mark_branch_unsat() {
    branch_found_unsat_[active_branch_] = true;
  }
  /**
   * Total number (i.e., sum of false and true) branches of this
   * literal decision in the subtree.
   *
   * @return Total model count for both the positive and negative subtrees.
   */
  const mpz_class getTotalModelCount() const {
    return branch_model_count_[0] + branch_model_count_[1];
  }
  /**
   * Gets the model count for the active branch.  If the active branch is false,
   * it will return the model count for negated version of the literal.  Otherwise,
   * it return the positive (unnegated) model count.
   *
   * @return Model count for the currently active literal state.
   */
  const mpz_class& getActiveModelCount() const {
    if (!active_branch_)
      return branch_model_count_[0];
    else
      return branch_model_count_[1];
  }
  /**
   * Gets the multiplier count for when running the solution sampler.
   * It is used for the case when variables become free due to the residual
   * formula no longer containing them.
   *
   * @return Solution count multiplier for the sampler.
   */
  const mpz_class& getSamplerSolutionMultiplier() const {
    if (!active_branch_)
      return stack_solution_count_multiplier_[0];
    else
      return stack_solution_count_multiplier_[1];
  }
//  /**
//   * Calculates the total descendent model count for this stack level.  It is equal to:
//   *
//   * FalseModelCount * FalseCountMultiplier + TrueModelCount * TrueCountMultiplier
//   *
//   * @return Total descendent model count.
//   */
//  const mpz_class GetDescendentModelCount() const {
//    mpz_class total_model_count = 0;
//    for (int i = 0; i < 2; i++) {
//      if (stack_solution_count_multiplier_[i] == 1)
//        total_model_count += branch_model_count_[i];
//      else
//        total_model_count += branch_model_count_[i] * stack_solution_count_multiplier_[i];
//    }
//    return total_model_count;
//  }
  /**
   * Checks whether the associated decision level has any descendent models.
   *
   * @return True if both the true and false descendents are UNSAT.
   */
  const bool HasNoDescendentModels() {
    return branch_found_unsat_[0] && branch_found_unsat_[1];
  }
  /**
   * Emancipated Variable List Accessor
   *
   * This function is used to access the freed variable list for a specific
   * variable branch.  This is used for variable assignments.
   *
   * @return List of freed variables for the specified active branch.
   */
  inline const std::vector<VariableIndex> &emancipated_vars() const {
    if (!active_branch_)
      return stack_emancipated_vars_[0];
    else
      return stack_emancipated_vars_[1];
  }
  /**
   * Cached Component ID List Accessor
   *
   * This function accesses the cached entry identification numbers for this stack level.
   *
   * @return Cached component IDs
   */
  inline const std::vector<CacheEntryID> &cached_comp_ids() const {
    if (!active_branch_)
      return cache_comp_ids_[0];
    else
      return cache_comp_ids_[1];
  }
  /**
   * Component Split State Accessor
   *
   * Accessor for whether this stack level corresponds with a component split.
   *
   * It only applies when performing random sampling.  Otherwise, it merely
   * returns false.
   *
   * @return True if component split and random sampling is running.
   */
  inline const bool isComponentSplit() const {
    return is_component_split_ && config_->perform_random_sampling_;
  }
  /**
   * Marks the decision level as a component split.
   *
   * It only applies when performing random sampling.
   */
  inline void setAsComponentSplit() {
    if (!config_->perform_random_sampling_)
      return;

    if (!is_component_split_) {
      component_split_depth_++;
      if (component_split_depth_ > statistics_->max_component_split_depth_)
        statistics_->max_component_split_depth_ = component_split_depth_;
    }
    is_component_split_ = true;
    first_component_ = true;
  }
  /**
   * Unmarks the decision level as a component split.
   *
   * This function only has an effect when performing random sampling.
   */
  inline void unsetAsComponentSplit() {
    if (!config_->perform_random_sampling_)
      return;

    if (is_component_split_)
      component_split_depth_--;
    assert(component_split_depth_>= 0);
    is_component_split_ = false;
    first_component_ = false;
  }
  /**
   * Component Split Depth Accessor
   *
   * Accesses the global component split depth information.
   *
   * @return Current component split depth.
   */
  inline static VariableIndex componentSplitDepth() {
    return component_split_depth_;
  }
  /**
   * Returns whether the first component in a component split has been
   * processed.  This is only relevant when performing random sampling.
   *
   * @return true if the first component has not been processed
   */
  inline const bool isFirstComponent() const {
    assert(is_component_split_);
    return first_component_ && config_->perform_random_sampling_;
  }
  /**
   * This is used for marking when the first component has been completed.
   * We use this for sample stitching.
   */
  inline void markFirstComponentComplete() {
    assert(is_component_split_);
    first_component_ = false;
  }
  /**
   * Configuration Updater
   *
   * Stores the solver configuration.  This is used so that global configuration
   * settings (e.g., verbose) can be used for printing.
   *
   * @param config Solver configuration
   */
  inline static void set_solver_config_and_statistics(SolverConfiguration &config,
                                                      DataAndStatistics &statistics) {
    config_ = &config;
    statistics_ = &statistics;
  }
  /**
   * Cached Assignment Updater
   *
   * Updates the cached assignment associated with this stack level.
   *
   * @param cached_assn Updates the cached assignment information
   * for this StackLevel.
   */
  inline void set_cached_assn(CachedAssignment &cached_assn) {
    assert(config_->perform_random_sampling_);
    cached_assn_ = cached_assn;
  }
  inline const CachedAssignment& cached_assn() const { return cached_assn_; }
  /**
   * Removes all information associated with the cached assignment
   * including the number of components it contains as well
   * as the assignment itself.
   */
  inline void ClearCachedAssn() { cached_assn_.clear(); }
  /**
   * Store the cache sample assignment.  It will be automatically
   * associated with the right active branch.
   *
   * @param cache_sample Sample that will be used for caching this
   * sample.
   */
  inline void set_cache_sample(const SampleAssignment & cache_sample) {
    if (!active_branch_)
      cache_sample_[0] = cache_sample;
    else
      cache_sample_[1] = cache_sample;
  }
  /**
   * Randomly select the cache sample between the two branches.
   *
   * @return A sample that can be stored in the cache for this
   * stack level.
   */
  SampleAssignment random_cache_sample() const;
};

/**
 * Vector based decision stack.  It is also used
 * during implicit BCP (i.e., failed literal testing).
 */
class DecisionStack: public std::vector<StackLevel> {
  VariableIndex failed_literal_test_active = 0;

 public:
  DecisionStack() : std::vector<StackLevel>() {}
  /**
   * Sets a flag to begin implicit BCP
   */
  void startFailedLitTest() {
//    failed_literal_test_active = true;  // Thurley's old code. Relies on casting
    failed_literal_test_active = 1;
  }
  /**
   * Sets a flag to halt implicit BCP.
   */
  void stopFailedLitTest() {
//    failed_literal_test_active = false;  // Thurley's old code. Relies on casting
    failed_literal_test_active = 0;
  }
  // end for implicit BCP
  /**
   * Accessor for the top of the decision stack.
   *
   * @return Top of the decision stack.
   */
  StackLevel &top() {
    return const_cast<StackLevel&>(top_const());  // Strip off the const
  }
  /**
   * Accessor for the top of the decision stack.
   *
   * @return Top of the decision stack.
   */
  const StackLevel &top_const() const {
    assert(!empty());
    return back();
  }
  /**
   * Get the element below the top, making it on deck in baseball
   * parlance (i.e., it would be up next after the top).
   *
   * @return One below the top of the stack.
   */
  StackLevel &on_deck() {
    assert(size() >= 2);
    return (*this)[size() - 2];
  }
  /**
   * Height of the decision stack.  This represents "d" in CDCL.
   *
   * If a failed literal test is underway, a failed literal
   * is pushed onto the stack during testing.
   *
   * @return Height of the stack.
   */
  DecisionLevel get_decision_level() const {
    assert(!empty());
    return size() - 1 + failed_literal_test_active;
  }  // 0 means pre-1st-decision
};

#endif  // STACK_H_
