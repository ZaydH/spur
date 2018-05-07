/**
 * solver.h
 *
 * Purpose: Primary function is to define the Solver() class which is what runs for the entire
 * program.  It also
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

#ifndef SOLVER_H_
#define SOLVER_H_

#include <sys/time.h>

#include <algorithm>
#include <utility>
#include <vector>
#include <string>

#include "statistics.h"
#include "instance.h"
#include "component_management.h"

#include "solver_config.h"
#include "model_sampler.h"
#include "cached_assignment.h"
//#include "top_tree_sampler.h"


class StopWatch {
  /**
   * Used to track time events in the execution of the solver.
   *
   * Built using the internal C++ timer (i.e., timeval()).  Resolution down to microseconds
   * via the "tv_usec" property.
   *
   * @see timeval
   */
 public:
  StopWatch();
  /**
   * Checks whether the time bound has been exceeded.
   *
   * @return true if the time bound has been exceeded.
   */
  bool timeBoundBroken() {
    timeval actual_time = (struct timeval) {0};
    gettimeofday(&actual_time, nullptr);
    return actual_time.tv_sec - start_time_.tv_sec > time_bound_;
  }

  bool start() {
    auto ret = static_cast<bool>(gettimeofday(&last_interval_start_, nullptr));
    start_time_ = stop_time_ = last_interval_start_;
    return !ret;
  }
  /**
   * Records and saves the stop time independent of the time zone.
   *
   * @return 0 for success, or -1 for failure (in which case errno is set appropriately).
   */
  bool stop() {
    return gettimeofday(&stop_time_, nullptr) == 0;
  }

  double getElapsedSeconds() {
    timeval r = getElapsedTime();
    return r.tv_sec + static_cast<double>(r.tv_usec) / 1000000;
  }

  bool interval_tick() {
    timeval actual_time = (struct timeval) {0};
    gettimeofday(&actual_time, nullptr);
    if (actual_time.tv_sec - last_interval_start_.tv_sec
        > interval_length_.tv_sec) {
      gettimeofday(&last_interval_start_, nullptr);
      return true;
    }
    return false;
  }

  /**
   * Updates the stopwatch's time bound.
   *
   * @param seconds New time bound.
   */
  void setTimeBound(uint64_t seconds) {
    time_bound_ = seconds;
  }
  // ZH Appears to neither be implemented or used.
  //long int getTimeBound();

 private:
  timeval start_time_;
  timeval stop_time_;

  uint64_t time_bound_ = UINT64_MAX;

  /**
   * Interval used to separate events.  It is set by default to 60s.
   *
   * For instance, this is used to set how often
   */
  timeval interval_length_;
  timeval last_interval_start_;

  /**
   * if we have started and then stopped the watch, this returns
   * the elapsed time.  Otherwise, time elapsed from start_time_
   * till now is returned
   */
  timeval getElapsedTime();
};

/**
 * Main object managing the execution of both the sampler and
 * sharpSAT.
 */
class Solver: public Instance {
 public:
  /**
   * Constructs a Solver object from command line input arguments.
   *
   * @param argc Number of arguments in the  argv array
   * @param argv Set of input arguments
   */
  Solver(int argc, char *argv[]);

  Solver(SolverConfiguration &config, DataAndStatistics &statistics,
         SampleSize num_samples = 0, bool two_pass = true)
      : config_(config), final_samples_(num_samples, config_) {
    LinkConfigAndStatistics();

    statistics_.input_file_ = statistics.input_file_;

    // ToDo modify to support caching multiple samples
    //config_.num_samples_to_cache_ = num_samples;
    config_.num_samples_to_cache_ = 1;
    sampler_stopwatch_.setTimeBound(config_.time_bound_seconds);  // Initialize the time bound.
    config_.samples_output_file = "";

    config_.disable_samples_write_ = true;
    config_.num_samples_ = num_samples;
    time(&statistics_.start_time_);
    config_.perform_two_pass_sampling_ = two_pass;
  }
  /**
   * Copy constructor
   *
   * @param other Another solver used as the basic of this solver.
   */
  Solver(const Solver &other)
      : config_(other.config_), final_samples_(other.final_samples_.num_samples(), config_) {
    // Bring over the configuration and statistics
    statistics_ = other.statistics_;
    LinkConfigAndStatistics();
    config_.num_samples_to_cache_ = 1;
    config_.disable_samples_write_ = true;

    // Data structures used in the createFromFile Function
    literal_pool_ = other.literal_pool_;
    variables_ = other.variables_;
    literal_values_ = other.literal_values_;
    unit_clauses_ = other.unit_clauses_;
    unused_vars_ = other.unused_vars_;

    independent_support = other.independent_support;
    has_independent_support_ = other.has_independent_support_;
    conflict_clauses_ = other.conflict_clauses_;
    occurrence_lists_ = other.occurrence_lists_;
    literals_ = other.literals_;

    // Initialize other variables as needed
    statistics_.start_time_ = other.statistics_.start_time_;  // Needed to ensure timeout
    sampler_stopwatch_.setTimeBound(config_.time_bound_seconds);  // Initialize the time bound.
    original_lit_pool_size_ = other.original_lit_pool_size_;
  }
  /**
   * Performs the main solver execution.
   *
   * Additional features include printing timing and other notes
   * during execution.  Also writes the final results to an output
   *
   * This solve method can override the stored input file with a
   * custom specified one.
   *
   * @param partial_assn Initial partial assignment.
   * @param results_output_file Output file to write the results
   */
  void solve(const PartialAssignment &partial_assn = PartialAssignment());
  /**
   * Selects uniformly at random a set of satisfying models from the satisfying
   * set for the implicit Boolean formula.
   */
  void sample_models();
  /**
   * Accessor for the Solver's configuration.
   *
   * @return Reference to the solver's internal configuration object.
   */
  const SolverConfiguration &config() {
    return config_;
  }
  /**
   * Modify the Solver()'s configuration to no longer perform top tree sampling.
   */
  void DisableTopTreeSampling() {
    assert(config_.perform_top_tree_sampling);
    config_.perform_top_tree_sampling = false;
  }
  /**
   * Accessor for the Solver()'s statistics information.
   *
   * @return Solver()'s statistics information.
   */
  const DataAndStatistics &statistics() {
    return statistics_;
  }
  /**
   * Creates a new (i.e., empty) SamplesManager for the solver and pushes it onto the stack.
   */
  void PushNewSamplesManagerOnStack() {
    if (config_.store_sampled_models())
      samples_stack_.emplace_back(SamplesManager(config_.num_samples_, config_));
  }
//  /**
//   * Literal Stack Location Printer
//   *
//   * Debug only tool.  Often, a CNF formula has hundreds of thousands
//   * of variables on the literal stack.  Finding where a specific variable
//   * is on the stack (if at all) can be time consuming.  This function
//   * prints whether the specified variable is on the stack.
//   *
//   * @param var Variable number of interest.
//   */
//  void PrintLiteralStackLocation(VariableIndex var);

  /////////////////////////////////////////////
  //  End Sampler Objects
  /////////////////////////////////////////////

 private:
  /**
   * Stores the solver configuration.  It is filled when the
   * solver is initially created based off command line parameters.
   */
  SolverConfiguration config_;
  /**
   * Decision stack represents binary variable assignments
   * only.  It does not capture BCP transactions.  Each stack
   * level also contains information about the associated
   * component(s).
   */
  DecisionStack stack_;
  /**
   * Stores the current partial assignment.  Literals on the stack
   * can be from branch variable selection in the function
   * @see Solver:decideLiteral or from BCP.
   */
  std::vector<LiteralID> literal_stack_;
  /**
   * Stores the set of components that remain to be analyzed.
   */
  ComponentManager comp_manager_ = ComponentManager(config_, statistics_, literal_values_);
  /**
   * Stopwatch used for measuring elapsed time specifically when
   * doing uniformly random sampling of the solution set. - ZH
   */
  StopWatch sampler_stopwatch_;
  /**
   * Stores the last time the conflict clauses were deleted.
   *
   * This is based off the number of variable decisions (i.e., variable
   * branching) during the flow of the algorithm.  It is NOT based
   * off actual elapsed time.
   */
  uint64_t last_ccl_deletion_time_ = 0;
  /**
   * Stores the last time the conflict clause storage was compacted.
   *
   * This is based off the number of variable decisions (i.e., variable
   * branching) during the flow of the algorithm.  It is NOT based
   * off actual elapsed time.
   */
  uint64_t last_ccl_cleanup_time_ = 0;
  /**
   * Vector of samples used to represent a stack of samples.
   * A stack is needed because of component splits.
   */
  std::vector<SamplesManager> samples_stack_;
  /**
   * Contains the final samples manager.  It is set at the end of
   * phase #1 and phase #2 of the solver.
   */
  SamplesManager final_samples_;
  /**
   * If preprocessing is enabled, this function:
   *   * Performs BCP
   *   * Hardwires and compacts components
   *
   *   Since the formula may be modified by this function,
   *   it also reruns the initialization of the variables
   *   @see Solver:initStack
   *
   * @return False if the formula has been found to be UNSAT.
   */
  bool simplePreProcess();

  bool prepFailedLiteralTest();
  /**
   * Performs the two pass sampling where it builds partial assignments in the first pass and
   * completes any partial formulas in the second pass.
   *
   * @param partial_assn Partial assignment to constrain the input CNF formula.
   */
  void reservoir_sample_models(const PartialAssignment &partial_assn, Solver &temp_solver);
//  /**
//   * MT - we assert that the formula is consistent and has not been
//   * found UNSAT yet hard wires all assertions in the literal stack
//   * into the formula removes all set variables and essentially
//   * reinitializes all further data.
//   */
//  void HardWireAndCompact();
  /**
   * Performs the actual model counting for the Boolean formula.
   *
   * If the execution time bound is encountered, the function will halt
   * and return a TIMEOUT error.
   *
   * @return Solver termination condition.
   */
  SolverExitState countSAT();
  /**
   * Makes a single decision (i.e., selects a branch LITERAL - either positive
   * or negative).  A certain number of decisions also causes a
   * literal score decay.
   *
   * This function is where the num_decisions_ field is incremented.
   */
  void decideLiteral();
  /**
   * Executed inside the @see Solver:countSAT function.
   *
   * It performs Boolean constraint propagation (BCP) to
   * eliminate variables via the unit clause rule.
   *
   * @return false is there is a conflict after the BCP.
   */
  bool bcp();
  /**
   * For each literal in the passed component, this function decays the score
   * of the LITERAL by cutting its score in half.
   *
   * @param comp Formula component
   */
  void decayActivitiesOf(Component & comp) {
    for (auto it = comp.varsBegin(); *it != varsSENTINEL; it++) {
      literal(LiteralID(*it, true)).activity_score_ *=0.5;
      literal(LiteralID(*it, false)).activity_score_ *=0.5;
    }
  }
  /**
   * Performs the failed literal test online to identify better branch
   * variables.
   */
  bool implicitBCP();
  /**
   * this is the actual BCP algorithm
   * starts propagating all literal in literal_stack_
   * beginning at offset start_at_stack_ofs
   */
  bool BCP(VariableIndex start_at_stack_ofs);
  /**
   * Goes to the next active branch of the literal in the stack. This may entail going
   * up many levels (i.e., backjumping in the case of a conflict.  It will also merge the solution
   * counts from the left and right side of a branch.
   *
   * This function also is used to write a components model count into the cache.
   *
   * It also merges the multiplicative solution counts when there is a component decomposition. If
   * an UNSAT was found, it also cleans all polluted cache entries.
   *
   * @return PROCESS_COMPONENT if the item on the top of the stack has another component to process.
   * EXIT if the entire run is complete.  RESOLVED indicates explore second active branch of a
   * variable.
   */
  SolverNextAction backtrack();
  /**
   * Depending on the Solver Configuration @see config_, this function either stores
   * just the model count (or
   */
  void ProcessCacheStore();
  /**
   * This function process a component split backtrack.  There are two primary cases
   * that needed to be handled.
   *
   * <ul>
   *    <li>Stitching - This handles combining the solutions of different
   *        recipes into a single sample.</li>
   *    <li>Merging - When a component split has completed processing,
   *        we backtrack through the tree.  This merges the component split subtree
   *        with the rest of the tree.</li>
   * </ul>
   */
  void ProcessSampleComponentSplitBacktrack();
  /**
   * Conflict Resolver
   *
   * @return If on the current decision level a second
   * branch can be visited, RESOLVED is returned. Otherwise
   * returns BACKTRACK
   */
  SolverNextAction resolveConflict();
  /**
   * Prints to the console and a file the results of the sample
   */
  void PrintFinalSamplerResults();

  /////////////////////////////////////////////
  //  BEGIN small helper functions
  /////////////////////////////////////////////

  /**
   * Calculates the score of the variable.  It includes the component score,
   * and positive negative literal scores.
   *
   * @param v Formula variable
   *
   * @return Variable's heuristic score for branch variable selection.
   */
  float scoreOf(VariableIndex v) {
    float score = comp_manager_.scoreOf(v);
    score += 10.0 * literal(LiteralID(v, true)).activity_score_;
    score += 10.0 * literal(LiteralID(v, false)).activity_score_;
//    score += (10*stack_.get_decision_level()) * literal(LiteralID(v, true)).activity_score_;
//    score += (10*stack_.get_decision_level()) * literal(LiteralID(v, false)).activity_score_;
    return score;
  }
  /**
   * If the literal is UNASSIGNED, it is assigned to true
   * and assigned to the current decision level with the specified
   * antecedent.
   *
   * If the literal is being set as part of preprocessing BCP,
   * it will not have an antecedent.
   *
   * @param lit Literal of interest (either negated or positive)
   * @param ant Antecedent of the literal (if any).
   *
   * @return False if the literal is already assigned and true if
   * it was assigned by this statement.
   */
  bool setLiteralIfFree(LiteralID lit, Antecedent ant = Antecedent(NOT_A_CLAUSE)) {
    if (literal_values_[lit] != X_TRI)
      return false;
    var(lit).decision_level = stack_.get_decision_level();
    var(lit).ante = ant;
    literal_stack_.push_back(lit);
    if (ant.isAClause() && ant.asCl() != NOT_A_CLAUSE)
      getHeaderOf(ant.asCl()).increaseScore();
    literal_values_[lit] = T_TRI;
    literal_values_[lit.neg()] = F_TRI;
    return true;
  }
  /**
   * If the solver is not running in quiet mode, this function will print the current run time
   * stats such as the cache size, cache hit rate, BCP hit rate, etc.
   */
  void printOnlineStats();
//  /**
//   * Literal Stack Contents Checker
//   *
//   * Checks whether the specified variable is in the literal stack.
//   *
//   * @param var Variable of interest
//   * @return true if the specified variable @see var is in the literal stack.
//   */
//  bool IsVarInLiteralStack(VariableIndex var);

  void setConflictState(LiteralID litA, LiteralID litB) {
    violated_clause.clear();
    violated_clause.push_back(litA);
    violated_clause.push_back(litB);
  }

  void setConflictState(ClauseOfs cl_ofs) {
    getHeaderOf(cl_ofs).increaseScore();
    violated_clause.clear();
    for (auto it = beginOf(cl_ofs); *it != SENTINEL_LIT; it++)
      violated_clause.push_back(*it);
  }

  std::vector<LiteralID>::const_iterator TOSLiteralsBegin() {
    return literal_stack_.begin() + stack_.top().literal_stack_ofs();
  }
  /**
   * Clears the stack_ and literal_stack_ objects.
   *
   * Pushes the base object onto the decision stack and switches its branch.  It also reallocates
   * and re-initializes the literal stack.
   *
   * @param resSize Size of the literal stack.
   */
  void initStack(unsigned long resSize = 0) {
    stack_.clear();
    if (resSize != 0)
      stack_.reserve(resSize);
    literal_stack_.clear();
    if (resSize != 0)
      literal_stack_.reserve(resSize);
    // initialize the stack to contain at least level zero
    stack_.push_back(StackLevel(1, 0, 2));
    stack_.back().changeBranch();
    // Reset the samples stack.
    if (config_.store_sampled_models()) {
      samples_stack_.clear();
      samples_stack_.reserve(30);
    }
    set_variable_depth_ = 0;
    statistics_.max_component_split_depth_ = 0;
  }

  const LiteralID &TOS_decLit() const {
    assert(stack_.top_const().literal_stack_ofs() < literal_stack_.size());
    return literal_stack_[stack_.top_const().literal_stack_ofs()];
  }
  /**
   * This function is called in the @see Solver:backtrack and @see Solver:resolveConflict
   * methods.  It clears a set of variable information in the state.  That includes:
   *
   * <ul>
   *   <li>Changing the top of the literal stack.</li>
   *   <li>Reseting the remaining components on the decision stack</li>
   *   <li>Cleans the component manager</li>
   * </ul>
   */
  void reactivateTOS() {
    for (auto it = TOSLiteralsBegin(); it != literal_stack_.end(); it++)
      unSet(*it);
    comp_manager_.cleanRemainingComponentsOf(stack_.top());
    literal_stack_.resize(stack_.top().literal_stack_ofs());
    stack_.top().resetRemainingComps();
  }

  bool fail_test(LiteralID lit) {
    VariableIndex sz = literal_stack_.size();
    // we increase the decLev artificially
    // s.t. after the tentative BCP call, we can learn a conflict clause
    // relative to the assn_ of *jt
    stack_.startFailedLitTest();
    setLiteralIfFree(lit);

    assert(!hasAntecedent(lit));

    bool bSucceeded = BCP(sz);
    if (!bSucceeded)
      recordAllUIPCauses();

    stack_.stopFailedLitTest();

    while (literal_stack_.size() > sz) {
      unSet(literal_stack_.back());
      literal_stack_.pop_back();
    }
    return bSucceeded;
  }
  /////////////////////////////////////////////
  //  BEGIN conflict analysis
  /////////////////////////////////////////////

  // if the state name is CONFLICT,
  // then violated_clause contains the clause determining the conflict;
  std::vector<LiteralID> violated_clause;
  // this is an array of all the clauses found
  // during the most recent conflict analysis
  // it might contain more than 2 clauses
  // but always will have:
  //      uip_clauses_.front() the 1UIP clause found
  //      uip_clauses_.back() the lastUIP clause found
  //  possible clauses in between will be other UIP clauses
  std::vector<std::vector<LiteralID>> uip_clauses_;

  // the assertion level of uip_clauses_.back()
  // or (if the decision variable did not have an antecedent
  // before) then assertionLevel_ == DL;
  int64_t assertion_level_ = 0;

  // build conflict clauses from most recent conflict
  // as stored in state_.violated_clause
  // solver state must be CONFLICT to work;
  // this first method record only the last UIP clause
  // so as to create clause that asserts the current decision
  // literal
  void recordLastUIPCauses();
  void recordAllUIPCauses();

  void minimizeAndStoreUIPClause(LiteralID uipLit,
                                 std::vector<LiteralID> & tmp_clause,
                                 const bool seen[]);

  // Commented out by ZH as not implemented
//  void storeUIPClause(LiteralID uipLit, std::vector<LiteralID> & tmp_clause);

//  int getAssertionLevel() const {
//    return assertion_level_;
//  }

  /////////////////////////////////////////////
  //  END conflict analysis
  /////////////////////////////////////////////

  /**
   * Set variable depth only includes true branch variable settings.  It does
   * not include variables implied via unit propagation.
   */
  VariableIndex set_variable_depth_ = 0;
  /**
   * Solver Initializer and Preprocessor
   *
   * Initializes the countSAT Solver and the sampler by creating
   * internal data structures using the formula in the specified CNF
   * file.
   *
   * It also runs a preprocessor on the file to perform basic BCP, check
   * immediately for UNSAT, etc.
   *
   * @param partial_assn Partial assignment to optionally constrain the solver.
   */
  bool InitializeSolverAndPreprocess(const PartialAssignment &partial_assn);

  void applyPartialAssignment(const PartialAssignment &partial_assn);
  /**
   * sharpSAT Results Reporter
   *
   * Updates the final statistics at the conclusion of sharpSAT.  It also
   * writes these final run statistics to the specified text file.
   */
  void ReportSharpSatResults();
  /**
   * Initial Sampler
   *
   * Performs the initial sampling.  When only a single sample is
   * requested, this completes all the sampling.  In contrast, when
   * multiple samples are requested, this function performs the initial
   * sampling and builds partial sampling based on cache hits.
   */
  void PerformInitialSampling();
  /**
   * Partial Assignment Filler
   *
   * Fill the partial assignments and make them complete assignments.  This will
   * entail storing partial assignments in cache.
   */
  inline void FillPartialAssignments(Solver &temp_solver);
  /**
   * Checks whether the sample stack size is valid at the end of running model counting.
   *
   * @return true if the samples stack at the end of @see Solver:countSAT is valid.
   */
  inline bool IsEndSamplesStackSizeValid() const {
    return samples_stack_.size() == 1
           || (samples_stack_.size() == 2 && samples_stack_[0].model_count() == 0);
  }
//  /**
//   * Checks whether the current point in the solver execution is valid for storing a top
//   * tree sample.
//   *
//   * @return true If there none of the current literal assignments led to a component split,
//   * the current depth is below the maximum and top tree mode is enabled.
//   */
//  inline bool IsValidTopTreeNodeStorePoint() {
//    return config_.perform_top_tree_sampling && set_variable_depth_ < config_.max_top_tree_depth_;
////           && HasNoUpperComponentSplit();
//  }
  /**
   * Link the current solver's configuration and statistics to the static representations
   * in the stack and samples manager.
   */
  inline void LinkConfigAndStatistics() {
    // ToDo If the solver is multithreaded, static for config and statistics will need to change.
    StackLevel::set_solver_config_and_statistics(config_, statistics_);
//    SamplesManager::set_solver_config(config_);
  }
};

#endif // SOLVER_H_
