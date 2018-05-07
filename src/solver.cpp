/**
 * solver.cpp
 *
 * Purpose: Defines the methods for the Solver() class.
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

#include <deque>
#include <string>
#include <utility>
#include <vector>
#include <sstream>
#include <unordered_map>

#include "primitive_types.h"
#include "solver.h"
//#include "top_tree_sampler.h"
#include "solver_config.h"
#include "sampler_tools.h"

StopWatch::StopWatch() {
  interval_length_.tv_sec = 60;
  gettimeofday(&last_interval_start_, nullptr);
  start_time_ = stop_time_ = last_interval_start_;
}


timeval StopWatch::getElapsedTime() {
  timeval other_time = stop_time_;
  if (stop_time_.tv_sec == start_time_.tv_sec
      && stop_time_.tv_usec == start_time_.tv_usec)
    gettimeofday(&other_time, nullptr);
  long int ad = 0;
  unsigned int bd = 0;

  if (other_time.tv_usec < start_time_.tv_usec) {
    ad = 1;
    bd = 1000000;
  }
  timeval r = (struct timeval) {0};
  r.tv_sec = other_time.tv_sec - ad - start_time_.tv_sec;
  r.tv_usec = other_time.tv_usec + bd - start_time_.tv_usec;
  return r;
}


bool Solver::simplePreProcess() {
  if (!config_.perform_pre_processing)
    return true;
  assert(literal_stack_.empty());
  // BEGIN process unit clauses
  for (auto lit : unit_clauses_)
    setLiteralIfFree(lit);
  // END process unit clauses
  VariableIndex start_ofs = 0;
  bool succeeded = BCP(start_ofs);

  if (succeeded)
    succeeded &= prepFailedLiteralTest();

  // ToDo major slowdown possible without HardWireAndCompact.  Need to make a different version.
//  if (succeeded)
//    HardWireAndCompact();
  return succeeded;
}

bool Solver::prepFailedLiteralTest() {
  unsigned long last_size;
  do {
    last_size = literal_stack_.size();
    for (unsigned v = 1; v < variables_.size(); v++) {
      if (isActive(v)) {
        unsigned long sz = literal_stack_.size();
        setLiteralIfFree(LiteralID(v, true));
        bool res = BCP(sz);
        while (literal_stack_.size() > sz) {
          unSet(literal_stack_.back());
          literal_stack_.pop_back();
        }

        if (!res) {
          sz = literal_stack_.size();
          setLiteralIfFree(LiteralID(v, false));
          if (!BCP(sz))
            return false;
        } else {
          sz = literal_stack_.size();
          setLiteralIfFree(LiteralID(v, false));
          bool resb = BCP(sz);
          while (literal_stack_.size() > sz) {
            unSet(literal_stack_.back());
            literal_stack_.pop_back();
          }
          if (!resb) {
            sz = literal_stack_.size();
            setLiteralIfFree(LiteralID(v, true));
            if (!BCP(sz))
              return false;
          }
        }
      }
    }
  } while (literal_stack_.size() > last_size);

  return true;
}

//void Solver::HardWireAndCompact() {
//  compactClauses();
//  compactVariables();
//  literal_stack_.clear();
//
//  for (auto l = LiteralID(1, false); l != literals_.end_lit(); l.inc()) {
//    literal(l).activity_score_ = literal(l).binary_links_.size() - 1;
//    literal(l).activity_score_ += occurrence_lists_[l].size();
//  }
//
//  statistics_.num_unit_clauses_ = unit_clauses_.size();
//
//  statistics_.num_original_binary_clauses_ = statistics_.num_binary_clauses_;
//  statistics_.num_original_unit_clauses_ = statistics_.num_unit_clauses_ = unit_clauses_.size();
//  initStack(num_variables());
//  original_lit_pool_size_ = literal_pool_.size();
//}

void Solver::sample_models() {
  if (!createfromFile(statistics_.input_file_)) {
    PrintFinalSamplerResults();
    return;
  }

  SampleAssignment::set_num_var(statistics_.num_variables_);

  Solver temp_solver(config_, statistics_);
  temp_solver.createfromFile(statistics_.input_file_);
  LinkConfigAndStatistics();

  reservoir_sample_models(PartialAssignment(), temp_solver);

  ReportSharpSatResults();
  PrintFinalSamplerResults();
}

void Solver::reservoir_sample_models(const PartialAssignment &partial_assn, Solver &temp_solver) {
  if (!InitializeSolverAndPreprocess(partial_assn))
    return;

  PerformInitialSampling();
  statistics_.sampler_pass_1_time_ = sampler_stopwatch_.getElapsedSeconds();

  // The UNSAT case can be reached in two different ways either by the
  // preprocessor or by running sharpSAT's main flow.
  if (statistics_.final_solution_count_ > 0 && config_.perform_two_pass_sampling_)
    FillPartialAssignments(temp_solver);

  // Print the final results after all sampling then quit.
  sampler_stopwatch_.stop();
  statistics_.sampler_time_elapsed_ = sampler_stopwatch_.getElapsedSeconds();
  statistics_.sampler_pass_2_time_ = statistics_.sampler_time_elapsed_
                                                              - statistics_.sampler_pass_1_time_;
}


void Solver::PrintFinalSamplerResults() {
  if (!config_.quiet)
    std::cout << "\nTotal Sampler Execution Time: " << statistics_.sampler_time_elapsed_ << "s\n\n";
  // Print any helper messages
  if (!config_.quiet && statistics_.final_solution_count_ < config_.num_samples_) {
      std::stringstream ss;
      ss << "Only " << statistics_.final_solution_count_ << " models exist but "
         << config_.num_samples_ << " samples were requested." << std::endl;
      PrintWarning(ss.str());
  }

  // Export to a file and the console.
  if (!config_.quiet)
    final_samples_.exportFinal(std::cout, statistics_, config_);
  if (!config_.samples_output_file.empty()) {
    std::ofstream samples_file(config_.samples_output_file);
    final_samples_.exportFinal(samples_file, statistics_, config_);
    samples_file.close();
  }
////  assert(final_samples_.VerifySolutions(statistics_.input_file_,
////                                        config_.skip_partial_assignment_fill));  // DebugZH
//  // ToDo Remove Debug Verification
//  bool valid_samples = final_samples_.VerifySolutions(statistics_.input_file_,
//                                                      config_.skip_partial_assignment_fill);
//  if (!valid_samples)
//    ExitWithError("INVALID SAMPLES", EX_SOFTWARE);
//  if (!final_samples_.IsComplete())
//    ExitWithError("INCOMPLETE SAMPLES", EX_SOFTWARE);
}


void Solver::PerformInitialSampling() {
  if (!config_.perform_two_pass_sampling_)
    config_.EnableSampleCaching();

  if (!config_.quiet) {
    if (config_.perform_two_pass_sampling_)
      std::cout << "STAGE #1: ";
    std::cout << "Build the initial partial assignments" << std::endl;
  }
  // If the solver has not been shown to be UNSAT by the preprocess,
  // run the basic sampler and model count
  statistics_.exit_state_ = countSAT();

  unused_vars_ = stack_.top().emancipated_vars();
  statistics_.set_final_solution_count(stack_.top().getTotalModelCount());
  statistics_.num_long_conflict_clauses_ = num_conflict_clauses();

  // Clean up the component so it is not an issue during pass two
  if (stack_.back().isComponentSplit())
    stack_.back().unsetAsComponentSplit();

  // If UNSAT, report it and exit.
  if (statistics_.final_solution_count_ == 0)
    return;

  if (unused_vars_.empty()) {
    if (config_.verbose)
      std::cout << "No unused variables.  Continuing..." << std::endl;
  } else {
    if (config_.verbose)
      std::cout << "Number of unused variables: " << unused_vars_.size() << std::endl;
    samples_stack_.back().AddEmancipatedVars(unused_vars_);
  }
  final_samples_ = samples_stack_.back();
  assert(IsEndSamplesStackSizeValid());
  assert(final_samples_.model_count() == statistics_.final_solution_count_);

  if (!config_.quiet) {
    if (config_.perform_two_pass_sampling_)
      std::cout << "STAGE #1: ";
    std::cout << "COMPLETED building initial partial assignments" << std::endl;
  }
}


inline void Solver::FillPartialAssignments(Solver &temp_solver) {
  if (config_.skip_partial_assignment_fill)
    return;

  if (!config_.quiet)
    std::cout << "STAGE #2 - Filling in partial assignments..." << std::endl;

  SamplesManager updated_samples = SamplesManager(config_.num_samples_, config_);

  // Free the memory to be used by the descendant solver's cache
//  unit_clauses_.clear();
  comp_manager_.initialize(literals_, literal_pool_, config_.quiet);
  deleteConflictClauses(true);

  // Group the samples so that same cache ids are tested together
  std::vector<std::pair<std::vector<CacheEntryID>, long>> grouped_samples;
  std::vector<ListOfSamples> list_elements;
  for (auto &sample : final_samples_.samples()) {
    auto new_key = sample.cache_comp_ids();
    SampleSize i = 0;
    for (i = 0; i < grouped_samples.size(); i++)
      if (grouped_samples[i].first == new_key)
        break;
    if (i < grouped_samples.size()) {
      grouped_samples[i].second += sample.sample_count();
      list_elements[i].emplace_back(sample);
    } else {
      grouped_samples.emplace_back(std::pair<std::vector<CacheEntryID>, long>(new_key,
                                                                            sample.sample_count()));
      list_elements.emplace_back();
      list_elements.back().emplace_back(sample);
    }
  }

  SampleSize group_cnt = 0, tot_count = 0;
  ListOfSamples * samples_list;
  auto sample_itr = final_samples_.samples().begin();
  while (sample_itr != final_samples_.samples().end()) {
    auto key = sample_itr->cache_comp_ids();
    long new_sample_count = 0;
    SampleSize group_idx;
    for (group_idx = 0; group_idx < grouped_samples.size(); group_idx++) {
      if (grouped_samples[group_idx].first == key) {
        new_sample_count = grouped_samples[group_idx].second;
        samples_list = &list_elements[group_idx];
        break;
      }
    }
    if (new_sample_count < 0) {
      if (sample_itr->cache_comp_ids().empty())
        ++sample_itr;
      else
        sample_itr = final_samples_.samples().erase(sample_itr);
      continue;
    }
    grouped_samples[group_idx].second = -new_sample_count;
    tot_count += new_sample_count;
    ++group_cnt;
    if (sample_itr->IsComplete()) {
      ++sample_itr;
      statistics_.numb_second_pass_vars_.push_back(0);
      if (!config_.quiet)
        std::cout << "Sample #" << group_cnt << " of " << grouped_samples.size()
                  << " is already a complete assignment.  Continuing..." << std::endl;
      continue;
    }

    // If the assignment is complete, then go to next sample
    VariableIndex num_unset_vars = sample_itr->num_unset_vars();
    statistics_.numb_second_pass_vars_.push_back(num_unset_vars);
    assert(num_unset_vars < statistics_.num_variables_);
    if (!config_.quiet)
      std::cout << "Completing sample #" << group_cnt  << " of " << grouped_samples.size()
                << " which has " << num_unset_vars << " variables unset and " << new_sample_count
                << " samples." << std::endl;

    // Build an intermediary solver
    Solver new_solver = temp_solver;
    new_solver.config_.quiet = new_solver.config_.quiet || !config_.verbose;
    PartialAssignment partial_assn;
    sample_itr->GetPartialAssignment(partial_assn);
    new_solver.config_.num_samples_ = static_cast<SampleSize>(new_sample_count);
    if (new_solver.config_.num_samples_ == 1)
      new_solver.config_.EnableSampleCaching();
    new_solver.reservoir_sample_models(partial_assn, temp_solver);

    assert(new_solver.IsEndSamplesStackSizeValid());
    assert(new_solver.stack_.top_const().getTotalModelCount()
           == new_solver.final_samples_.model_count());
    assert(new_solver.stack_.top_const().getTotalModelCount() > 0);
    assert(new_solver.final_samples_.IsComplete());

    // Take the updated sample and store it.
    if (samples_list->size() > 1)
      new_solver.final_samples_.TransferVariableAssignments(*samples_list);

    // Insert the new elements into the samples and delete the old one.
    final_samples_.samples().splice(sample_itr, new_solver.final_samples_.samples());
    sample_itr = final_samples_.samples().erase(sample_itr);
  }
  if (tot_count != config_.num_samples_)
    ExitWithError("Not all samples tested", EX_SOFTWARE);

  if (!config_.quiet)
    std::cout << "STAGE #2 - COMPLETE" << std::endl;
  // Relink the configuration and statistics objects which may have been unlinked during the
  // dependent tasks.
  LinkConfigAndStatistics();
}



bool Solver::InitializeSolverAndPreprocess(const PartialAssignment &partial_assn) {
  statistics_.set_final_solution_count(0);  // Zero out initial model count

  applyPartialAssignment(partial_assn);

  initStack(num_variables());

  if (!config_.quiet) {
    if (!config_.perform_random_sampling_)
      std::cout << "Performing Exact Model Counting..." << std::endl;
    else
      std::cout << "Performing Uniform Model Sampling..." << std::endl;
    std::cout << "Input File:  " << statistics_.input_file_ << std::endl;
    if (config_.perform_random_sampling_)
      std::cout << "Output File: " << config_.samples_output_file << std::endl;
  }

  if (!config_.quiet)
    std::cout << "\nPreprocessing ..." << std::flush;
  bool notfoundUNSAT = simplePreProcess();
  if (!config_.quiet)
    std::cout << " DONE" << std::endl;

  if (!notfoundUNSAT) {
    statistics_.exit_state_ = SolverExitState::SUCCESS;
    statistics_.final_solution_count_ = 0;
    if (!config_.quiet)
      std::cout << "\nFOUND UNSAT DURING PREPROCESSING " << std::endl;
    return notfoundUNSAT;
  }

  // If preprocessor did not find the formula to be unsatisfiable, then
  // get ready for model counting
  if (!config_.quiet)
    statistics_.printShortFormulaInfo();

  last_ccl_deletion_time_ = last_ccl_cleanup_time_ = statistics_.getTime();

  violated_clause.reserve(num_variables());

  comp_manager_.initialize(literals_, literal_pool_, config_.quiet);
  return notfoundUNSAT;
}


void Solver::applyPartialAssignment(const PartialAssignment &partial_assn) {
  // Optionally constrain the solution with unit clauses that match the partial assignment
  if (partial_assn.empty())
    return;

  assert(partial_assn.size() == num_variables() + FIRST_VAR);
  std::vector<LiteralID> literals;
  literals.emplace_back();
  for (VariableIndex var_num = FIRST_VAR; var_num <= num_variables(); var_num++) {
    AssignmentEncoding var_assn = partial_assn[var_num];
    if (var_assn != ASSN_U) {
      literals[0] = LiteralID((var_assn == ASSN_F) ? (-1 * static_cast<int>(var_num))
                                                   : static_cast<int>(var_num));
      statistics_.incorporateClauseData(literals);
      addClause(literals);
    }
  }
  statistics_.num_original_unit_clauses_ = statistics_.num_unit_clauses_ = unit_clauses_.size();
}


void Solver::solve(const PartialAssignment &partial_assn) {
  bool notfoundUNSAT = InitializeSolverAndPreprocess(partial_assn);

  if (notfoundUNSAT) {
    statistics_.exit_state_ = countSAT();
    statistics_.set_final_solution_count(stack_.top().getTotalModelCount());

    statistics_.num_long_conflict_clauses_ = num_conflict_clauses();
  }

  ReportSharpSatResults();
}

void Solver::ReportSharpSatResults() {
  statistics_.sampler_time_elapsed_ = sampler_stopwatch_.getElapsedSeconds();
  comp_manager_.gatherStatistics();
//  if (!results_output_file.empty())
//    statistics_.writeToFile(results_output_file);
  if (!config_.quiet)
    statistics_.printShort();
}

SolverExitState Solver::countSAT() {
  SolverNextAction state = SolverNextAction::RESOLVED;
  // Put the initial item on the top of the samples stack
  PushNewSamplesManagerOnStack();
  while (true) {
    while (comp_manager_.findNextRemainingComponentOf(stack_.top(), literal_stack_,
                                                      samples_stack_.back())) {
      decideLiteral();
      if (sampler_stopwatch_.timeBoundBroken()) {
        if (!config_.quiet)
          PrintError("TIMEOUT");
        exit(EXIT_TIMEOUT);
      }

      if (sampler_stopwatch_.interval_tick())
        printOnlineStats();

      while (!bcp()) {
        state = resolveConflict();
        if (state == SolverNextAction::BACKTRACK)
          break;
      }
      if (state == SolverNextAction::BACKTRACK)
        break;
    }
    state = backtrack();
    if (state == SolverNextAction::EXIT)
      return SolverExitState::SUCCESS;
    while (state != SolverNextAction::PROCESS_COMPONENT && !bcp()) {
      state = resolveConflict();
      if (state == SolverNextAction::BACKTRACK) {
        state = backtrack();
        if (state == SolverNextAction::EXIT)
          return SolverExitState::SUCCESS;
      }
    }
  }
}


void Solver::decideLiteral() {
  // Store the decision level info for determining what to do with the samples
  StackLevel *prev_top = &stack_.top();

  // establish another decision stack level
  StackLevel newLevel = StackLevel(stack_.top().currentRemainingComponent(),
                                   literal_stack_.size(),
                                   comp_manager_.component_stack_size());
  if (config_.perform_random_sampling_ && !stack_.top().isComponentSplit())
    newLevel.configureNewLevel(stack_.top(), stack_.get_decision_level());

  stack_.push_back(newLevel);

  // Manage the passing of the samples
  if (config_.perform_random_sampling_) {
    // Each component splits necessitate formula merging so create blanks
    if (prev_top->isComponentSplit()) {
      if (prev_top->isFirstComponent()) {
        PushNewSamplesManagerOnStack();
        prev_top->markFirstComponentComplete();
      }
      PushNewSamplesManagerOnStack();
    }
  }

  float max_score = -1;
  unsigned max_score_var = 0;
  // Select the variable with the highest score as the branch variable
  for (auto it = comp_manager_.superComponentOf(stack_.top()).varsBegin();
       *it != varsSENTINEL; it++) {
    float score = scoreOf(*it);
    if (score > max_score) {
      max_score = score;
      max_score_var = *it;
    }
  }
  // this assert should always hold,
  // if not then there is a bug in the logic of countSAT();
  assert(max_score_var != 0);

  // Create the literal to assign
  // Select either the negated or unnegated literal depending on
  // which form of the literal is most active.
  LiteralID theLit(max_score_var,
                   literal(LiteralID(max_score_var, true)).activity_score_
                   > literal(LiteralID(max_score_var, false)).activity_score_);

  setLiteralIfFree(theLit);
  statistics_.num_decisions_++;
  set_variable_depth_++;
  if (set_variable_depth_ > statistics_.max_branch_var_depth_)
    statistics_.max_branch_var_depth_ = set_variable_depth_;

  if (statistics_.num_decisions_ % 128 == 0)
//    if (statistics_.num_conflicts_ % 128 == 0)
    decayActivities();
  // decayActivitiesOf(comp_manager_.superComponentOf(stack_.top()));
  if (config_.verbose)
    std::cout << "Literal Stack Location #" << literal_stack_.size() - 1 << ": Variable #"
              << literal_stack_.back().var() << " assigned to "
              << ((literal_stack_.back().sign()) ? "TRUE" : "FALSE") << std::endl;
  assert(stack_.top_const().remaining_components_ofs() <= comp_manager_.component_stack_size());
}


SolverNextAction Solver::backtrack() {
  assert(stack_.top_const().remaining_components_ofs() <= comp_manager_.component_stack_size());
  do {
    if (stack_.top().branch_found_unsat())
      comp_manager_.removeAllCachePollutionsOf(stack_.top());
    else if (stack_.top().anotherCompProcessible())
      return SolverNextAction::PROCESS_COMPONENT;

    if (!stack_.top().isSecondBranch()) {
      LiteralID aLit = TOS_decLit();
      assert(stack_.get_decision_level() > 0);
      if (stack_.top().isComponentSplit())
        ProcessSampleComponentSplitBacktrack();
      // Must close branch after processing backtracking as it will affect the flow otherwise
      stack_.top().changeBranch();
      reactivateTOS();
      setLiteralIfFree(aLit.neg(), Antecedent(NOT_A_CLAUSE));
      if (config_.verbose)
        std::cout << "Literal Stack Location #" << literal_stack_.size() - 1 << ": Variable #"
                  << literal_stack_.back().var() << " switched to "
                  << (literal_stack_.back().sign() ? "TRUE" : "FALSE") << std::endl;
      return SolverNextAction::RESOLVED;
    }

    if (stack_.get_decision_level() <= 0)
      break;  // Bottomed out the decision stack so program execution is complete.
    reactivateTOS();

    assert(stack_.size() >= 2);
    // Merge the solution count up the stack.
    (stack_.end() - 2)->includeSolution(stack_.top().getTotalModelCount());
    if (config_.perform_random_sampling_)
      ProcessSampleComponentSplitBacktrack();

    // OTHERWISE:  backtrack further
    // Store the model count AND POTENTIALLY the assignment in cache.
    ProcessCacheStore();

    // Process top-tree sampling for max depth.
    // Must be after processing component splits to ensure that an immediately dependent
    // component split is cleared.
    if (config_.perform_top_tree_sampling && set_variable_depth_ == config_.max_top_tree_depth_) {
//      if (!stack_.top().HasNoDescendentModels() && HasNoUpperComponentSplit()) {
      if (!stack_.top().HasNoDescendentModels()) {
        mpz_class multiplier;
        if (stack_.on_deck().isComponentSplit())
          multiplier = 1;
        else
          multiplier = stack_.on_deck().getSamplerSolutionMultiplier();
      }
    }
    set_variable_depth_--;
    assert(set_variable_depth_ >= 0 && set_variable_depth_ <= literal_stack_.size());

    stack_.pop_back();
    // step to the next component not yet processed
    stack_.top().nextUnprocessedComponent();

    assert(stack_.top_const().remaining_components_ofs() < comp_manager_.component_stack_size()+1);
  } while (stack_.get_decision_level() >= 0);
  return SolverNextAction::EXIT;
}


void Solver::ProcessCacheStore() {
  if (stack_.size() == 1 || stack_.top().getTotalModelCount() == 0
      || !config_.perform_sample_caching()) {
    // Just store the model count as normal
    comp_manager_.cacheModelCountOf(stack_.top().super_component(),
                                    stack_.top().getTotalModelCount());
  } else {
    // Bundle the model count with a partial assignment
    StackLevel top = stack_.top();
    Component top_comp = comp_manager_.component(top.super_component());
    SampleAssignment cache_sample = top.random_cache_sample();
    comp_manager_.cacheModelCountAndAssignment(top.super_component(), top.getTotalModelCount(),
                                               cache_sample, top_comp);
    stack_.on_deck().set_cache_sample(cache_sample);
  }
}

void Solver::ProcessSampleComponentSplitBacktrack() {
  // Combine the results of a component split
  if (stack_.top().isComponentSplit()) {
    // Merge the component split with its parent
    if (config_.verbose) {
      std::stringstream ss;
      ss << "Component " << stack_.top().super_component() << " merging component split "
         << "at depth " << StackLevel::componentSplitDepth() << ".";
      PrintInColor(ss, COLOR_MAGENTA);
    }

    if (config_.store_sampled_models()) {
      SampleAssignment cached_sample;
      VariableIndex on_deck = samples_stack_.size() - 2;
      samples_stack_[on_deck].merge(samples_stack_.back(),
                                    stack_.top().getSamplerSolutionMultiplier(),
                                    stack_.top().emancipated_vars(),
                                    stack_.top().cached_comp_ids(),
                                    stack_.top().cached_assn(), cached_sample);
      samples_stack_.pop_back();
      if (config_.perform_sample_caching()) {
        // ToDo modify to support multiple sample caching
        stack_.top().ClearCachedAssn();
        stack_.top().set_cache_sample(cached_sample);
      }
//      assert(samples_stack_.back().VerifySolutions(statistics_.input_file_, true));  // DebugZH
    }
    stack_.top().unsetAsComponentSplit();
  }
  // Only close a component branch after testing both true and false paths.
  // If there is a back to back component split, may need to close both the component
  // split and the branch in a single round.
  if (stack_.size() <= 1)
    return;
  if (stack_.top().isSecondBranch() && stack_.on_deck().isComponentSplit()) {
    if (config_.verbose) {
      std::stringstream ss;
      ss << "Component branch #" << stack_.top().super_component() << " split end reached at depth "
         << StackLevel::componentSplitDepth() << ".  Stitching the sample.";
      PrintInColor(ss, COLOR_BLUE);
    }
    if (config_.store_sampled_models()) {
      samples_stack_[samples_stack_.size() - 2].stitch(samples_stack_.back());
      samples_stack_.pop_back();
//      assert(samples_stack_.back().VerifySolutions(statistics_.input_file_, true));  // DebugZH
    }
  }
}

SolverNextAction Solver::resolveConflict() {
  recordLastUIPCauses();

  if (statistics_.num_clauses_learned_ - last_ccl_deletion_time_
      > statistics_.clause_deletion_interval()) {
    deleteConflictClauses();
    last_ccl_deletion_time_ = statistics_.num_clauses_learned_;
  }

  if (statistics_.num_clauses_learned_ - last_ccl_cleanup_time_ > 100000) {
    compactConflictLiteralPool();
    last_ccl_cleanup_time_ = statistics_.num_clauses_learned_;
  }

  statistics_.num_conflicts_++;

  assert(stack_.top_const().remaining_components_ofs() <= comp_manager_.component_stack_size());

  assert(uip_clauses_.size() == 1);

  // DEBUG
  if (uip_clauses_.back().empty() && !config_.quiet)
    std::cerr << "EMPTY CLAUSE FOUND" << std::endl;
  // END DEBUG

  stack_.top().mark_branch_unsat();
  // BEGIN Backtracking
  // maybe the other branch had some solutions
  if (stack_.top().isSecondBranch()) {
    return SolverNextAction::BACKTRACK;
  }

  Antecedent ant(NOT_A_CLAUSE);
  // this has to be checked since using implicit BCP
  // and checking literals there not exhaustively
  // we cannot guarantee that uip_clauses_.back().front() == TOS_decLit().neg()
  // this is because we might have checked a literal
  // during implict BCP which has been a failed literal
  // due only to assignments made at lower decision levels
  if (uip_clauses_.back().front() == TOS_decLit().neg()) {
    assert(TOS_decLit().neg() == uip_clauses_.back()[0]);
    var(TOS_decLit().neg()).ante = addUIPConflictClause(
        uip_clauses_.back());
    ant = var(TOS_decLit()).ante;
  }
//  // RRR
//  else if (var(uip_clauses_.back().front()).decision_level
//      < stack_.get_decision_level()
//      && assertion_level_ <  stack_.get_decision_level()) {
//         stack_.top().set_both_branches_unsat();
//         return BACKTRACK;
//  }
//
//
//  // RRR
  assert(stack_.get_decision_level() > 0);
  assert(stack_.top_const().branch_found_unsat());

  // we do not have to remove pollutions here,
  // since conflicts only arise directly before
  // remaining components are stored
  // hence
  assert(
      stack_.top_const().remaining_components_ofs() == comp_manager_.component_stack_size());

  stack_.top().changeBranch();
  LiteralID lit = TOS_decLit();
  reactivateTOS();
  setLiteralIfFree(lit.neg(), ant);
  // END Backtracking
  return SolverNextAction::RESOLVED;
}

bool Solver::bcp() {
// the asserted literal has been set, so we start
// bcp on that literal
  VariableIndex start_ofs = literal_stack_.size() - 1;

// BEGIN process unit clauses
  for (auto lit : unit_clauses_)
    setLiteralIfFree(lit);
// END process unit clauses

  bool bSucceeded = BCP(start_ofs);

  if (config_.perform_failed_lit_test && bSucceeded) {
    bSucceeded = implicitBCP();
  }
  return bSucceeded;
}

bool Solver::BCP(VariableIndex start_at_stack_ofs) {
  for (VariableIndex i = start_at_stack_ofs; i < literal_stack_.size(); i++) {
    LiteralID unLit = literal_stack_[i].neg();
    // BEGIN Propagate Bin Clauses
    for (auto bt = literal(unLit).binary_links_.begin();
         *bt != SENTINEL_LIT; bt++) {
      if (isResolved(*bt)) {
        setConflictState(unLit, *bt);
        return false;
      }
      setLiteralIfFree(*bt, Antecedent(unLit));
    }
    // END Propagate Bin Clauses
    for (auto itcl = literal(unLit).watch_list_.rbegin();
         *itcl != SENTINEL_CL; itcl++) {
      bool isLitA = (*beginOf(*itcl) == unLit);
      auto p_watchLit = beginOf(*itcl) + 1 - isLitA;
      auto p_otherLit = beginOf(*itcl) + isLitA;

      if (isSatisfied(*p_otherLit))
        continue;
      auto itL = beginOf(*itcl) + 2;
      while (isResolved(*itL))
        itL++;
      // either we found a free or satisfied lit
      if (*itL != SENTINEL_LIT) {
        literal(*itL).addWatchLinkTo(*itcl);
        std::swap(*itL, *p_watchLit);
        *itcl = literal(unLit).watch_list_.back();
        literal(unLit).watch_list_.pop_back();
      } else {
        // or p_unLit stays resolved
        // and we have hence no free literal left
        // for p_otherLit remain poss: Active or Resolved
        if (setLiteralIfFree(*p_otherLit, Antecedent(*itcl))) {  // implication
          if (isLitA)
            std::swap(*p_otherLit, *p_watchLit);
        } else {
          setConflictState(*itcl);
          return false;
        }
      }
    }
  }
  return true;
}

//bool Solver::implicitBCP() {
//  static std::vector<LiteralID> test_lits(num_variables());
//  static LiteralIndexedVector<unsigned char> viewed_lits(num_variables() + 1,
//      0);
//
//  unsigned stack_ofs = stack_.top().literal_stack_ofs();
//  while (stack_ofs < literal_stack_.size()) {
//    test_lits.clear();
//    for (auto it = literal_stack_.begin() + stack_ofs;
//        it != literal_stack_.end(); it++) {
//      for (auto cl_ofs : occurrence_lists_[it->neg()])
//        if (!isSatisfied(cl_ofs)) {
//          for (auto lt = beginOf(cl_ofs); *lt != SENTINEL_LIT; lt++)
//            if (isActive(*lt) && !viewed_lits[lt->neg()]) {
//              test_lits.push_back(lt->neg());
//              viewed_lits[lt->neg()] = true;
//
//            }
//        }
//    }
//
//    stack_ofs = literal_stack_.size();
//    for (auto jt = test_lits.begin(); jt != test_lits.end(); jt++)
//      viewed_lits[*jt] = false;
//
//    statistics_.num_failed_literal_tests_ += test_lits.size();
//
//    for (auto lit : test_lits)
//      if (isActive(lit)) {
//        unsigned sz = literal_stack_.size();
//        // we increase the decLev artificially
//        // s.t. after the tentative BCP call, we can learn a conflict clause
//        // relative to the assn_ of *jt
//        stack_.startFailedLitTest();
//        setLiteralIfFree(lit);
//
//        assert(!hasAntecedent(lit));
//
//        bool bSucceeded = BCP(sz);
//        if (!bSucceeded)
//          recordAllUIPCauses();
//
//        stack_.stopFailedLitTest();
//
//        while (literal_stack_.size() > sz) {
//          unSet(literal_stack_.back());
//          literal_stack_.pop_back();
//        }
//
//        if (!bSucceeded) {
//          statistics_.num_failed_literals_detected_++;
//          sz = literal_stack_.size();
//          for (auto it = uip_clauses_.rbegin(); it != uip_clauses_.rend();
//              it++) {
//            setLiteralIfFree(it->front(), addUIPConflictClause(*it));
//          }
//          if (!BCP(sz))
//            return false;
//        }
//      }
//  }
//  return true;
//}

// this is IBCP 30.08
bool Solver::implicitBCP() {
  static std::vector<LiteralID> test_lits(num_variables());
  static LiteralIndexedVector<unsigned char> viewed_lits(num_variables() + 1, 0);

  VariableIndex stack_ofs = stack_.top().literal_stack_ofs();
  while (stack_ofs < literal_stack_.size()) {
    test_lits.clear();
    for (auto it = literal_stack_.begin() + stack_ofs;
         it != literal_stack_.end(); it++) {
      for (auto cl_ofs : occurrence_lists_[it->neg()])
        if (!isSatisfied(cl_ofs)) {
          for (auto lt = beginOf(cl_ofs); *lt != SENTINEL_LIT; lt++)
            if (isActive(*lt) && !viewed_lits[lt->neg()]) {
              test_lits.push_back(lt->neg());
              viewed_lits[lt->neg()] = true;
            }
        }
    }
    VariableIndex num_curr_lits = literal_stack_.size() - stack_ofs;
    stack_ofs = literal_stack_.size();
    for (auto jt = test_lits.begin(); jt != test_lits.end(); jt++)
      viewed_lits[*jt] = false;

    std::vector<float> scores;
    scores.clear();
    for (auto &test_lit : test_lits) {
      scores.push_back(literal(test_lit).activity_score_);
    }
    sort(scores.begin(), scores.end());
    num_curr_lits = 10 + num_curr_lits / 20;
    float threshold = 0.0;
    if (scores.size() > num_curr_lits) {
      threshold = scores[scores.size() - num_curr_lits];
    }

    statistics_.num_failed_literal_tests_ += test_lits.size();

    for (auto lit : test_lits)
      if (isActive(lit) && threshold <= literal(lit).activity_score_) {
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

        if (!bSucceeded) {
          statistics_.num_failed_literals_detected_++;
          sz = literal_stack_.size();
          for (auto it = uip_clauses_.rbegin();
               it != uip_clauses_.rend(); it++) {
            // DEBUG
            if (it->empty() && !config_.quiet)
              PrintError("EMPTY CLAUSE FOUND");
            // END DEBUG
            setLiteralIfFree(it->front(), addUIPConflictClause(*it));
          }
          if (!BCP(sz))
            return false;
        }
      }
  }

  // BEGIN TEST
//  float max_score = -1;
//  float score;
//  unsigned max_score_var = 0;
//  for (auto it =
//      component_analyzer_.superComponentOf(stack_.top()).varsBegin();
//      *it != varsSENTINEL; it++)
//    if (isActive(*it)) {
//      score = scoreOf(*it);
//      if (score > max_score) {
//        max_score = score;
//        max_score_var = *it;
//      }
//    }
//  LiteralID theLit(max_score_var,
//      literal(LiteralID(max_score_var, true)).activity_score_
//          > literal(LiteralID(max_score_var, false)).activity_score_);
//  if (!fail_test(theLit.neg())) {
//    std::cout << ".";
//
//    statistics_.num_failed_literals_detected_++;
//    unsigned sz = literal_stack_.size();
//    for (auto it = uip_clauses_.rbegin(); it != uip_clauses_.rend(); it++) {
//      setLiteralIfFree(it->front(), addUIPConflictClause(*it));
//    }
//    if (!BCP(sz))
//      return false;
//
//  }
  // END
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// BEGIN module conflictAnalyzer
///////////////////////////////////////////////////////////////////////////////////////////////

void Solver::minimizeAndStoreUIPClause(LiteralID uipLit,
                                       std::vector<LiteralID> &tmp_clause,
                                       const bool seen[]) {
  static std::deque<LiteralID> clause;
  clause.clear();
  assertion_level_ = 0;
  for (auto lit : tmp_clause) {
    if (existsUnitClauseOf(lit.var()))
      continue;
    bool resolve_out = false;
    if (hasAntecedent(lit)) {
      resolve_out = true;
      if (getAntecedent(lit).isAClause()) {
        for (auto it = beginOf(getAntecedent(lit).asCl()) + 1;
             *it != SENTINEL_CL; it++)
          if (!seen[it->var()]) {
            resolve_out = false;
            break;
          }
      } else if (!seen[getAntecedent(lit).asLit().var()]) {
        resolve_out = false;
      }
    }

    if (!resolve_out) {
      // uipLit should be the sole literal of this Decision Level
      if (var(lit).decision_level >= assertion_level_) {
        assertion_level_ = var(lit).decision_level;
        clause.push_front(lit);
      } else {
        clause.push_back(lit);
      }
    }
  }

  if (uipLit.var())
    assert(var_const(uipLit).decision_level == stack_.get_decision_level());

//  assert(uipLit.var() != 0);
  if (uipLit.var() != 0)
    clause.push_front(uipLit);
  uip_clauses_.emplace_back(std::vector<LiteralID>(clause.begin(), clause.end()));
}

void Solver::recordLastUIPCauses() {
// note:
// variables of lower dl: if seen we dont work with them anymore
// variables of this dl: if seen we incorporate their
// antecedent and set to unseen
  bool seen[num_variables() + 1];
  memset(seen, false, sizeof(bool) * (num_variables() + 1));

  static std::vector<LiteralID> tmp_clause;
  tmp_clause.clear();

  assertion_level_ = 0;
  uip_clauses_.clear();

  unsigned long lit_stack_ofs = literal_stack_.size();
  long DL = stack_.get_decision_level();
  unsigned lits_at_current_dl = 0;

  for (auto l : violated_clause) {
    if (var(l).decision_level == 0 || existsUnitClauseOf(l.var()))
      continue;
    if (var(l).decision_level < DL)
      tmp_clause.push_back(l);
    else
      lits_at_current_dl++;
    literal(l).increaseActivity();
    seen[l.var()] = true;
  }

  LiteralID curr_lit;
  while (lits_at_current_dl) {
    assert(lit_stack_ofs != 0);
    curr_lit = literal_stack_[--lit_stack_ofs];

    if (!seen[curr_lit.var()])
      continue;

    seen[curr_lit.var()] = false;

    if (lits_at_current_dl-- == 1) {
      // perform UIP stuff
      if (!hasAntecedent(curr_lit)) {
        // this should be the decision literal when in first branch
        // or it is a literal decided to explore in failed literal testing
//        assert(stack_.TOS_decLit() == curr_lit);
//        std::cout << "R" << curr_lit.toInt() << "S"
//                  << var(curr_lit).ante.isAnt() << " "  << std::endl;
        break;
      }
    }

    assert(hasAntecedent(curr_lit));

//    std::cout << "{" << curr_lit.toInt() << "}";
    if (getAntecedent(curr_lit).isAClause()) {
      updateActivities(getAntecedent(curr_lit).asCl());
      assert(curr_lit == *beginOf(antecedent(curr_lit).asCl()));

      for (auto it = beginOf(getAntecedent(curr_lit).asCl()) + 1;
           *it != SENTINEL_CL; it++) {
        if (seen[it->var()] || (var(*it).decision_level == 0)
            || existsUnitClauseOf(it->var()))
          continue;
        if (var(*it).decision_level < DL)
          tmp_clause.push_back(*it);
        else
          lits_at_current_dl++;
        seen[it->var()] = true;
      }
    } else {
      LiteralID alit = getAntecedent(curr_lit).asLit();
      literal(alit).increaseActivity();
      literal(curr_lit).increaseActivity();
      if (!seen[alit.var()] && var(alit).decision_level != 0
          && !existsUnitClauseOf(alit.var())) {
        if (var(alit).decision_level < DL)
          tmp_clause.push_back(alit);
        else
          lits_at_current_dl++;
        seen[alit.var()] = true;
      }
    }
    curr_lit = NOT_A_LIT;
  }

//  std::cout << "T" << curr_lit.toInt() << "U "
//            << var(curr_lit).decision_level << ", " << stack_.get_decision_level()
//            << "\n"
//            << "V"  << var(curr_lit).ante.isAnt() << " "  << std::endl;
  minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause, seen);

//  if (var(curr_lit).decision_level > assertion_level_)
//    assertion_level_ = var(curr_lit).decision_level;
}

void Solver::recordAllUIPCauses() {
// note:
// variables of lower dl: if seen we dont work with them anymore
// variables of this dl: if seen we incorporate their
// antecedent and set to unseen
  bool seen[num_variables() + 1];
  memset(seen, false, sizeof(bool) * (num_variables() + 1));

  static std::vector<LiteralID> tmp_clause;
  tmp_clause.clear();

  assertion_level_ = 0;
  uip_clauses_.clear();

  unsigned long lit_stack_ofs = literal_stack_.size();
  long DL = stack_.get_decision_level();
  unsigned lits_at_current_dl = 0;

  for (auto l : violated_clause) {
    if (var(l).decision_level == 0 || existsUnitClauseOf(l.var()))
      continue;
    if (var(l).decision_level < DL)
      tmp_clause.push_back(l);
    else
      lits_at_current_dl++;
    literal(l).increaseActivity();
    seen[l.var()] = true;
  }
  unsigned n = 0;
  LiteralID curr_lit;
  while (lits_at_current_dl) {
    assert(lit_stack_ofs != 0);
    curr_lit = literal_stack_[--lit_stack_ofs];

    if (!seen[curr_lit.var()])
      continue;

    seen[curr_lit.var()] = false;

    if (lits_at_current_dl-- == 1) {
      n++;
      if (!hasAntecedent(curr_lit)) {
        // this should be the decision literal when in first branch
        // or it is a literal decided to explore in failed literal testing
        //assert(stack_.TOS_decLit() == curr_lit);
        break;
      }
      // perform UIP stuff
      minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause, seen);
    }

    assert(hasAntecedent(curr_lit));

    if (getAntecedent(curr_lit).isAClause()) {
      updateActivities(getAntecedent(curr_lit).asCl());
      assert(curr_lit == *beginOf(getAntecedent(curr_lit).asCl()));

      for (auto it = beginOf(getAntecedent(curr_lit).asCl()) + 1;
           *it != SENTINEL_CL; it++) {
        if (seen[it->var()] || (var(*it).decision_level == 0)
            || existsUnitClauseOf(it->var()))
          continue;
        if (var(*it).decision_level < DL)
          tmp_clause.push_back(*it);
        else
          lits_at_current_dl++;
        seen[it->var()] = true;
      }
    } else {
      LiteralID alit = getAntecedent(curr_lit).asLit();
      literal(alit).increaseActivity();
      literal(curr_lit).increaseActivity();
      if (!seen[alit.var()] && var(alit).decision_level != 0
          && !existsUnitClauseOf(alit.var())) {
        if (var(alit).decision_level < DL)
          tmp_clause.push_back(alit);
        else
          lits_at_current_dl++;
        seen[alit.var()] = true;
      }
    }
  }
  if (!hasAntecedent(curr_lit)) {
    minimizeAndStoreUIPClause(curr_lit.neg(), tmp_clause, seen);
  }
//  if (var(curr_lit).decision_level > assertion_level_)
//    assertion_level_ = var(curr_lit).decision_level;
}

void Solver::printOnlineStats() {
  if (config_.quiet)
    return;

  std::cout << "\ntime elapsed: " << sampler_stopwatch_.getElapsedSeconds() << "s" << std::endl;
  if (config_.verbose) {
    std::cout << "conflict clauses (all / bin / unit) \t"
              << num_conflict_clauses()
              << "/" << statistics_.num_binary_conflict_clauses_ << "/"
              << unit_clauses_.size() << "\n"
              << "failed literals found by implicit BCP \t "
              << statistics_.num_failed_literals_detected_ << "\n";

    std::cout << "implicit BCP miss rate \t "
              << statistics_.implicitBCP_miss_rate() * 100 << "%\n";

    comp_manager_.gatherStatistics();

    std::cout << "cache size " << statistics_.cache_MB_memory_usage() << "MB" << "\n"
              << "components (stored / hits) \t\t"
              << statistics_.cached_component_count() << "/"
              << statistics_.cache_hits() << "\n"
              << "avg. variable count (stored / hits) \t"
              << statistics_.getAvgComponentSize() << "/"
              << statistics_.getAvgCacheHitSize() << "\n"
              << "cache miss rate " << statistics_.cache_miss_rate() * 100 << "%"
              << std::endl;
  }
}


//bool Solver::IsVarInLiteralStack(const VariableIndex var) {
//  for (VariableIndex i = 0; i < literal_stack_.size(); i++) {
//    if (literal_stack_[i].var() == var) {
//      std::stringstream ss;
//      ss << "Var #" << var << " is in the literal stack at level " << i
//          << " with val " << literal_stack_[i].sign();
//      PrintError(ss);
//      return true;
//    }
//  }
//  return false;
//}
//
//
//void Solver::PrintLiteralStackLocation(VariableIndex var) {
//  for (VariableIndex i = 0; i < literal_stack_.size(); i++) {
//    if (literal_stack_[i].var() == var) {
//      std::stringstream ss;
//      ss << "Var #" << var << " is located in slot " << i << " and has sign "
//         << (literal_stack_[i].sign() ? "POS" : "NEG");
//      PrintError(ss);
//    }
//  }
//}

Solver::Solver(int argc, char *argv[]) : final_samples_(0, config_) {
  // Store the configuration for visibility by the stack.
  LinkConfigAndStatistics();

  time(&statistics_.start_time_);
  config_.num_samples_ = 0;  // By default initialize the model count to zero
  for (int i = 1; i < argc; i++) {
//    if (strcmp(argv[i], "-noCC") == 0)
//      config_.perform_component_caching = false;
//    else if (strcmp(argv[i], "-noIBCP") == 0)
//      config_.perform_failed_lit_test = false;
//    else if (strcmp(argv[i], "-noPP") == 0)
//      config_.perform_pre_processing = false;
//    else if (strcmp(argv[i], "-q") == 0)
    if (strcmp(argv[i], "-q") == 0) {
      config_.quiet = true;
    } else if (strcmp(argv[i], "-v") == 0) {
      config_.verbose = true;
    } else if (strcmp(argv[i], "-d") == 0) {
      config_.debug_mode = true;
    } else if (strcmp(argv[i], "-tp") == 0) {
      config_.perform_two_pass_sampling_ = true;
    } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "-out") == 0) {
      if (argc <= i + 1 || !config_.perform_random_sampling_)
        ExitInvalidParam("Invalid parameters for sampling");
      if (strcmp(argv[i], "-s") == 0) {
        long num_samples = strtoul(argv[++i], nullptr, STR_DECIMAL_BASE);
        config_.num_samples_ = static_cast<SampleSize>(num_samples);
        if (config_.num_samples_ < 1)
          ExitInvalidParam("Must sample at least one model");
      } else {
        config_.samples_output_file = argv[++i];
      }
    } else if (strcmp(argv[i], "-no-sample-write") == 0) {
      config_.disable_samples_write_ = true;
    } else if (strcmp(argv[i], "-count-only") == 0) {
      config_.perform_random_sampling_ = false;
      if (config_.num_samples_ > 0 || !config_.samples_output_file.empty())
        ExitInvalidParam("Invalid parameters for counting and sampling");
    } else if (strcmp(argv[i], "-t") == 0) {
      if (argc <= i + 1)
        ExitInvalidParam("Time bound missing");
      config_.time_bound_seconds = strtoul(argv[++i], nullptr, STR_DECIMAL_BASE);
    } else if (strcmp(argv[i], "-cs") == 0) {
      if (argc <= i + 1)
        ExitInvalidParam("No cache size specified");
      statistics_.maximum_cache_size_bytes_ = strtoul(argv[++i], nullptr, STR_DECIMAL_BASE)
                                              * (uint64_t) 1000000;
    } else if (strcmp(argv[i], "-cnf") == 0) {
      if (argc <= i + 1)
        ExitInvalidParam("No CNF file specified");
      statistics_.input_file_ = argv[++i];
    } else if (strcmp(argv[i], "-no-partial-fill") == 0) {
      // This is a debug only feature. Functionality is not guaranteed.
      config_.skip_partial_assignment_fill = true;
    } else if (strcmp(argv[i], "-top-tree") == 0) {
      config_.perform_top_tree_sampling = true;
    } else if (strcmp(argv[i], "-top-tree-depth") == 0) {
      long top_tree_depth = strtoul(argv[++i], nullptr, STR_DECIMAL_BASE);
      if (top_tree_depth <= 1)
        ExitInvalidParam("The top tree depth must be greater than 1.");
      config_.max_top_tree_depth_ = static_cast<TreeNodeIndex>(top_tree_depth);
    } else if (strcmp(argv[i], "-max-leaf-size") == 0) {
      long max_leaf_size = strtoul(argv[++i], nullptr, STR_DECIMAL_BASE);
      if (max_leaf_size <= 1)
        ExitInvalidParam("The top tree leaf size must be greater than 1.");
      config_.max_top_tree_leaf_sample_count = static_cast<SampleSize>(max_leaf_size);
    } else {
      ExitInvalidParam(static_cast<std::string>("Unknown parameter found \"") + argv[i] + "\"");
    }
  }

  // Perform additional cross-checking of input parameters
  if (config_.quiet && config_.verbose)
    ExitInvalidParam("Invalid combination of verbose and quiet");
  if (config_.num_samples_ < 1 && config_.perform_random_sampling_)
    ExitInvalidParam("If sampling is enabled, a sample count greater than or equal \n"
                      "to 1 must be specified.");
  if (config_.perform_two_pass_sampling_ && !config_.perform_random_sampling_)
    ExitInvalidParam("The two pass sampling flag is only applicable when sampling is enabled.\n");
  if (config_.perform_top_tree_sampling && !config_.perform_random_sampling_)
    ExitInvalidParam("Top tree sampling cannot be enabled if random sampling is disabled.");
  if (config_.disable_samples_write_ && !config_.perform_random_sampling_)
    ExitInvalidParam("Disabling sample write cannot be selected if sampling is disabled.");

  if (config_.samples_output_file.empty() && config_.perform_random_sampling_) {
    std::string input_file = statistics_.input_file_;
    uint64_t os_sep_loc = input_file.rfind('\\');
    os_sep_loc = (os_sep_loc != std::string::npos) ? os_sep_loc : input_file.rfind('/');

    std::string out_file;
    if (os_sep_loc != std::string::npos)
      out_file = statistics_.input_file_.substr(0, os_sep_loc + 1);
    // Prepend "results_" to the samples filename
    out_file += "samples_";

    // Append the filename and change extension to ".txt"
    uint64_t filename_start_loc = os_sep_loc + 1;
    uint64_t file_ext_start = input_file.rfind('.');
    if (file_ext_start <= filename_start_loc || file_ext_start == std::string::npos)
      file_ext_start = statistics_.input_file_.size();
    uint64_t filename_len = file_ext_start - filename_start_loc;
    out_file += statistics_.input_file_.substr(filename_start_loc, filename_len);
    out_file += ".txt";
    config_.samples_output_file = out_file;
    if (!config_.quiet)
      PrintWarning(static_cast<std::string>("No sample results file specified.\n")
                   + "Using default filename: \"" + out_file + "\"");
  }
  if (config_.perform_two_pass_sampling_ && config_.num_samples_ > 1 && !config_.quiet)
    PrintWarning("The two pass sampling flag only has an effect\n"
                     "when the number of samples equals 1. Ignoring...");

  // Initialize the time bound.
  sampler_stopwatch_.setTimeBound(config_.time_bound_seconds);

  // Any time the same count is larger than one, always do two pass sampling
  if (config_.num_samples_ > 1)
    config_.perform_two_pass_sampling_ = true;
  // Must parse at the end in case debug mode is selected.
  Random::init(&config_);
}
