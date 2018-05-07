/**
 * model_sampler.cpp
 *
 * Purpose: Classes for storing and managing partial and complete SAT assignments.
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

#include <gmpxx.h>
#include <fstream>

// Used for precision string printing
#include <iomanip>
#include <sstream>
#include <algorithm>

#include "stack.h"
#include "model_sampler.h"
//#include "top_tree_sampler.h"


// Satisfy the linker by initializing the static variables
VariableIndex SampleAssignment::num_var_ = 0;
VariableIndex SampleAssignment::var_vec_len_ = 0;

void SamplesManager::exportFinal(std::ostream &out,
                                 const DataAndStatistics& statistics,
                                 const SolverConfiguration& config) {
  time_t ltime = time(&ltime);
  tm cur_time;
  localtime_r(&ltime, &cur_time);

  // Create the output file header
  char buf[50];
  out << "#START_HEADER" << "\n"
      << "start_time," << asctime_r(&cur_time, buf)  // Implicit new line.
      << "formula_file," << statistics.input_file_ << "\n"
      << "num_vars," << statistics.num_variables_ << "\n"
      << "independent_support_size," << statistics.num_indep_support_variables_ << "\n"
      << "num_clauses," << statistics.num_clauses() << "\n"
      << "tot_num_models," << statistics.final_solution_count().get_str() << "\n"
      << "max_component_split_depth," << statistics.max_component_split_depth_ << "\n"
      << "max_branch_var_depth," << statistics.max_branch_var_depth_ << "\n"
      << "num_samples," << config.num_samples_ << "\n"
      << "num_second_pass_groups," << statistics.numb_second_pass_vars_.size() << "\n"
      << "num_second_pass_vars";

  VariableIndex average_rem_var_cnt = 0;
  for (auto rem_var_cnt : statistics.numb_second_pass_vars_) {
    out << "," << rem_var_cnt;
    average_rem_var_cnt += rem_var_cnt;
  }
  std::stringstream stream;
  stream << std::fixed << std::setprecision(2)
         << static_cast<double>(average_rem_var_cnt)/statistics.numb_second_pass_vars_.size();
  out << "\n";
//  TopTreeSampler::printTreeSolverStatistics(out);
  out << "avg_second_pass_var," << stream.str() << "\n"
      << "execution_time," << statistics.sampler_time_elapsed_ << "\n"
      << "pass1_time," << statistics.sampler_pass_1_time_ << "\n"
      << "pass2_time," << statistics.sampler_pass_2_time_ << "\n"
      << "#END_HEADER\n";

  // Do not write the solutions to the console.
  if (&out == &std::cout || config.disable_samples_write_)
    return;

  // Write the samples (if any
  out << "\n\n#START_SAMPLES\n";
  if (statistics.final_solution_count_ > 0) {
    for (const auto &sample : samples_)
      out << sample.sample_count() << "," << sample.ToString() << "\n";
  } else {
    out << "UNSAT\n";
  }
  out << "#END_SAMPLES" << std::endl;
}


void SamplesManager::reservoirSample(const Component * active_comp,
                                     const std::vector<LiteralID> &literal_stack,
                                     const mpz_class &solution_weight,
                                     const mpz_class &weight_multiplier,
                                     const AltComponentAnalyzer &ana,
                                     const VariableIndex literal_stack_ofs,
                                     const std::vector<VariableIndex> &freed_vars,
                                     const std::vector<CacheEntryID> &cached_comp_ids,
                                     const CachedAssignment &cached_assn,
                                     SampleAssignment& cached_sample) {
  assert(num_samples() > 0);
  if (solution_weight == 0)
    return;  // UNSAT so move along.

  mpz_class sample_weight = solution_weight;
  if (weight_multiplier != 1)
    sample_weight *= weight_multiplier;
  solution_count_ += sample_weight;

  SampleSize num_samples_to_replace;
  if (solution_count_ == sample_weight) {
    // First solution always gets weight equivalent to the number of samples.
    assert(samples_.empty());

    num_samples_to_replace = num_samples();
    if (config_->verbose) {
      std::stringstream ss;
      ss << "\t\tInitial sample auto-selected.";
      PrintInColor(ss, COLOR_GREEN);
    }
  } else {
    // Use a binomial random variable to determine how many samples to replace.
    std::vector<SampleSize> samples_to_replace;
    GenerateSamplesToReplace(sample_weight, samples_to_replace);
    num_samples_to_replace = static_cast<SampleSize>(samples_to_replace.size());
    // Nothing to replace in this case.
    if (num_samples_to_replace == 0) {
      if (config_->perform_sample_caching()) {
        BuildSample(cached_sample, active_comp, literal_stack,
                    literal_stack_ofs, ana, freed_vars, cached_comp_ids);
        cached_sample.IncorporateCachedAssignment(cached_assn);
        assert(cached_sample.VerifyEmancipatedVars());
      }
      return;
    }

    RemoveSamples(samples_to_replace);
    if (config_->verbose) {
      std::stringstream ss;
      ss << "\t\tReservoir sampling accepted " << samples_to_replace.size() << " samples.";
      PrintInColor(ss, COLOR_GREEN);
    }
  }

  samples_.emplace_back(SampleAssignment(num_samples_to_replace));
  BuildSample(samples_.back(), active_comp, literal_stack,
              literal_stack_ofs, ana, freed_vars, cached_comp_ids);

  // Handle the special requirements of sample caching.
  if (config_->perform_sample_caching()) {
    samples_.back().IncorporateCachedAssignment(cached_assn);
    cached_sample = SampleAssignment(config_->num_samples_to_cache_, samples_.back());
    assert(cached_sample.VerifyEmancipatedVars());
    assert(num_samples() == cached_sample.sample_count());
  }
}

void SamplesManager::GenerateSamplesToReplace(const mpz_class &new_sample_weight,
                                              std::vector<SampleSize> &samples_to_replace) const {
  assert(new_sample_weight <= solution_count_);
  SampleSize num_samples_to_replace = Random::Mpz::binom(num_samples(), solution_count_,
                                                         new_sample_weight);
  samples_to_replace.clear();
  if (num_samples_to_replace == 0)
    return;

  Random::SelectRangeInts(num_samples(), num_samples_to_replace, samples_to_replace);
  // This sorts in DESCENDING ORDER.  Hence, it would be <5, 4, 3, 2, 1>.  This is for
  // Simplified removal later.
  std::sort(samples_to_replace.rbegin(), samples_to_replace.rend());
}



bool SamplesManager::VerifySolutions(const std::string &input_file_path,
                                     const bool skip_unassigned) const {
  std::vector<std::vector<signed long>> clauses;
  buildCnfClauseLiterals(input_file_path, clauses);
  bool cls_sat;
  for (const auto &sample : samples_) {
    for (SampleSize i = 0; i < sample.sample_count(); i++) {
      PartialAssignment assn;
      sample.BuildRandomizedPartialAssignment(assn);
      for (auto cls : clauses) {
        cls_sat = false;
        for (auto lit : cls) {
          // Any unassigned literals automatically fail.
          if (assn[abs(lit)] == ASSN_U) {
            if (skip_unassigned) {// If true, ignore all clauses with unassigned variables
              cls_sat = true;
              break;
            }
            std::cerr << "Error: Variable #" << abs(lit) << " Unassigned bit in final sample\n";
            return false;
          }
          // One true literal satisfies the clause so we can break
          if ((lit < 0 && assn[abs(lit)] == ASSN_F)
              || (lit > 0 && assn[abs(lit)] == ASSN_T)) {
            cls_sat = true;
            break;
          }
        }
        if (cls_sat)
          continue;
        std::stringstream ss;
        ss << "Assignment is invalid.  Clause is ";
        for (auto lit : cls)
          ss << lit << ", ";
        PrintError(ss);
        return false;
      }
    }
  }
  return true;
}

void SamplesManager::BuildSample(SampleAssignment &new_sample,
                                 const Component *active_comp,
                                 const std::vector<LiteralID> &literal_stack,
                                 const VariableIndex last_branch_lit,
                                 const AltComponentAnalyzer &ana,
                                 const std::vector<VariableIndex> &freed_vars,
                                 const std::vector<CacheEntryID> &cached_comp_ids) {
  assert(last_branch_lit < literal_stack.size() || literal_stack.empty());
  assert(new_sample.num_set_vars() == 0);

  // Literal stack values are constrained so explicitly set them
  for (auto lit : literal_stack)
    new_sample.setVarAssignment(lit.var(), lit.sign() ? ASSN_T : ASSN_F);
  new_sample.addEmancipatedVars(freed_vars);
  new_sample.addCachedCompIds(cached_comp_ids);

  // Get all the free/unassigned variables
  std::vector<VariableIndex> new_free_vars;
  VariableIndex i = 0;
  const auto &component_data = active_comp->getData();
  // Run to the end of the variables
  while (component_data[i] != varsSENTINEL) {
    VariableIndex var_num = component_data[i];
    i++;
    // Ignore all non-free variables.
    // Sometimes variables set in the last or in BCP may appear free but are not
    if (ana.scoreOf(var_num) > 0)
      continue;
    if (SamplesManager::isVarInLitStack(var_num, literal_stack, last_branch_lit))
      continue;
    new_free_vars.push_back(var_num);
    // DebugZH - Make sure literal was not missed in the literal stack.
    assert(!SamplesManager::isVarInLitStack(var_num, literal_stack));
  }
  new_sample.addEmancipatedVars(new_free_vars);
}

void SamplesManager::buildCnfClauseLiterals(const std::string &input_file_path,
                                            std::vector<std::vector<signed long>> &clauses) {
  // BEGIN File input
  std::ifstream input_file(input_file_path);
  if (!input_file)
    ExitWithError("Cannot open file " + input_file_path + " Solution validation failed.", EX_IOERR);

  // Remove the comments header
  VariableIndex nVars;
  ClauseIndex nCls;
  const uint64_t max_ignore = UINT64_MAX;  // Max number of characters on line to ignore.
  std::string idstring;
  char c;
  while (input_file >> c && c != 'p')
    input_file.ignore(max_ignore, '\n');
  if (!(input_file >> idstring && idstring == "cnf" && input_file >> nVars && input_file >> nCls))
    ExitWithError("Invalid CNF file\n" "Solution validation failed.", EX_DATAERR);

  // Analyze each clause and add literals/variables as appropriate.
  long lit;
  clauses.clear();
  while ((input_file >> c) && clauses.size() < nCls) {
    // Ignore comment inline
    if (c == 'c') {
      input_file.ignore(max_ignore, '\n');
      continue;
    }
    clauses.emplace_back();  // Create a new clause
    input_file.unget();  // Extracted a non-space character to determine if a clause, so put it back
    if ((c == '-') || isdigit(c)) {
      while ((input_file >> lit) && lit != 0) {
        clauses.back().push_back(lit);
      }
    }
    input_file.ignore(max_ignore, '\n');
    if (clauses.back().empty())
      clauses.pop_back();
  }
  input_file.close();
  assert(clauses.size() == nCls);
}

void SamplesManager::merge(SamplesManager &other, const mpz_class &other_multiplier,
                           const std::vector<VariableIndex> &freed_vars,
                           const std::vector<CacheEntryID> &cached_comp_ids,
                           const CachedAssignment & cached_assn,
                           SampleAssignment& cached_sample) {
  // No solutions in the component split so move on.
  if (other.solution_count_ == 0 || other_multiplier == 0)
    return;

  // Determine which samples to take from the other component split.
  mpz_class other_weighted = other.solution_count_;
  if (other_multiplier > 1)
    other_weighted *= other_multiplier;
  solution_count_ += other_weighted;

  other.AddEmancipatedVars(freed_vars);
  other.AddCachedCompIds(cached_comp_ids);
  // Handle the case of sample caching.
  if (config_->perform_sample_caching()) {
    assert(config_->num_samples_to_cache_ == 1);
    other.samples_.front().IncorporateCachedAssignment(cached_assn);
    cached_sample = other.samples_.front();
  }

  std::vector<SampleSize> samples_to_replace;
  GenerateSamplesToReplace(other_weighted, samples_to_replace);
  if (samples_to_replace.empty())
    return;
  RemoveSamples(samples_to_replace);

  // What is removed from the existing set is what is kept from the other set.
  other.KeepSamples(samples_to_replace);

  if (config_->verbose) {
    std::stringstream ss;
    ss << "\t\tMerged " << samples_to_replace.size() << " samples from the component branches.";
    PrintInColor(ss, COLOR_GREEN);
  }
  append(other);
  assert(GetActualSampleCount() == num_samples());
}


//bool IsVarInLiteralStack(const std::vector<LiteralID> &literal_stack, VariableIndex var) {
//  for (unsigned i = 0; i < literal_stack.size(); i++)
//    if (literal_stack[i].var() == var) {
//      std::stringstream ss;
//      ss << "Var #" << var << " is in the literal stack at level " << i << " with val "
//         << (literal_stack[i].sign());
//      PrintInColor(std::cerr, ss.str(), COLOR_RED);
//      return true;
//    }
//  return false;
//}
