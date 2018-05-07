/**
 * model_sampler.h
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

#ifndef SHARPSAT_SOLUTION_RECIPE_H
#define SHARPSAT_SOLUTION_RECIPE_H

#include <gmpxx.h>

#include <random>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>

#include "structures.h"
#include "component_types/component.h"
#include "statistics.h"
#include "alt_component_analyzer.h"
#include "solver_config.h"
#include "cached_assignment.h"
#include "rand_distributions.h"

/**
 * Defines the encoding for variable states in the assignment.
 */

typedef std::vector<uint32_t> AssignmentContainer;

class SampleAssignment {
  friend class SamplesManager;
 private:
  /**
   * Creates a sample assignment with the specified count but that is a copy of the other
   * specified SampleAssignment.
   *
   * @param sample_count Number of samples this sample contains.
   * @param other Another sample object whose contents (other than sample count) will be duplicated
   */
  explicit SampleAssignment(SampleSize sample_count, const SampleAssignment &other) :
      sample_count_(sample_count), assn_(other.assn_), num_vars_set_(other.num_vars_set_),
      emancipated_vars_(other.emancipated_vars_), cache_comp_ids_(other.cache_comp_ids_) { }
  /**
   * Num of samples stored.
   */
  SampleSize sample_count_;
  /**
   * Recipes are stored as a set bit vectors.  This will simplify concatenation
   * and minimizes space and time overhead.
   */
  AssignmentContainer assn_;
  /**
   * Number of variables set.
   */
  VariableIndex num_vars_set_ = 0;
  /**
   * List of variables that are unset and not emancipated.
   */
  std::vector<VariableIndex> remaining_vars_;
  /**
   * Variables in the assignment that are free and can be set to any value.
   */
  std::vector<VariableIndex> emancipated_vars_;
  /**
   * Cache Identification numbers for all cached components in this sample.
   */
  std::vector<CacheEntryID> cache_comp_ids_;
  /**
   * Key for uniquely identifying the components cached in the object.  Memoized to reduce the
   * need for recomputation.
   */
  std::string comp_cache_key_;
  /**
   * Number of variables in the original file
   */
  static VariableIndex num_var_;
  /**
   * Number of words in the vector representing the sample model.
   */
  static VariableIndex var_vec_len_;
  /**
   * Number of bits used in the assignment word encoding.
   */
  static const unsigned int size_of_word_ = sizeof(typename AssignmentContainer::value_type)
                                            * BITS_PER_BYTE;
  /**
   * Number of bits used to encode each bit in the assignment structure.
   *
   * Other values for the bits per var are not supported and would require additional
   * logic added to the code to handle variables that straddle two words.
   */
  static const VariableIndex bits_per_var_ = 2;
  /**
   * Bit mask used in bitwise addition when extracting a variable value
   * from the assignment structure.
   */
  static const AssignmentContainer::value_type var_mask_ = 0x3;
  /**
   * Recipe Stitching Compatibility Checker
   *
   * Given the implicit and explicitly partial assignment parameters, this function
   * verifies that there are no corresponding variable pairings are assigned to 0
   * and 1.
   *
   * This is generally a debug-only verification function.
   *
   * @param other A second partial assignment that is to be stitched with the
   * implicit parameter.
   *
   * @return True if the partial assignments are compatible.
   */
  bool VerifyStitchingCompatibility(const SampleAssignment &other) const {
    assert(assn_.size() == other.assn_.size());
    assert(sample_count_ <= other.sample_count_);  // Less than or equal for simplified code
                                                  // when stitching the SamplesManager objects

    AssignmentEncoding this_val, other_val;
    for (VariableIndex i = FIRST_VAR; i <= num_var_; i++) {
      this_val = this->var_assignment(i);
      other_val = other.var_assignment(i);
      if (this_val != other_val && this_val != ASSN_U && other_val != ASSN_U) {
        std::cerr << "Stitching error on variable #" << i << std::endl
                  << "Implicit Value = " << this_val << ". Other Value = "
                  << other_val << std::endl;
        return false;
      }
    }
    return true;
  }
  /**
   * Helper function to standardize calculating the location of a variable in a
   * assignment.
   *
   * @param var_numb Variable identification number from 1 to num_var inclusive.
   * @param word_num Word number
   * @param bit_num Starting bit of the variable's value
   */
  inline static void calculateWordAndBitNumbers(const VariableIndex &var_num,
                                                VariableIndex &word_num, VariableIndex &bit_num) {
    assert(var_num >= FIRST_VAR && var_num <= num_var_);
    word_num = (var_num * bits_per_var_) / size_of_word_;
    assert(word_num >= 0 && word_num < var_vec_len_);
    bit_num = (var_num * bits_per_var_) % size_of_word_;
  }
  /**
   * Checks that the assignment has no duplicated emancipated variables.
   *
   * @return True if there are no duplicate emancipated variables.
   */
  const bool VerifyEmancipatedVars() const {
    if (emancipated_vars_.empty())
      return true;

    std::vector<VariableIndex> all_vars = emancipated_vars_;
    std::sort(all_vars.begin(), all_vars.end());
    if (var_assignment(all_vars[0]) != ASSN_U)
      return false;
    for (VariableIndex i = 1; i < all_vars.size(); i++) {
      if (all_vars[i] == all_vars[i-1] || var_assignment(all_vars[i]) != ASSN_U)
        return false;
    }
    return true;
  }
  /**
   * Checks that the assignment has no duplicated emancipated variables.
   *
   * @return True if there are no duplicate emancipated variables.
   */
  template <typename T>
  static const bool VerifyNoDuplicates(const std::vector<T> &vec) {
    if (vec.empty())
      return true;

    std::vector<T> vec_copy = vec;
    std::sort(vec_copy.begin(), vec_copy.end());
    for (VariableIndex i = 1; i < vec_copy.size(); i++) {
      if (vec_copy[i] == vec_copy[i-1])
        return false;
    }
    return true;
  }
  /**
   * Accessor for the unemancipated variables.
   *
   * @return List of emancipated variables.
   */
  const std::vector<VariableIndex>& emancipated_vars() const { return emancipated_vars_; }
  /**
   * Deletes the specified variables from the emancipated variables.
   *
   * @param vars_to_delete List of variables to delete.
   */
  void DeleteEmancipatedVars(const std::vector<VariableIndex> &vars_to_delete) {
    for (auto &var_id : vars_to_delete) {
      auto itr = std::find(emancipated_vars_.begin(), emancipated_vars_.end(), var_id);
      emancipated_vars_.erase(itr);
    }
  }
  /**
   * Adds a set of emancipated variables to an object.
   *
   * @param new_emancipated_vars The new emancipated variables to add.
   */
  void addEmancipatedVars(const std::vector<VariableIndex> &new_emancipated_vars) {
    if (new_emancipated_vars.empty())
      return;

    emancipated_vars_.insert(emancipated_vars_.end(), new_emancipated_vars.begin(),
                             new_emancipated_vars.end());
//    for (auto var_id : new_emancipated_vars) {
//      if (!IsVarEmancipated(var_id))
//        emancipated_vars_.push_back(var_id);
//    }
    assert(VerifyEmancipatedVars());
  }
  /**
   * Adds one set of cached component IDs to the current object.
   *
   * @param cache_comp_ids
   */
  void addCachedCompIds(const std::vector<CacheEntryID> &cache_comp_ids) {
    comp_cache_key_.clear();
    cache_comp_ids_.insert(cache_comp_ids_.end(), cache_comp_ids.begin(), cache_comp_ids.end());
    assert(VerifyNoDuplicates<CacheEntryID>(cache_comp_ids_));
  }
  /**
   * Decrease the number of samples associated with this sample.
   *
   * @param dec_sample_count Number of samples by which this sample will be reduced
   */
  void DecreaseSampleCount(SampleSize dec_sample_count) {
    assert(dec_sample_count > 0 && dec_sample_count < sample_count_);
    sample_count_ -= dec_sample_count;
  }
  /**
   * Sets the sample count for the SampleAssignment to zero.
   */
  void zeroSampleCount() { sample_count_ = 0; }
  /**
   * Number of Variables Accessor
   *
   * Returns the total number of variables in the assignment.
   *
   * @return Variable count.
   */
  inline const VariableIndex num_var() const { return num_var_; }
  /**
   * Cached Assignment Incorporator
   *
   * Incorporates the cached assignment into the sample assignment.
   *
   * @param cached_assn A cached variable assignment.
   */
  void IncorporateCachedAssignment(const CachedAssignment & cached_assn) {
    for (auto lit : cached_assn.literals()) {
      assert(var_assignment(lit.var()) == ASSN_U);  // Only assign to unassigned variables
      setVarAssignment(lit.var(), lit.sign() ? ASSN_T : ASSN_F);
    }
    if (!cached_assn.emancipated_vars().empty())
      emancipated_vars_.insert(emancipated_vars_.end(), cached_assn.emancipated_vars().begin(),
                               cached_assn.emancipated_vars().end());
    assert(VerifyEmancipatedVars());
  }
  /**
   * Stitches together two solution recipes.  It merges the other explicitly supplied
   * recipe and uses it to update the implicit recipe.
   *
   * @param other Recipe with which to stitch.
   */
  inline void stitch(const SampleAssignment &other) {
    assert(VerifyStitchingCompatibility(other));
    for (VariableIndex i = 0; i < assn_.size(); i++)
      assn_[i] &= other.assn_[i];

    addEmancipatedVars(other.emancipated_vars_);
    addCachedCompIds(other.cache_comp_ids_);
    num_vars_set_ += other.num_vars_set_;

    remaining_vars_.clear();
    // Make sure no variables miraculously materialized.
    assert(num_set_vars_const() + emancipated_vars_.size() <= num_var());
  }
  /**
   * Splits the existing samples object in two with the existing object's sample count
   * decreased by the specified sample count and the sample's count being specified.  Other than
   * that it is an identical copy of the previous assignment.
   *
   * @param new_assn_sample_count Number of samples in the new split off sample.  This is the amount
   * the implicit object's sample count will be reduced by.
   */
  inline SampleAssignment split(SampleSize new_assn_sample_count) {
    DecreaseSampleCount(new_assn_sample_count);
    return SampleAssignment(new_assn_sample_count, *this);
  }
  /**
   * Updates the sample count to the specified value.
   *
   * @param new_sample_count New sample count for the sample object.
   */
  inline void set_sample_count(SampleSize new_sample_count) { sample_count_ = new_sample_count; }
  /**
   * Builds a partial assignment with the emancipated variables set for verification purposes.
   *
   * @param all_vars Partial assignment
   */
  inline void BuildRandomizedPartialAssignment(PartialAssignment &all_vars) const {
    GetPartialAssignment(all_vars);
    for (auto var : emancipated_vars_)
      all_vars[var] = (Random::uniform(0, 1)) ? ASSN_F : ASSN_T;
  }
  /**
   * Helper function to create a temporary object that can be used to unset the assigned variables
   * for another SampleAssignment() object.
   *
   * @param vars_to_keep List of variables that would not be zeroed out by this assignment.
   *
   * @return An object that can be used to unset the variable assignments of another object.
   */
  static SampleAssignment buildUnsetterAssignment(std::vector<VariableIndex> vars_to_keep) {
    SampleAssignment new_assn;
    for (auto &var : vars_to_keep)
      new_assn.setVarAssignment(var, ASSN_F);
    return new_assn;
  }
  /**
   * Unsets multiple variable assignments in the associated assignment.
   *
   * @param unsetter Object that contains the information on which variables to unset.
   */
  void unsetVariableAssignments(SampleAssignment &unsetter) {
    assert(assn_.size() == unsetter.assn_.size());
    for (VariableIndex i = 0; i < assn_.size(); i++)
      assn_[i] |= unsetter.assn_[i];
  }
  /**
   * Gets the variables that remain to be set or emancipated.
   *
   * @return Variables that remain to be set or emancipated.
   */
  std::vector<VariableIndex> GetRemainingVariables() {
    if (remaining_vars_.empty())
      num_set_vars();  // Builds the unset variables list.
    return remaining_vars_;
  }

 public:
  /**
   * Creates a blank assignment with all variables unset.
   *
   * @param sample_count Number of samples this sample will represent.
   */
  explicit SampleAssignment(SampleSize sample_count) : sample_count_(sample_count) {
    assn_.resize(var_vec_len_, static_cast<typename AssignmentContainer::value_type>(-1));
  }
  /**
   * Empty constructor with zero samples.  It is not valid for general use.
   */
  SampleAssignment() : SampleAssignment(0) {}
  /**
   * Converts the sample to a string,
   *
   * @return String representation of the sample.
   */
  inline std::string ToString() const {
    std::stringstream ss;
    for (VariableIndex i = FIRST_VAR; i <= num_var(); i++) {
      AssignmentEncoding var_val = var_assignment(i);
      switch (var_val) {
        case ASSN_F: ss << '0'; break;
        case ASSN_T: ss << '1'; break;
        case ASSN_U: ss << '*'; break;
      }
    }
    return ss.str();
  }
  /**
   * Variable Assignment Accessor.
   *
   * Extracts the value of a single variable from the partial assignment
   *
   * @param var Identification number of the variable from 1 to num_var
   *
   * @return Value of the specified value number
   */
  inline const AssignmentEncoding var_assignment(const VariableIndex &var) const {
    VariableIndex word_num, bit_num;
    calculateWordAndBitNumbers(var, word_num, bit_num);
    return (AssignmentEncoding) ((assn_[word_num] >> bit_num) & var_mask_);
  }
  /**
   * Extracts the current assignment in vector form.
   *
   * @return Assignment as a vector
   */
  inline void GetPartialAssignment(PartialAssignment &all_vars) const {
    all_vars.clear();
    all_vars.resize(num_var_ + FIRST_VAR);
    for (VariableIndex i = FIRST_VAR; i <= num_var_; i++)
      all_vars[i] = var_assignment(i);
  }
  /**
   * Variable Value Setter
   *
   * Modifies the bit values for a specific variable.
   *
   * @param var Variable identification number between 1 and num_var
   * @param val Value to encode the variable as.
   * @param lit Literal ID
   */
  inline void setVarAssignment(const VariableIndex var, const AssignmentEncoding &val) {
    assert((val == ASSN_F || val == ASSN_T) && var_assignment(var) == ASSN_U);

    VariableIndex word_num, bit_num;
    calculateWordAndBitNumbers(var, word_num, bit_num);

    // Updates the bits of interest only by masking then setting them.
    assn_[word_num] = (assn_[word_num] & ~(var_mask_ << bit_num)) | (val << bit_num);
    num_vars_set_++;
    remaining_vars_.clear();
  }
  /**
   * Defines the number of variables in the Boolean formula.  This is used to
   * encode the recipe.
   *
   * The number of variables dictates the size of the
   *
   * @param num_var Number of variables in the equation.
   */
  static void set_num_var(const VariableIndex num_var) {
    assert(num_var_ == 0);  // This function should onyl be called once.
    num_var_ = num_var;
    var_vec_len_ = ((num_var + 1) * bits_per_var_ / size_of_word_)+ 1;
  }
  /**
   * Complete Assignment Checker
   *
   * Determines whether the sample model is partial or complete.
   *
   * @return true if the sample model is a complete assignment.
   */
  const bool IsComplete() const {
//    return cache_comp_ids().empty();
    return num_set_vars_const() + emancipated_vars_.size() == num_var_;
  }
//  /**
//   * Complete Assignment Checker
//   *
//   * Determines whether the sample model is partial or complete.
//   *
//   * @return true if the sample model is a complete assignment.
//   */
//  bool IsComplete() const {
//    return num_set_vars_const() + emancipated_vars_.size() == num_var_;
//  }
//  /**
//   * Accessor for the number of cached component IDs in this object.
//   * @return Number of cached component IDs in this object
//   */
//  inline const uint64_t cached_comp_count() const { return cache_comp_ids_.size(); }
  /**
   * Accessor for the number of samples in this object.
   *
   * @return Number of samples represented by the object.
   */
  const SampleSize sample_count() const { return sample_count_; }
  /**
   * Accesor for the number of variables set in this assignment.  Note that this does not include
   * any emancipated variables since those are free and not set.
   *
   * @return Number of variables set in this assignment.
   */
  const VariableIndex num_set_vars() {
    if (!remaining_vars_.empty())
      return num_vars_set_;

    // ToDo once variable setting is efficient, just return the set variable count
    num_vars_set_ = 0;
    for (VariableIndex i = FIRST_VAR; i <= num_var(); i++) {
      if (var_assignment(i) != ASSN_U) {
        num_vars_set_++;
      } else if (!IsVarEmancipated(i)) {
        remaining_vars_.emplace_back(i);
      }
    }
    return num_vars_set_;
  }
  /**
   * Accessor for the number of variables set in this assignment.  Note that this does not include
   * any emancipated variables since those are free and not set.
   *
   * @return Number of variables set in this assignment.
   */
  const VariableIndex num_set_vars_const() const {
    // ToDo once variable setting is efficient, just return the set variable count
    SampleSize num_vars_set = 0;
    for (VariableIndex i = FIRST_VAR; i <= num_var(); i++)
      if (var_assignment(i) != ASSN_U)
        num_vars_set++;
    return num_vars_set;
  }
//  /**
//   * Updates the assignment of the implicit assignment with that of the specified one.  It does
//   * NOT update the sample count
//   *
//   * @param other Another SampleAssignment.
//   */
//  void updateAssignmentOnly(const SampleAssignment &other) {
//    this->assn_ = other.assn_;
//    this->emancipated_vars_ = other.emancipated_vars_;
//  }
  /**
   * Gets the number of remaining variables that are unset.  This does NOT include emancipated
   * variables.
   *
   * @return Number of unset unemancipated variables.
   */
  const VariableIndex num_unset_vars() {
    return num_var_ - num_set_vars() - emancipated_vars_.size();
  }
//  /**
//   * Builds and returns the set of unconstrained variables in this sample.
//   *
//   * @return Identification number of the unset variables.
//   */
//  const std::vector<VariableIndex> GetUnsetConstrainedVars() const {
//    std::vector<VariableIndex> unset_vars;
//    for (VariableIndex i = FIRST_VAR; i <= num_var(); i++)
//      if (var_assignment(i) == ASSN_U && !IsVarEmancipated(i))
//        unset_vars.push_back(i);
//    return unset_vars;
//  }
  /**
   * Checks if the specified variable is emancipated.
   *
   * @param var Variable identification number
   * @return True if the variable is emancipated.
   */
  const bool IsVarEmancipated(VariableIndex var) const {
    return std::find(emancipated_vars_.begin(), emancipated_vars_.end(), var)
           != emancipated_vars_.end();
  }
  /**
   * Cached component ID accessor.
   *
   * @returns All cached components identification numbers in this object.
   */
  const std::vector<CacheEntryID>& cache_comp_ids() const {
    return cache_comp_ids_;
  }
//  /**
//   * Generate a unique key for the cached components in the sample assignment.
//   *
//   * @return Key string for the cached components.
//   */
//  std::string GetCachedCompKey() {
//    if (!comp_cache_key_.empty())
//      return comp_cache_key_;
//    std::sort(cache_comp_ids_.begin(), cache_comp_ids_.end());
//    std::stringstream ss;
//    for (auto cached_comp : cache_comp_ids_) {
//      if (cached_comp != cache_comp_ids_[0])
//        ss << ",";
//      ss << cached_comp;
//    }
//    comp_cache_key_ = ss.str();
//    return comp_cache_key_;
//  }
};


typedef std::list<SampleAssignment> ListOfSamples;

class SamplesManager {
 private:
  /**
   * Total number of samples observed since this SampleManager's creation.
   */
  mpz_class solution_count_ = 0;
  ListOfSamples samples_;
//  /**
//   * Stores the final expected number of samples to be returned at the
//   * end of sampling.
//   */
//  static SampleSize final_num_samples_;
//  /**
//   * Stores the current number of samples being built by the sampler.
//   */
//  static SampleSize samples_manager_vector_size_;
  /**
   * Number of samples stored by this object.
   */
  SampleSize tot_num_samples_;
  /**
   * Configuration of the solver that created this object.
   */
  SolverConfiguration *config_;
//  /**
//   * Used in the random number generator.  Stores the random bits
//   * to be extracted.
//   */
//  static int random_bits_;
//  /**
//   * Next random bit to be used. It is between [0,NUM_INT_BITS-1).
//   */
//  static int next_rand_bit_;
//  static const int NUM_INT_BITS_;
  /**
   * New Sample Builder
   *
   * Based off the current state of the solver, this function creates a set of
   * new samples and stores them in the "new_samples" object.
   *
   * @param new_samples New sample that will be created
   * @param active_comp Currently active component
   * @param literal_stack Stack of assigned literals
   * @param last_branch_lit Last literal used as a branch variable (i.e., not a BCP variable).
   *        This parameter determines where in the literal stack to search for
   *        already assigned variables.
   * @param ana Component analyzer for selecting solutions.
   * @param freed_vars Variables that became free (i.e., not in any clauses
   *        during assignment, and remain unset.
   */
  static void BuildSample(SampleAssignment &new_sample,
                          const Component *active_comp,
                          const std::vector<LiteralID> &literal_stack,
                          VariableIndex last_branch_lit,
                          const AltComponentAnalyzer &ana,
                          const std::vector<VariableIndex> &freed_vars,
                          const std::vector<CacheEntryID> &cached_comp_ids);
  /**
   * Literal Stack Contents Checker
   *
   * Checks whether the specified variable is in the literal stack.
   *
   * @param var_num Variable number of interest
   * @param literal_stack List of assigned literals
   * @param start_ofs Location in the literal stack to begin the search
   * @return True if the variable is in the literal stack.
   */
  static bool isVarInLitStack(const VariableIndex var_num,
                              const std::vector<LiteralID> &literal_stack,
                              const VariableIndex start_ofs = 0) {
    for (auto itr = literal_stack.begin() + start_ofs; itr != literal_stack.end(); itr++)
      if ((*itr).var() == var_num)
        return true;
    return false;
  }
  /**
   * Clause Literal List Builder
   *
   * Builds a list of clauses of literals.  Each clause is an element in the other
   * list.  A clause consists of a set of literals.  Each literal is the variable
   * number with the signing representing whether the literal is negated or
   * not.
   *
   * @param input_file_path Path to the CNF file
   */
  static void buildCnfClauseLiterals(const std::string &input_file_path,
                                     std::vector<std::vector<signed long>> &clauses);
  /**
   * Splits a sample object into two for use during stitching.
   *
   * @param itr Iterator pointing to the sample object to be split.  This will be changed to
   * point to the new inserted element.
   *
   * @param new_assn_sample_count Size of the new sample object BEFORE the split.  It should be
   * greater than zero and less than the number of samples in the object.
   */
  void splitSampleAndInsert(std::list<SampleAssignment>::iterator &itr,
                            const SampleSize &new_assn_sample_count) {
    assert(new_assn_sample_count > 0 && new_assn_sample_count < itr->sample_count());
    // Need to increment then decrement since this method inserts the new element before itr
    SampleAssignment new_node = itr->split(new_assn_sample_count);
    samples_.insert(itr, new_node);
    --itr;
  }

 public:
  /**
   * Initialize the sample manager.  It creates complete blank samples.
   */
  SamplesManager(SampleSize num_samples, SolverConfiguration &config) : config_(&config) {
    tot_num_samples_ = num_samples;
  }
//  /**
//   * Copy constructor.
//   */
//  SamplesManager(const SamplesManager &other)
//      : SamplesManager(other.num_samples(), *other.config_) {}
//  /**
//   * Equality operator.
//   *
//   * @param other Object to which the implicit object will be set.
//   * @return Reference to the new SamplesManager created.  This allows for chaining equality
//   * operators.
//   */
//  SamplesManager& operator=(const SamplesManager &other) {
//    this->tot_num_samples_ = other.tot_num_samples_;
//    this->samples_ = other.samples_;
//    this->config_ = other.config_;
//    return *this;
//  }
  /**
   * Sample Exporter
   *
   * Exports the set of samples to an output stream.  This can be either
   * a file or to the console.
   *
   * @param output_file_path Path for export file
   * @param statistics CNF and runtime statistics structure.
   * @param config Solver configuration
   */
  void exportFinal(std::ostream &out, const DataAndStatistics &statistics,
                   const SolverConfiguration& config);
  /**
   * Performs reservoir sampling
   *
   * @param active_comp Currently active component.
   * @param literal_stack Current assigned literal stack
   * @param sample_weight Weight of the current set of samples to assign.
   */
  void reservoirSample(const Component * active_comp,
                       const std::vector<LiteralID> & literal_stack,
                       const mpz_class &solution_weight,
                       const mpz_class &weight_multiplier,
                       const AltComponentAnalyzer &ana,
                       VariableIndex literal_stack_ofs,
                       const std::vector<VariableIndex>& freed_vars,
                       const std::vector<CacheEntryID> &cached_comp_ids,
                       const CachedAssignment& cached_assn,
                       SampleAssignment& cached_sample);
  /**
   * Sample Replacement List Generator
   *
   * Builds a list of the samples that will be replaced from the assignment
   * collection.  The selection of whether a sample will be replaced
   * is based off their relative weighting of their sample counts.
   *
   * @param new_sample_weight Model count weight for the new samples.
   * @param samples_to_replace Contains the indices of the samples in the list that will be
   * replaced.
   */
  inline void GenerateSamplesToReplace(const mpz_class &new_sample_weight,
                                       std::vector<SampleSize> &samples_to_replace) const;
//  /**
//   * Sample Variable Assignment Accessor
//   *
//   * Gets the value of the variable assignment for a specific sample in the list of samples
//   *
//   * Debug function.
//   *
//   * @param sample_num Sample number between 0 (inclusive) and num_samples (exclusive)
//   * @param var Variable number
//   *
//   * @return Variable's assigned value.
//   */
//  AssignmentEncoding sample_var_val(const SampleSize sample_num, const VariableIndex var) const {
//    assert(var >= FIRST_VAR && var <= samples_[sample_num].num_var());
//    return samples_[sample_num].var_assignment(var);
//  }
  /**
   * Samples Manager Stitcher
   *
   * Stitches together to sample managers.  This is used when combining
   * the results after a component split.
   *
   * @param other Component split to be merged with.
   */
  inline void stitch(SamplesManager &other) {
    assert(this->tot_num_samples_ == other.tot_num_samples_);
    if (solution_count_ == 0) {
      this->samples_ = other.samples_;
      solution_count_ = other.solution_count_;
      return;
    } else {
      solution_count_ *= other.solution_count_;
    }
    // Handle the UNSAT case
    if (other.solution_count_ == 0) {
      this->samples_.clear();
      return;
    }

    // Depending on the number of sample objects, different stitching approaches are faster
    StitchShuffledArray(other);

    // After the size normalization, perform the stitching sample by sample.
    assert(verifyPostStitchingCorrectness(other));
  }
  /**
   * Perform stitching permutation generation by creating an vector of size tot_num_samples_ (|S|)
   * and shuffling it via a Fisher-Yates shuffle which has running time Theta(|S|).
   *
   * @param other SamplesManager that will be stitched to the implicit SamplesManager.
   */
  inline void StitchShuffledArray(SamplesManager &other) {
    assert(this->GetActualSampleCount() == other.GetActualSampleCount());

    std::vector<ListOfSamples::iterator> other_samples_itrs;
    other_samples_itrs.reserve(other.samples_.size());
    std::vector<SampleSize> other_sample_order;
    other_sample_order.reserve(other.tot_num_samples_);

    // Store a reference to each element in other's samples
    // These will be used for simplifying the stitching look-up.
    SampleSize i = 0;
    for (auto other_itr = other.samples_.begin(); other_itr != other.samples_.end(); ++other_itr) {
      other_samples_itrs.emplace_back(other_itr);
      for (SampleSize j=0; j < other_itr->sample_count(); j++)
        other_sample_order.emplace_back(i);
      ++i;
    }
    Random::shuffle<SampleSize>(other_sample_order.begin(), other_sample_order.end());

    // Split and merge to create the permutations
    SampleSize sample_offset = 0;
    for (auto this_itr = samples_.begin(); this_itr != samples_.end(); ) {
      SampleSize sample_end = sample_offset + this_itr->sample_count();
//      // Sort the subindices for each "this" sample to make it easier to count and organize
//      if (this_itr->sample_count() > 1) {
//        auto itr_begin = other_sample_order.begin(), itr_end = other_sample_order.begin();
//        std::advance(itr_begin, sample_offset);
//        std::advance(itr_end, sample_end);
//        std::sort(itr_begin, itr_end);
//      }
      std::vector<SampleSize> samples_per_element(other_sample_order.size(), 0);
      for (SampleSize sample_cnt = sample_offset; sample_cnt < sample_end; ++sample_cnt)
        samples_per_element[other_sample_order[sample_cnt]]++;

      // Indices from the same "other" sample are grouped together so split and stitch.
      for (SampleSize other_cnt = 0; other_cnt < samples_per_element.size(); ++other_cnt) {
        if (samples_per_element[other_cnt] == 0)
          continue;
        SplitAndStitch(this_itr, other_samples_itrs[other_cnt], samples_per_element[other_cnt]);
      }

      // Update all pointer
      sample_offset = sample_end;
    }
  }
  /**
   * Splits the sample point to by \p this_itr.  The new split sample has \p num_new_samples. This
   * new sample is then stitched with the sample pointed to by \p other_itr.  After performing the
   * stitching, the sample count of \p other_itr is decreased by \p num_new_samples.
   *
   * **Note**: that if \p num_new_samples and the sample count of \p this_itr are equal, the sample
   * pointed to by this_itr is not split.  All other steps proceed normally.
   *
   * @param this_itr Sample to be split and then merged.
   *
   * @param other_itr Object that will be stitched with the new split sample.
   *
   * @param num_new_samples Sample size for the new stitched object.  It must be greater than zero
   * and less than or equal to the sample count of \p this_itr.
   */
  void SplitAndStitch(ListOfSamples::iterator &this_itr, ListOfSamples::iterator &other_itr,
                      SampleSize &num_new_samples) {
    assert(num_new_samples > 0 && num_new_samples <= this_itr->sample_count()
           && num_new_samples <= other_itr->sample_count());

    // If applicable add a new node
    if (num_new_samples != this_itr->sample_count())
      splitSampleAndInsert(this_itr, num_new_samples);

    this_itr->stitch(*other_itr);
    ++this_itr;
    if (other_itr->sample_count() == num_new_samples)
      other_itr->zeroSampleCount();
    else
      other_itr->DecreaseSampleCount(num_new_samples);
  }
  /**
   * Debug helper function to verify that post stitching, the number of samples is correct.
   *
   * @param other Set of samples that there were stitched to the implicit sample set.
   *
   * @return True if the stitching appears correct.
   */
  const bool verifyPostStitchingCorrectness(SamplesManager &other) const {
    for (auto &sample : other.samples_) {
      if (sample.sample_count() != 0) {
        PrintInColor("An other_sample had non-zero size.", PrintColor::COLOR_RED);
        return false;
      }
    }

    if (!this->verifySampleCount()) {
      PrintInColor("The sample count of the stitched object is incorrect.", PrintColor::COLOR_RED);
      return false;
    }

    return true;
  }
  /**
   * Samples Manager Merger
   *
   * After a component split, a new descendant SamplesManager is created.  This function
   * will merge the existing SamplesManager with the original one before the component
   * split.
   *
   * @param other Sample manager created for the descendants of the component split.
   *
   * @param other_multiplier Used to scale the model count of the @see other
   * model count.
   *
   * @param freed_vars List of emanicipated variables that can be set to either
   * true or false.
   *
   * @param cached_assn Assignment stored in the cached that is used when
   * storing samples in cache.
   *
   * @param cached_sample Sample output when performing sample caching.  This will
   * be used to update the top of the decision stack.
   */
  void merge(SamplesManager &other, const mpz_class &other_multiplier,
             const std::vector<VariableIndex> &freed_vars,
             const std::vector<CacheEntryID> &cached_comp_ids,
             const CachedAssignment & cached_assn,
             SampleAssignment& cached_sample);
  /**
   * Sample Verifier
   *
   * Verifies that all samples generated by the program are valid and actually
   * satisfies the input file.  This is mostly for debug and should NEVER
   * return false.  If false is ever returned, something is wrong.
   *
   * @param input_file_path Path to a file in CNF formula
   * @param skip_unassigned If true, skip verification of all clauses with
   *        an unassigned literal.
   * @return true if all samples are valid.
   */
  bool VerifySolutions(const std::string &input_file_path,
                       bool skip_unassigned = false) const;
  /**
   * Model Count Extractor
   *
   * Gets the current number of solutions tracked by the manager object
   *
   * @return Model count for the samples manager object
   */
  const mpz_class &model_count() const { return solution_count_; }
  /**
   * Emancipated and Unused Variable Adder
   *
   * Unused variables are any variables that do not appear at all in a given
   * formula.  This sets them at the end of the countSAT execution.
   *
   * @param emancipated_vars List of the IDs of unused variables.
   */
  void AddEmancipatedVars(const std::vector<VariableIndex> &emancipated_vars) {
    for (auto &sample : samples_)
      sample.addEmancipatedVars(emancipated_vars);

    // Multiply by 2^num_unused_vars since those represent a parallel cylinder
    mpz_mul_2exp(solution_count_.get_mpz_t(), solution_count_.get_mpz_t(),
                 emancipated_vars.size());
  }
  /**
   * Adds the specified cached component identification numbers to all samples in the
   * collection.
   *
   * @param cached_comp_ids Cached component identification numbers for all samples.
   */
  void AddCachedCompIds(const std::vector<CacheEntryID> &cached_comp_ids) {
    for (auto &sample : samples_)
      sample.addCachedCompIds(cached_comp_ids);
  }
  void TransferVariableAssignments(ListOfSamples &others) {
    std::vector<VariableIndex> unset_vars = others.front().GetRemainingVariables();

    assert(!unset_vars.empty());
    auto unsetter = SampleAssignment::buildUnsetterAssignment(unset_vars);
    // Delete redundant emancipated variables.
    for (auto &sample : samples_) {
      sample.unsetVariableAssignments(unsetter);
      sample.DeleteEmancipatedVars(others.front().emancipated_vars());
    }

    SampleSize sample_count = 0;
    for (auto &other : others)
      sample_count += other.sample_count();
    assert(sample_count == this->num_samples());
    SamplesManager others_manager(sample_count, *config_);
    others_manager.solution_count_ = 1;
    others_manager.samples_ = others;

    this->stitch(others_manager);
  }
//  /**
//   * Solver Configuration Storer
//   *
//   * This function is used to store a reference to the solver's configuration.
//   * This is useful in case any of the run parameters are used.
//   * @param config Solver's configuration.
//   */
//  static void set_solver_config(SolverConfiguration &config) { config_ = &config; }
  /**
   * Accessor to a sample manager's variable assignment.
   *
   * @return Reference to the sample managers recipe assignments.
   */
  ListOfSamples &samples() { return samples_; }
//  /**
//   * Sample Setter
//   *
//   * This function replaces the sample (based off the specified number)
//   * with the new model passed to the function.
//   *
//   * @param sample_num Sample number to be set
//   * @param new_model New sample model to be stored
//   */
//  inline void set_sample(SampleSize sample_num,
//                         const SampleAssignment &new_model) {
//    assert(sample_num >= 0 && sample_num < samples_.size());
//    samples_[sample_num] = new_model;
//  }
//  /**
//   * Sample Accessor
//   *
//   +   * Extracts a reference to the specified sample from the manager.
//   *
//   +   * @param sample_num Number of the sample to access - base 0
//   *
//   +   * @return Sample at the specified number
//   */
//  inline const SampleAssignment &sample(SampleSize sample_num) const {
//    assert(sample_num >= 0 && sample_num < samples_.size());
//    return samples_[sample_num];
//  }
  /**
   * Checks whether all the samples in the set are complete.
   *
   * @return true if all contained samples are complete.
   */
  inline bool IsComplete() const {
    for (const auto &sample : samples_)
      if (!sample.IsComplete())
        return false;
    return true;
  }
  /**
   * Adds the samples from one SamplesManager() to that of annother.  This is different from merging
   * and stitching as the implicit SamplesManager()'s samples remain unchanged.  The only difference
   * is the samples from \p other are appended to the collections of samples in the implicit
   * parameter.
   *
   * @param other Another SampleManager() object whose samples will be copied to the implicit
   * parameter.
   */
  inline void append(SamplesManager &other) {
    append(other.samples_);
  }
  /**
   * Adds the samples from one SamplesManager() to that of annother.  This is different from merging
   * and stitching as the implicit SamplesManager()'s samples remain unchanged.  The only difference
   * is the samples from \p other are appended to the collections of samples in the implicit
   * parameter.
   *
   * @param other A list of samples to add.
   */
  inline void append(ListOfSamples &other) {
    samples_.splice(samples_.end(), other);
    assert(GetActualSampleCount() <= tot_num_samples_);
  }
  /**
   * Accessor for the number of samples actually stored by the SamplesManager() object.
   *
   * @return Number of samples actually stored by the SamplesManager().
   */
  inline const SampleSize num_samples() const {
//    assert(verifySampleCount());
    return tot_num_samples_;
  }
  /**
   * Removes the samples specified at the indices in the \p samples_to_remove input.
   *
   * @param samples_to_remove Identification numbers of the samples to remove.
   */
  void RemoveSamples(std::vector<SampleSize> &samples_to_remove) {
    if (samples_to_remove.empty())
      return;
    if (num_samples() == samples_to_remove.size()) {
      samples_.clear();
      return;
    }

    // Delete from back to front to prevent deletion affecting counts.
    assert(!samples_.empty());
    auto sample_itr = --(samples_.end());
    SampleSize cur_sample_start_count = num_samples() - sample_itr->sample_count();
    SampleSize num_to_remove = 0;
    for (auto &sample_to_remove : samples_to_remove) {
      assert(sample_to_remove >= 0 && sample_to_remove < num_samples());
      // Skip to the next sample to remove.
      while (sample_to_remove < cur_sample_start_count) {
        if (num_to_remove > 0) {
          if (sample_itr->sample_count() == num_to_remove)
            // If no samples are left then remove the object
            sample_itr = samples_.erase(sample_itr);
          else
            sample_itr->DecreaseSampleCount(num_to_remove);
          num_to_remove = 0;
        }
        --sample_itr;
        cur_sample_start_count -= sample_itr->sample_count();
      }
      num_to_remove++;
    }

    // Handle the last element that requires removal
    if (sample_itr->sample_count() == num_to_remove)
      sample_itr = samples_.erase(sample_itr);
    else
      sample_itr->DecreaseSampleCount(num_to_remove);
    assert(GetActualSampleCount() == tot_num_samples_ - samples_to_remove.size());
  }
  /**
   * Specifies which samples from this SampleManager should be kept.  All others are discarded.
   * This is used as a compliment to the RemoveSamples() function since when there is a merge,
   * one set keeps samples and the other discards them.
   *
   * @param samples_to_keep Identification number of the samples that will be kept.
   */
  void KeepSamples(std::vector<SampleSize> &samples_to_keep) {
    if (samples_to_keep.empty()) {
      samples_.clear();
      return;
    }
    if (num_samples() == samples_to_keep.size())
      return;

    // Delete from back to front to prevent deletion affecting counts.
    assert(!samples_.empty());
    auto sample_itr = --(samples_.end());
    SampleSize cur_sample_start_count = num_samples() - sample_itr->sample_count();
    SampleSize num_to_keep = 0;
    for (auto &sample_to_keep : samples_to_keep) {
      assert(sample_to_keep >= 0 && sample_to_keep < num_samples());
      // Skip to the next sample to check
      while (sample_to_keep < cur_sample_start_count) {
        if (num_to_keep > 0) {
          sample_itr->set_sample_count(num_to_keep);
          num_to_keep = 0;
        } else {
          sample_itr = samples_.erase(sample_itr);
        }
        --sample_itr;
        cur_sample_start_count -= sample_itr->sample_count();
      }
      num_to_keep++;
    }

    // Handle the last node to keep
    sample_itr->set_sample_count(num_to_keep);
    // Any nodes never encountered have to be removed.
    if (sample_itr != samples_.begin())
      samples_.erase(samples_.begin(), sample_itr);
    assert(GetActualSampleCount() == samples_to_keep.size());
  }
  /**
   * Verifies that the SampleManager actually has the number of samples it is supposed to.
   *
   * @return True if the sample count is correct.
   */
  const bool verifySampleCount() const {
    return GetActualSampleCount() == tot_num_samples_;
  }
  /**
   * Gets the actual number of samples that are contained in this manager. In MOST (but not all)
   * cases, this should equal the expected number of samples.  However, that is not necessarily
   * the case in some rare but useful exceptions.
   *
   * @return Actual number of samples contained in this object.
   */
  const SampleSize GetActualSampleCount() const {
    SampleSize actual_sample_count = 0;
    for (auto &sample : samples_) {
      assert(sample.sample_count() > 0);  // Make sure no dead samples
      actual_sample_count += sample.sample_count();
    }
    return actual_sample_count;
  }
  /**
   * Eliminate all existing samples.
   */
  void clear() { samples_.clear(); }
};


///**
// * Recipe Structure Initializer
// *
// * Initializes recipe static structures.  This requires a dedicated function
// * because of the scope requirements of static objects.
// *
// * @param num_var Number of variables in the Boolean formula.
// */
//void InitializeSamplerStructures(VariableIndex num_var);
///**
// * Sample Size Initializer
// *
// * Initializes the number of samples to by collected by the algorithm.
// *
// * @param sample_count Number of samples
// */
//void InitializeSampleCount(SampleSize sample_count);
//
//
//bool IsVarInLiteralStack(const std::vector<LiteralID> &literal_stack, VariableIndex var);

#endif //SHARPSAT_SOLUTION_RECIPE_H
