/**
 *  cached_assignment.h
 *
 *  Purpose: Stores a variable assignment retrieved from the cache.
 *
 *  @author Zayd Hammoudeh <zayd@ucsc.edu>
 *  @version 0.00.00
 *
 *  Copyright (C) 2018 Zayd Hammoudeh.
 *  All rights reserved.
 *
 * This software may be modified and distributed under the terms of the MIT license.  See the
 * LICENSE file for details.
 *
 * Original Author: Marc Thurley.
 */

#ifndef SHARPSAT_CACHED_ASSIGNMENT_H
#define SHARPSAT_CACHED_ASSIGNMENT_H

#include <algorithm>
#include <vector>
#include "structures.h"

#define CACHED_VARIABLE_LEN 2 // In Bits

class CachedAssignment{
 public:
  CachedAssignment() {
    assigned_literals_.reserve(1);
    emancipated_vars_.reserve(1);
  }
  /**
   * Sort the list of literals for easier management.
   */
  inline void Sort() {
    std::sort(assigned_literals_.begin(), assigned_literals_.end());
    std::sort(emancipated_vars_.begin(), emancipated_vars_.begin());
  }
  /**
   * Increases the size of the literals list by the specified amount.
   *
   * @param size_adder The amount the literals list will be increased.
   */
  inline void IncreaseSize(VariableIndex size_adder) {
    VariableIndex new_size = assigned_literals_.capacity() + size_adder;
    assigned_literals_.reserve(new_size);
  }
//  /**
//   * Add the specified literal to the list of literals in the cache
//   * assignment.
//   *
//   * @param lit New literal value.
//   */
//  inline void AddLiteral(LiteralID lit) {assigned_literals_.push_back(lit);}
  /**
   * Variable List Extractor
   *
   * Builds a list of variables in the cached assignment.
   *
   * @return List of variable numbers in this assignment.
   */
  const std::vector<VariableIndex> vars() {
    std::vector<VariableIndex> cached_vars;
    cached_vars.reserve(assigned_literals_.size());

    for (auto lit : assigned_literals_)
      cached_vars.push_back(lit.var());
    cached_vars.insert(cached_vars.end(), emancipated_vars_.begin(), emancipated_vars_.end());
    return cached_vars;
  }
  /**
   * Stores the total number of components that make up
   * this cached assignment.
   *
   * @return Number of different components that make up this
   * assignment.
   */
  const VariableIndex num_components() const { return num_cached_components; }
  /**
   * An empty assignment has no associated components.
   *
   * @return true if the assignment is empty.
   */
  const bool empty() const {return num_cached_components == 0; }
  /**
   * Processes a cached component and incorporates its information into
   * the cached assignment.
   *
   * @param comp Component to be processed
   *
   * @param model_count_and_assn Component unshifted model count.  It contains the
   * component assignment in the lower <i>n</i> bits where <i>n</i> is the
   * number of bits in the component.
   */
  void ProcessComponent(const Component * comp, mpz_class model_count_and_assn) {
    // Nothing to store in the UNSAT case
    if (model_count_and_assn == 0)
      return;
    IncreaseSize(comp->num_variables());  // Increase the capacity

    // Extract the literal assignments from the MPZ object
    std::vector<ComponentVarAndCls> comp_data = comp->getData();
    for (VariableIndex i = 0; i < comp->num_variables(); i++) {
      VariableIndex bit_index = CACHED_VARIABLE_LEN * i;
      // Check if the variable is emancipated
      if (mpz_tstbit(model_count_and_assn.get_mpz_t(), bit_index + 1) == 1) {
        emancipated_vars_.emplace_back(comp_data[i]);
      } else {
        bool sign = (mpz_tstbit(model_count_and_assn.get_mpz_t(), bit_index) == 1);
        assigned_literals_.emplace_back(LiteralID(comp_data[i], sign));
      }
    }
    num_cached_components++;
  }
  /**
   * Accessor for the literals in the cached assignment.
   *
   * @return List of literals in the cached assignment.
   */
  const std::vector<LiteralID>& literals() const { return assigned_literals_; }
  /**
   * Accessor for get the set of emancipated variabes in this cached assignment.
   *
   * @return Emancipated variables in the cached assignment.
   */
  const std::vector<VariableIndex>& emancipated_vars() const { return emancipated_vars_; }
  /**
   * Resets the component and deletes all associated component
   * assignment information including the number of components and
   * the assigned literals.
   */
  void clear() {
    assigned_literals_.clear();
    emancipated_vars_.clear();
    num_cached_components = 0;
  }
//  /**
//   * Component Split Accessor
//   *
//   * Accessor for whether this cached assignment corresponds with a component
//   * split.
//   *
//   * @return true if this cached assignment corresponds to a component split.
//   */
//  const bool is_component_split() const { return is_component_split_;}

 private:
  std::vector<VariableIndex> emancipated_vars_;
  std::vector<LiteralID> assigned_literals_;
  unsigned long num_cached_components = 0;
};

#endif // SHARPSAT_CACHED_ASSIGNMENT_H
