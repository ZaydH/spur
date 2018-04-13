/**
 * component_management.cpp
 *
 * Purpose: Defines methods for the ComponentManagement() class.
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


#include <algorithm>
#include "component_management.h"

void ComponentManager::initialize(LiteralIndexedVector<Literal> & literals,
                                  std::vector<LiteralID> &lit_pool, bool quiet) {
  ana_.initialize(literals, lit_pool);
  // BEGIN CACHE INIT
  CacheableComponent::adjustPackSize(ana_.max_variable_id(), ana_.max_clause_id());
  // Prevent a memory leak
  for (auto &i : component_stack_)
    delete i;
  component_stack_.clear();
  component_stack_.reserve(ana_.max_variable_id() + 2);
  component_stack_.push_back(new Component());
  component_stack_.push_back(new Component());
  assert(component_stack_.size() == 2);
  component_stack_.back()->createAsDummyComponent(ana_.max_variable_id(),
      ana_.max_clause_id());
  if (config_.perform_random_sampling_)
    cached_vars_.clear();

  cache_.init(*component_stack_.back(), quiet);
}


void ComponentManager::removeAllCachePollutionsOf(StackLevel &top) {
  // all processed components are found in
  // [top.currentRemainingComponent(), component_stack_.size())
  // first, remove the list of descendants from the father
  assert(top.remaining_components_ofs() <= component_stack_.size());
  assert(top.super_component() != 0);
  assert(cache_.hasEntry(super_component(top).id()));

  if (top.remaining_components_ofs() == component_stack_.size())
    return;

  for (unsigned u = top.remaining_components_ofs(); u < component_stack_.size();
      u++) {
    assert(cache_.hasEntry(component_stack_[u]->id()));
    cache_.cleanPollutionsInvolving(component_stack_[u]->id());
  }

#ifdef DEBUG
  cache_.test_descendantstree_consistency();
#endif
}

std::vector<VariableIndex> ComponentManager::buildFreedVariableList(const StackLevel & top,
                                                     const std::vector<LiteralID> &literal_stack) {
  Component * descendant_comp, *ref_comp = component_stack_[top.super_component()];
  const VariableIndex original_var_count = ref_comp->num_variables();

  VariableIndex num_unused_vars = original_var_count;
  std::vector<bool> var_in_descendants(original_var_count, false);
  std::vector<ComponentVarAndCls> descendant_vars, ref_vars = ref_comp->getData();

  for (unsigned comp_id=top.remaining_components_ofs();
       comp_id < top.unprocessed_components_end(); comp_id++) {
    descendant_comp = component_stack_[comp_id];
    descendant_vars = descendant_comp->getData();
    UpdateVarDescedantsList(ref_vars, ref_comp->num_variables(), var_in_descendants,
                            descendant_vars, descendant_comp->num_variables());
    num_unused_vars -= descendant_comp->num_variables();
  }

  // Build a list of all literals assigned in the last round
  const long list_stack_ofs = top.literal_stack_ofs();
  VariableIndex num_lits_assigned = literal_stack.size() - list_stack_ofs;

  std::vector<ComponentVarAndCls> assigned_stack_vars;
  assigned_stack_vars.reserve(num_lits_assigned);
  for (VariableIndex stack_i = 0; stack_i < num_lits_assigned; stack_i++)
    assigned_stack_vars.push_back(literal_stack[list_stack_ofs + stack_i].var());
  // Sort the vector and exclude the assigned literals
  sort(assigned_stack_vars.begin(), assigned_stack_vars.end());
  num_unused_vars -= UpdateVarDescedantsList(ref_vars, ref_comp->num_variables(),
                                             var_in_descendants, assigned_stack_vars,
                                             num_lits_assigned);

  // Skip the cached variables
  num_unused_vars -= UpdateVarDescedantsList(ref_vars, ref_comp->num_variables(),
                                             var_in_descendants,
                                             cached_vars_, cached_vars_.size());

  // If a variable does not appear in the descendant nor in the cache and is not on
  // top of the literal stack, it was freed.
  std::vector<VariableIndex> freed_vars;
  freed_vars.reserve(num_unused_vars);
  for (VariableIndex ref_i = 0; ref_i < ref_comp->num_variables(); ref_i++) {
    if (!var_in_descendants[ref_i])
      freed_vars.push_back(ref_vars[ref_i]);
  }
  assert(num_unused_vars == freed_vars.size());
  return freed_vars;
}


VariableIndex UpdateVarDescedantsList(const std::vector<ComponentVarAndCls> &ref_vars,
                                      VariableIndex ref_num_vars,
                                      std::vector<bool> &varInDescendants,
                                      std::vector<ComponentVarAndCls> descendant_vars,
                                      VariableIndex desc_num_vars) {
  if (desc_num_vars == 0)
    return 0;
  VariableIndex varsRemoved = 0;
  VariableIndex ref_i, desc_i = 0;
  for (ref_i = 0; ref_i < ref_num_vars; ref_i++) {
    // Reference list is always a superset and since both are sorted, go to the next
    if (ref_vars[ref_i] < descendant_vars[desc_i])
      continue;
    if (ref_vars[ref_i] > descendant_vars[desc_i])
      desc_i++;
    if (desc_i >= desc_num_vars)
      break;

    // If the descendant and reference match, the var is present
    if (ref_vars[ref_i] == descendant_vars[desc_i]) {
      // DebugZH - A variable should only be in a single descedent
      assert(!varInDescendants[ref_i]);
      varInDescendants[ref_i] = true;
      varsRemoved++;
    }
  }
  return varsRemoved;
}
