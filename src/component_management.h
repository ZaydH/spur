/**
 * component_management.h
 *
 * Purpose: Defines the ComponentManager() class for managing the Boolean formula components.
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

#ifndef COMPONENT_MANAGEMENT_H_
#define COMPONENT_MANAGEMENT_H_

#include <gmpxx.h>

#include <utility>
#include <algorithm>
#include <string>
#include <vector>

#include "component_types/component.h"
#include "component_cache.h"
#include "alt_component_analyzer.h"
#include "containers.h"
#include "stack.h"
#include "solver_config.h"
#include "model_sampler.h"
#include "top_tree_sampler.h"

typedef AltComponentAnalyzer ComponentAnalyzer;

class ComponentManager {
 public:
  ComponentManager(SolverConfiguration &config, DataAndStatistics &statistics,
                   LiteralIndexedVector<TriValue> & lit_values) :
    config_(config), /*statistics_(statistics),*/ cache_(statistics, config),
    ana_(statistics, lit_values) {
  }
  /**
   * Initializes the cache, component stack, and component analyzer.
   *
   * @param literals Set of valid literals
   * @param lit_pool
   * @param quiet Quiet parameter used to set whether console printing is enabled.
   */
  void initialize(LiteralIndexedVector<Literal> & literals,
                  std::vector<LiteralID> &lit_pool, bool quiet = false);

  unsigned scoreOf(VariableIndex v) {
    return ana_.scoreOf(v);
  }
  /**
   * Cache Model Count Storer
   *
   * If caching is enabled, this method stores the model count of the component
   * onto the stack.
   *
   * @param stack_comp_id Component ID of the stack element
   * @param value Model count of the element
   */
  void cacheModelCountOf(VariableIndex stack_comp_id, const mpz_class &value) {
    if (config_.perform_component_caching)
      cache_.storeValueOf(component_stack_[stack_comp_id]->id(), value);
  }
  /**
   * Cache Model Count and Assignment Encoder and Storer
   *
   * Prepend the component's MPZ value with the sample assignment.
   *
   * @param stack_comp_id Component ID of the stack element.
   * @param value Model count
   * @param assn Sample assignment to encode in the cache with the model count.
   * @param component Actual component to get the variables.
   */
  void cacheModelCountAndAssignment(const unsigned stack_comp_id, const mpz_class &value,
                                    const SampleAssignment &assn, const Component &component) {
    mpz_class combined_value = 0;
    mpz_class temp_value = 0;
    // If the component is UNSAT do not push any assignment
    if (value > 0) {
      VariableIndex last_shift_len = 0, offset_amount = 0;
      uint64_t cached_word = 0;
      const VariableIndex MAX_OFFSET_AMOUNT = 28;

      std::vector<ComponentVarAndCls> comp_data = component.getData();
      for (VariableIndex i = 0; i < component.num_variables(); i++) {
        AssignmentEncoding var_assn = assn.var_assignment(comp_data[i]);
        // Make sure the assignment is valid
        assert((var_assn != ASSN_U && !assn.IsVarEmancipated(comp_data[i]))
               || (var_assn == ASSN_U && assn.IsVarEmancipated(comp_data[i])));

        cached_word += var_assn << (offset_amount);
        if (offset_amount == MAX_OFFSET_AMOUNT) {
          mpz_import(temp_value.get_mpz_t(), 1, 1, sizeof(cached_word), 0, 0, &cached_word);
          combined_value += temp_value << last_shift_len;
          last_shift_len = CACHED_VARIABLE_LEN * (i + 1);
          offset_amount = 0;
          cached_word = 0;
        } else {
          offset_amount += CACHED_VARIABLE_LEN;
        }
      }
      if (offset_amount > 0) {
        mpz_import(temp_value.get_mpz_t(), 1, 1, sizeof(cached_word), 0, 0, &cached_word);
        combined_value += temp_value << last_shift_len;
      }
      combined_value += value << (CACHED_VARIABLE_LEN * component.num_variables());
    }
    cacheModelCountOf(stack_comp_id, combined_value);
  }
  /**
   * Gets a pointer to the actual super component that corresponds to
   * passed stack level.
   *
   * @param lev Level in the stack whose super component will be extracted
   *
   * @return Pointer to the super component of the stack item.
   */
  Component & superComponentOf(StackLevel &lev) {
    assert(component_stack_.size() > lev.super_component());
    return *component_stack_[lev.super_component()];
  }
  /**
   * Const version of the superComponentOf() method.
   *
   * @param lev Level in the stack whose super component will be extracted
   *
   * @return Pointer to the super component of the stack item.
   */
  const Component & super_component(StackLevel &lev) const {
    assert(component_stack_.size() > lev.super_component());
    return *component_stack_[lev.super_component()];
  }
  /**
   * Extracts the number of components currently in the component stack.
   *
   * @return Number of element in teh component stack.
   */
  VariableIndex component_stack_size() const {
    return static_cast<VariableIndex>(component_stack_.size());
  }

  void cleanRemainingComponentsOf(StackLevel &top) {
    while (component_stack_.size() > top.remaining_components_ofs()) {
      if (cache_.hasEntry(component_stack_.back()->id()))
        cache_.entry(component_stack_.back()->id()).set_deletable();
      delete component_stack_.back();
      component_stack_.pop_back();
    }
    assert(top.remaining_components_ofs() <= component_stack_.size());
  }
  /**
   *
   * @param comp_idx
   * @return
   */
  inline const Component &component(unsigned long comp_idx) const {
    return *component_stack_[comp_idx];
  }
//  Component & currentRemainingComponentOf(StackLevel &top) {
//    assert(component_stack_.size() > top.currentRemainingComponent());
//    return *component_stack_[top.currentRemainingComponent()];
//  }
  /**
   * Checks for the next yet to explore remaining component of top
   * returns true if a non-trivial non-cached component
   * has been found and is now stack_.TOS_NextComp()
   * returns false if all components have been processed;
   *
   * This function has been modified to support sampling.  That includes
   * passing the parameter "literal_stack".  This is neeeded to
   * build samples.
   *
   * @param top Top of the decision stack.
   * @param literal_stack_ Current assigned literal stack.
   * @param samples Current reservoir sample set.
   * @param is_valid_top_tree_store_point True if it is a currently valid point to store a top
   * tree element.
   * @return true if unprocessed components remain.
   */
  inline bool findNextRemainingComponentOf(StackLevel &top,
                                           const std::vector<LiteralID> &literal_stack,
                                           SamplesManager &samples,
                                           const bool &is_valid_top_tree_store_point);
  /**
   * Examines the top of the decision stack.  The function will
   * try to decompose the component if applicable.  For any new component,
   * the function will check if the new components exists in the cache;
   * if it does, the results are used to update the super component.
   *
   * Otherwise, the new component is added to the component stack (not
   * the decision stack).  The end of unprocessed components is also updated.
   *
   * @param top Element at the top of the decision stack.
   * @param literal_stack Current assigned literal stack.
   */
  inline void recordRemainingCompsFor(StackLevel &top,
                                      const std::vector<LiteralID> &literal_stack,
                                      const bool &is_valid_top_tree_store_point);
  /**
   * Build residual formula component.  It uses the contents of the literal stack
   * to decide how to construct remaining components.
   *
   * @param top Element at the top of the decision stack.
   */
  inline void buildResidualComponent(StackLevel &top);
  /**
   * Accessor for the top of the component stack.
   *
   * @return Pointer to the component on top of the component stack.
   */
  inline Component* top_stack() {
    return component_stack_.back();
  }

  inline void sortComponentStackRange(unsigned long start, unsigned long end);

  void gatherStatistics() {
//     statistics_.cache_bytes_memory_usage_ =
//       cache_.recompute_bytes_memory_usage();
    cache_.compute_byte_size_infrastructure();
  }

  void removeAllCachePollutionsOf(StackLevel &top);
  /**
   * Freed Variable List Builder
   *
   * When assigning variables in a formula it is common that
   * some variables become free (i.e., no longer appear in any of
   * the clauses).  These free variables will not appear on the literal
   * stack but must be tracked to create the sample.  This function stores all
   * such freed variables.
   *
   * @param literal_stack Current literal stack
   *
   * @return Vector of the variable indexes of all freed variables
   */
  std::vector<VariableIndex> buildFreedVariableList(const StackLevel & top,
                                                    const std::vector<LiteralID> &literal_stack);

 private:
  std::vector<VariableIndex> cached_vars_;
  std::vector<CacheEntryID> cached_comp_ids_;

  SolverConfiguration &config_;

//  DataAndStatistics &statistics_;

  std::vector<Component *> component_stack_;
  ComponentCache cache_;
  ComponentAnalyzer ana_;
};

/**
 * Sorts the components stack from start (inclusive) to end (exclusive).
 * The components with more variables are placed lower on the stack (i.e.,
 * further from the top) and processed later.
 *
 * @param start Start of component stack to be sorted (inclusive)
 * @param end End of the component stack to sort (exclusive)
 */
void ComponentManager::sortComponentStackRange(unsigned long start, unsigned long end) {
  assert(start <= end);
  // sort the remaining components for processing
  for (unsigned long i = start; i < end; i++)
    for (unsigned long j = i + 1; j < end; j++) {
      if (component_stack_[i]->num_variables()
          < component_stack_[j]->num_variables())
        std::swap(component_stack_[i], component_stack_[j]);
    }
}

bool ComponentManager::findNextRemainingComponentOf(StackLevel &top,
                                                    const std::vector<LiteralID> &literal_stack,
                                                    SamplesManager &samples,
                                                    const bool &is_valid_top_tree_store_point) {
  // record Remaining Components if there are none!
  if (component_stack_.size() <= top.remaining_components_ofs())
    recordRemainingCompsFor(top, literal_stack, is_valid_top_tree_store_point);

  assert(!top.branch_found_unsat());
  if (top.hasUnprocessedComponents()) {
    // Perform reservoir sampling if the model count is greater than zero
    if (config_.perform_random_sampling_ && top.getActiveModelCount() > 0) {
      if (!top.isComponentSplit() || top.isFirstComponent()) {
        top.includeSolutionSampleMultiplier(top.getActiveModelCount());
        std::vector<VariableIndex> freed_vars = buildFreedVariableList(top, literal_stack);
        top.addFreeVariables(freed_vars);
        top.addCachedCompIds(cached_comp_ids_);
      }
    }
    return true;
  }

  // if no component remains
  // make sure, at least that the current branch is considered SAT
  top.includeSolution(1);

  // Perform reservoir sampling if the model count is greater than zero
  if (config_.store_sampled_models()) {
    if (config_.verbose) {
      std::cout << "\t# Solutions Found: " << top.getActiveModelCount()
                << ". Total multiplied weight is "
                << top.getActiveModelCount() * top.getSamplerSolutionMultiplier()
                << std::endl;
    }
    SampleAssignment cached_sample(config_.num_samples_to_cache_);
    if (!top.cached_comp_ids().empty())
      cached_comp_ids_.insert(cached_comp_ids_.end(), top.cached_comp_ids().begin(),
                              top.cached_comp_ids().end());
    samples.reservoirSample(component_stack_[top.super_component()], literal_stack,
                            top.getActiveModelCount(), top.getSamplerSolutionMultiplier(),
                            ana_, top.literal_stack_ofs(), top.emancipated_vars(),
                            cached_comp_ids_,  top.cached_assn(), cached_sample);
    if (config_.perform_sample_caching()) {
      top.ClearCachedAssn();
      top.set_cache_sample(cached_sample);
    }
    // DebugZH
//    assert(samples.VerifySolutions(statistics_.input_file_, true));
  } else if (is_valid_top_tree_store_point) {
    assert(config_.perform_top_tree_sampling);
    // Store the top tree sample if a cylinder and no component split and above max depth.
    TopTreeSampler::StoreSample(top.getActiveModelCount(),
                                top.getSamplerSolutionMultiplier(), literal_stack,
                                top.emancipated_vars(), TopTreeNodeType::CYLINDER);
  }
  return false;
}
/**
 * Component at the top of the stack is analyzed.  If a new and unseen component
 * has never been seen, then the new component is placed on the component stack.
 *
 * The component stack is then sorted so components with fewer variables are
 * at the top of the stack.
 *
 * @param top Element at the top of the decision stack.
 * @param literal_stack_ Current assigned literal stack.
 * @param is_valid_top_tree_store_point True if it is a currently valid point to store a top
 * tree element.
 */
void ComponentManager::recordRemainingCompsFor(StackLevel &top,
                                               const std::vector<LiteralID> &literal_stack,
                                               const bool &is_valid_top_tree_store_point) {
  Component & super_comp = superComponentOf(top);
  VariableIndex new_comps_start_ofs = component_stack_.size();  // Location to store new comp if any

  // Clear the variables considered cached.
  cached_vars_.clear();
  cached_comp_ids_.clear();
  CachedAssignment cached_assn;

  // Initialize data structures for the component analyzer
  ana_.setupAnalysisContext(top, super_comp);

  // Go through each variable and checks its component to see if cached or new
  for (auto vt = super_comp.varsBegin(); *vt != varsSENTINEL; vt++) {
    if (ana_.isUnseenAndActive(*vt) && ana_.exploreRemainingCompOf(*vt)) {
      Component *p_new_comp = ana_.makeComponentFromArcheType();
      auto *packed_comp = new CacheableComponent(ana_.getArchetype().current_comp_for_caching_);
      CacheEntryID cache_entry_id;
      if (!cache_.manageNewComponent(top, *packed_comp, p_new_comp, cached_assn, cache_entry_id)) {
        component_stack_.push_back(p_new_comp);
        p_new_comp->set_id(cache_.storeAsEntry(*packed_comp, super_comp.id()));
      } else {
        if (config_.perform_random_sampling_) {
          // Store any cached vars in case of a split.
          std::vector<ComponentVarAndCls> cached_comp_data = p_new_comp->getData();
          VariableIndex comp_i = 0;
          // Store the variables in the component for proper processing of emancipated variables
          while (cached_comp_data[comp_i] != varsSENTINEL) {
            cached_vars_.push_back(cached_comp_data[comp_i]);
            comp_i++;
          }
          cached_comp_ids_.push_back(cache_entry_id);
        }
        // Delete component as it was found in the cache
        delete packed_comp;
        delete p_new_comp;
      }
    }
  }

  top.set_unprocessed_components_end(component_stack_.size());
  sortComponentStackRange(new_comps_start_ofs, component_stack_.size());
  if (config_.perform_random_sampling_) {
    long split_width = component_stack_.size() - new_comps_start_ofs;
    if (split_width > 1 || (split_width == 1 && !cached_assn.empty())) {
      top.setAsComponentSplit();
      if (is_valid_top_tree_store_point)
        TopTreeSampler::StartComponentSplit(literal_stack, top.emancipated_vars());
      if (config_.verbose) {
        std::stringstream ss;
        ss << "Component #" << (new_comps_start_ofs - 1) << " split into " << (split_width)
           << " pieces at depth " << StackLevel::componentSplitDepth() << ".";
        PrintInColor(ss, COLOR_CYAN);
      }
    }
    if (!cached_assn.empty())
      top.set_cached_assn(cached_assn);
  }
  // For simplicity, sort the cached variable list
  if (cached_vars_.size() > 1) {
    sort(cached_vars_.begin(), cached_vars_.end());
    sort(cached_comp_ids_.begin(), cached_comp_ids_.end());
  }
}


void ComponentManager::buildResidualComponent(StackLevel &top) {
  Component & super_comp = superComponentOf(top);

  std::vector<VariableIndex> comp_vars;
  std::vector<ClauseIndex> comp_cls;

  // Initialize data structures for the component analyzer
  ana_.setupAnalysisContext(top, super_comp);

  // Go through all components to be spawned and build a list of vars and comp
  for (auto vt = super_comp.varsBegin(); *vt != varsSENTINEL; vt++) {
    if (ana_.isUnseenAndActive(*vt) && ana_.exploreRemainingCompOf(*vt)) {
      Component *p_new_comp = ana_.makeComponentFromArcheType();
      auto comp_data = p_new_comp->getData();

      for (VariableIndex i = 0; i < p_new_comp->num_variables(); i++)
        comp_vars.push_back(comp_data[i]);
      for (ClauseIndex i = 0; i < p_new_comp->numLongClauses(); i++) {
        ClauseIndex cls = comp_data[p_new_comp->clauses_ofs() + i];
        comp_cls.push_back(cls);
      }
      delete p_new_comp;
    }
  }
  delete component_stack_.back();
  component_stack_.pop_back();

  // Sort the lists as components are build from only sorted lists
  std::sort(comp_vars.begin(), comp_vars.end());
  std::sort(comp_cls.begin(), comp_cls.end());

  auto * initial_component = new Component();
  for (auto var : comp_vars)
    initial_component->addVar(var);
  initial_component->closeVariableData();
  for (auto cls : comp_cls)
    initial_component->addCl(cls);
  initial_component->closeClauseData();
  component_stack_.push_back(initial_component);
}
/**
 * Helper function for determining freed variables.
 *
 * It goes through a list of variables and marks any that appears in the
 * descendant list.  If it is appears, then it is by definition not free.
 *
 * @param ref_vars Superset list of reference variables.
 * @param ref_num_vars Number of variables in @see ref_vars (exclusive)
 * @param varInDescendants One to one mapping of variables currently or previous
 * marked as non-free.  A variable should only ever be marked as non-free (i.e., true).
 * @param descendant_vars Subset of @see ref_vars.  This is a set of non-free vars.
 * @param desc_num_vars Number of variables in @see descendant_vars (exclusive)
 *
 * @return Number of descendant elements marked as present
 */
VariableIndex UpdateVarDescedantsList(const std::vector<ComponentVarAndCls> &ref_vars,
                                      VariableIndex ref_num_vars,
                                      std::vector<bool> &varInDescendants,
                                      std::vector<ComponentVarAndCls> descendant_vars,
                                      VariableIndex desc_num_vars);
//bool varInVector(long var_num, std::vector<unsigned int> vec);

#endif /* COMPONENT_MANAGEMENT_H_ */
