/**
 * component_cache.h
 *
 * Purpose: Defines the ComponentCache() class for storing components in the cache.
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

#ifndef COMPONENT_CACHE_H_
#define COMPONENT_CACHE_H_

#include <gmpxx.h>

#include <vector>
#include <sstream>

#include "component_types/cacheable_component.h"
#include "statistics.h"
#include "component_types/component.h"
#include "stack.h"
#include "solver_config.h"
#include "cached_assignment.h"

/**
 * Data structure that is used as the component cache for
 * the number of solutions for different components.
 */
class ComponentCache {
 public:
  ComponentCache(DataAndStatistics &statistics, SolverConfiguration &config);

  ~ComponentCache() {
    // debug_dump_data();
    for (auto &pentry : entry_base_) {
//      if (pentry != nullptr)
//        delete pentry;
      delete pentry;  // Do NOT need the check as deleting nullptr has no effect
    }
  }

  void init(Component &super_comp, bool quiet = false);

  // compute the size in bytes of the component cache from scratch
  // the value is stored in bytes_memory_usage_
  /**
   * Rebuild the size of the cache infrastructure (in bytes) from scratch by
   * computing the size of the individual component including:
   * <ul>
   *    <li>Size of the cache entry ID table (based on its capacity</li>
   *    <li>Size of the cacheable component area.</li>
   *    <li>Size of the cache entry slots used to store the list of free entries.</li>
   * </ul>
   *
   * It also updates the statistics for the formula including the
   * cache infrastructure bytes usage,
   *
   * This is called when the cache entries are deleted.
   *
   * @return Size of the cache infrastructure.
   */
  uint64_t compute_byte_size_infrastructure();
  /**
   * Access a reference to an entry in the component cache.  If the specified
   * entry ID does not exist, the program will crash. This is the constant version
   * of the function and is safe to use in assert statements.
   *
   * @param id Cache entry identification number
   *
   * @return Reference to the specified cache entry ID.
   */
  const CacheableComponent &entry_const(const CacheEntryID id) const {
    assert(entry_base_.size() > id);  // Verify the ID is value for the size.
    assert(entry_base_[id] != nullptr);
    return *entry_base_[id];
  }
  /**
   * Access a reference to an entry in the component cache.  If the specified
   * entry ID does not exist, the program will crash.
   *
   * @param id Cache entry identification number
   * @return Reference to the specified cache entry ID.
   */
  CacheableComponent &entry(const CacheEntryID id) {
    return const_cast<CacheableComponent&>(entry_const(id));
  }
//  /**
//   * Uses a cached component's CacheEntryID to extract the components
//   * cached information.
//   *
//   * @param comp Component to be extracted from the cache.
//   * @return Cache entry for the specified component.
//   */
//  CacheableComponent &entry(const Component& comp) {
//    return entry(comp.id());
//  }
  /**
   * Determines whether the
   *
   * @param id removeFromHashTable(CacheEntryID id);
   * @return
   */
  bool hasEntry(CacheEntryID id) const {
    assert(entry_base_.size() > id);
    return static_cast<bool>(entry_base_[id]);
  }
  /**
   * Removes the entry ID from the hash table but not from the entry base.
   *
   * @param id Identification number of the cache entry.
   */
  inline void removeFromHashTable(CacheEntryID id);

  // we delete the Component with ID id
  // and all its descendants from the cache
  inline void cleanPollutionsInvolving(CacheEntryID id);

  // creates a CCacheEntry in the entry base
  // which contains a packed copy of comp
  // returns the id of the entry created
  // stores in the entry the position of
  // comp which is a part of the component stack
  inline CacheEntryID storeAsEntry(CacheableComponent &ccomp,
                                   CacheEntryID super_comp_id);
  /**
   * New Component Cache Checker
   *
   * MT - Check quickly if the model count of the component is cached.
   * If so, incorporate it into the model count of top
   * If not, store the packed version of it in the entry_base of the cache
   *
   * @param top Top of the literal decision stack.
   * @param packed_comp New component to check whether it is in cache.
   * @param cached_assn Cached assignment to be modified
   * @return True if the new component is in cache.
   */
  bool manageNewComponent(StackLevel &top, CacheableComponent &packed_comp,
                          Component * unpacked_comp, CachedAssignment &cached_assn,
                          CacheEntryID &cache_entry_id) {
    statistics_.num_cache_look_ups_++;
    unsigned table_ofs =  packed_comp.hashkey() & table_size_mask_;

    CacheEntryID act_id = table_[table_ofs];
    while (act_id) {
      if (entry(act_id).equals(packed_comp)) {
        if (config_.verbose) {
          std::stringstream ss;
          ss << "\tCache hit for component #" << top.super_component() << " with model count "
             << entry(act_id).model_count() << " for cache ID #" << act_id;
          PrintInColor(ss, COLOR_YELLOW);
        }
//        if (config_.perform_random_sampling_ && config_.perform_sample_caching()) {
        if (config_.perform_sample_caching()) {
          // ToDo Add support for multiple samples in cache.
          assert(config_.num_samples_to_cache_ == 1);
          VariableIndex num_vars = unpacked_comp->num_variables();
          VariableIndex bit_shift = CACHED_VARIABLE_LEN * num_vars;
          mpz_class actual_model_count = entry(act_id).model_count() >> bit_shift;
          top.includeSolution(actual_model_count);
          cached_assn.ProcessComponent(unpacked_comp, entry(act_id).model_count());
        } else {
          top.includeSolution(entry(act_id).model_count());
        }
        cache_entry_id = act_id;  // Store the cache ID
        statistics_.incorporate_cache_hit(packed_comp);
        return true;
      }
      act_id = entry(act_id).next_bucket_element();
    }
    return false;
  }


  // unchecked erase of an entry from entry_base_
  void eraseEntry(CacheEntryID id) {
    statistics_.incorporate_cache_erase(*entry_base_[id]);
    delete entry_base_[id];
    entry_base_[id] = nullptr;
    free_entry_base_slots_.push_back(id);
  }

  // store the number in model_count as the model count of CacheEntryID id
  /**
   * Stores the model count of the specified CacheEntryId item in the
   * cache.
   *
   * If may invoke a cache resize if the new entry cannot fit in the cache.
   *
   * @param id  ID number of the cache entry.
   * @param model_count Number of models for the specified ID.
   */
  inline void storeValueOf(CacheEntryID id, const mpz_class &model_count);
  /**
   * Purges old entries from the cache to get the cache size to roughly half of its
   * current size.
   *
   * It then rehashes the table and updates the cache size information.
   *
   * @return Always true.
   */
  bool deleteEntries();

  // delete entries, keeping the descendants tree consistent
  inline void removeFromDescendantsTree(CacheEntryID id);

  // test function to ensure consistency of the descendant tree
  inline void test_descendantstree_consistency();

  void debug_dump_data();

 private:
  /**
   * Determine if the cache size merits it being rehashed.
   *
   * If it does, rehash the table.
   */
  void considerCacheResize() {
    if (entry_base_.size() > table_.size()) {
      reHashTable(2*table_.size());
    }
  }

  void reHashTable(unsigned long size) {
    table_.clear();
    table_.resize(size, 0);
    // we assert that table size is a power of 2
    // otherwise the table_size_mask_ doesn't work
    assert((table_.size() & (table_.size() - 1)) == 0);
    table_size_mask_ = table_.size() - 1;
    std::cout << "ZH - Rehash the table.\n";
    std::cout << "ts " << table_.size() << " " << table_size_mask_ << std::endl;
    unsigned collisions = 0;
    for (unsigned id = 2; id < entry_base_.size(); id++) {
      if (entry_base_[id] != nullptr) {
        entry_base_[id]->set_next_bucket_element(0);
        if (entry_base_[id]->modelCountFound()) {
          unsigned table_ofs = tableEntry(id);
          collisions += (table_[table_ofs] > 0 ? 1 : 0);
          entry_base_[id]->set_next_bucket_element(table_[table_ofs]);
          table_[table_ofs] = id;
        }
      }
    }
    std::cout << "ZH - Number of collisions.\n" << "coll " << collisions << std::endl;
  }

  unsigned tableEntry(CacheEntryID id) {
    return entry(id).hashkey() & table_size_mask_;
  }

  void add_descendant(CacheEntryID compid, CacheEntryID descendantid) {
    assert(descendantid != entry_const(compid).first_descendant());
    entry(descendantid).set_next_sibling(entry_const(compid).first_descendant());
    entry(compid).set_first_descendant(descendantid);
  }

  void remove_firstdescendantOf(CacheEntryID compid) {
    CacheEntryID desc = entry(compid).first_descendant();
    if (desc != 0)
      entry(compid).set_first_descendant(entry(desc).next_sibling());
  }

  std::vector<CacheableComponent *> entry_base_;
  std::vector<CacheEntryID> free_entry_base_slots_;

  // the actual hash table
  // by means of which the cache is accessed
  /**
   * Actual hash table storing the cache entries.
   * It is the means by which the cache is accessed.
   *
   * Vector of unsigned integers.
   */
  std::vector<CacheEntryID> table_;

  /**
   * Table size is a power of 2.  If its size is 2^n, then
   * the mask is 2^n - 1.
   */
  unsigned table_size_mask_ = INT_MAX;

  DataAndStatistics &statistics_;

  unsigned long my_time_ = 0;

  /**
   * Reference to the outer solver configuration.  This was added for use
   * by the sampler.
   */
  SolverConfiguration &config_;
};


#include "component_cache.inc"


#endif /* COMPONENT_CACHE_H_ */
