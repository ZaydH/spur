/**
 * top_tree_sampler.h
 *
 * Purpose: Defines the TopTreeSampler() class that can be used to select partial assignments
 * for top of the compact counting tree.
 *
 * @author Zayd Hammoudeh <zayd@ucsc.edu>
 * @version 0.00.00
 *
 * Copyright (C) 2018 Zayd Hammoudeh.
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms of the MIT license.  See the
 * LICENSE file for details.
 */

#ifndef TOP_TREE_SAMPLER_H
#define TOP_TREE_SAMPLER_H

#include <gmp.h>

#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <iomanip>

#include "primitive_types.h"
#include "solver_config.h"
#include "statistics.h"
#include "sampler_tools.h"
#include "rand_distributions.h"
#include "structures.h"

#define TOP_TREE_DATA_TYPE_SEP " "
#define TOP_TREE_FILE_ITEM_SEP ","

#define TOP_TREE_COMPONENT_SPLIT_START  "{"
#define TOP_TREE_COMPONENT_SPLIT_END    "}"
#define TOP_TREE_COMPONENT_START  "["
#define TOP_TREE_COMPONENT_END    "]"
#define TOP_TREE_CYLINDER_KEY "CYLINDER"
#define TOP_TREE_MAX_DEPTH_KEY "MAX_DEPTH"
// When only a trivial formula left at max depth, can have any empty literal string
#define TOP_TREE_EMPTY_LIT_STACK "NONE"


class TopTreeSampler {
 private:
  enum class TopTreeSampleNodeType {
    ROOT,
    MAX_DEPTH_NODE,
    CYLINDER_NODE,
    COMPONENT_SPLIT,
    COMPONENT,
    EMPTY_CONSTRUCTOR
  };

  struct TopTreeNode {
    TopTreeSampleNodeType node_type_;
    std::vector<TopTreeLiteral> assigned_lits_;
    std::vector<VariableIndex> emancipated_vars_;
    mpz_class model_count;
    /**
     * Number of samples associated with this node and all of its descendants
     */
    SampleSize num_samples = 0;
    /**
     * Fully explored nodes automatically get passed to the partial_outputs set.
     */
    bool fully_explored_ = true;
    mpz_class component_split_multiplier;
    std::vector<TopTreeNode> descendants;

    explicit TopTreeNode(TopTreeSampleNodeType node_type = TopTreeSampleNodeType::EMPTY_CONSTRUCTOR,
                         const std::vector<TopTreeLiteral> &assn_lits
                                                                    = std::vector<TopTreeLiteral>(),
                         const std::vector<VariableIndex> &emancipated_vars
                                                                     = std::vector<VariableIndex>(),
                         const mpz_class &model_count = 0)
        : node_type_(node_type), assigned_lits_(assn_lits), emancipated_vars_(emancipated_vars) {
      if (node_type == TopTreeSampleNodeType::MAX_DEPTH_NODE)
        fully_explored_ = false;
      descendants.reserve(5);
      this->model_count = model_count;
    }
    /**
     * Adds the assigned literals and emancipated variables of one top tree node to the implicit
     * parameter.
     *
     * @param other The TopTreeNode() whose emancipated vars and assigned literals will be added
     */
    void addAssignedLitsAndEmancipatedVars(const TopTreeNode &other) {
      if (!other.assigned_lits_.empty())
        assigned_lits_.insert(assigned_lits_.end(),
                              other.assigned_lits_.begin(), other.assigned_lits_.end());
      assert(verifyNoDuplicatesAssignedLits());

      if (!other.emancipated_vars_.empty())
        emancipated_vars_.insert(emancipated_vars_.end(),
                                 other.emancipated_vars_.begin(), other.emancipated_vars_.end());
    }
    /**
     * Merges a set of |s| component into the larger component split.  Before performing the merge,
     * it randomly shuffles the components.
     *
     * @param other A set of component objects.
     */
    void randomComponentMerge(TopTreeNode &other) {
      assert(descendants.size() == num_samples && num_samples == other.num_samples
             && other.node_type_ == TopTreeSampleNodeType::COMPONENT);

      // Build list of IDs to transfer
      std::vector<SampleSize> other_desc_ids;
      other_desc_ids.reserve(num_samples);
      for (SampleSize i = 0; i < other.descendants.size(); i++)
        for (SampleSize j = 0; j < other.descendants[i].num_samples; j++)
          other_desc_ids.emplace_back(i);
      Random::shuffle(other_desc_ids);

      for (SampleSize i = 0; i < descendants.size(); i++) {
        // Copy the assigned literals from the components
        std::vector<TopTreeLiteral> &this_desc = descendants[i].assigned_lits_;

        SampleSize desc_id = other_desc_ids[i];
        std::vector<TopTreeLiteral> &other_desc = other.descendants[desc_id].assigned_lits_;

        this_desc.insert(this_desc.end(), other_desc.begin(), other_desc.end());
        assert(descendants[i].verifyNoDuplicatesAssignedLits());
      }
    }
    /**
     * Compare two SampleSetTreeNode() objects and order them by the one having the most samples
     * associated with it.
     *
     * @param first A SampleSetTreeNode() object.
     * @param second A SampleSetTreeNode() object.
     * @return True if \p first has more sample than \p second.
     */
    static bool compare_num_samples(const TopTreeNode &first, const TopTreeNode &second) {
      return first.num_samples > second.num_samples;
    }
    /**
     * Sort the partial assignments such that the one with the most samples requested is first.
     */
    void SortBySampleCount() {
      std::sort(descendants.begin(), descendants.end(), TopTreeNode::compare_num_samples);
    }
    /**
     * Copies the descendants of one SampleSetTreeNode to another.  It only copies those
     * descendants that have at least one sample associated with them.
     *
     * @param other SampleSetTree() node object whose descendants will be transferred to the
     * implicit parameter.
     */
    void addDescendants(const TopTreeNode &other) {
      assert(!other.descendants.empty());
      for (const auto &desc : other.descendants) {
        assert(desc.num_samples > 0);
        descendants.emplace_back(desc);
        num_samples += desc.num_samples;
      }
//      descendants.insert(descendants.end(), other.descendants.begin(), other.descendants.end());
    }
    /**
     * Creates a new set of descendants
     *
     * @param num_new_descendants Number of new descendants to create.
     */
    void createEmptyDescendants(SampleSize num_new_descendants) {
      SampleSize previous_size = descendants.size();
      descendants.resize(descendants.size() + num_new_descendants);
      for (SampleSize i = previous_size; i < descendants.size(); i++)
        descendants[i].num_samples = 1;
    }
    /**
     * Checks whether the node has the same partial assignment as the specified node.
     *
     * @param other Another SampleSetTreeNode.
     * @return True if the partial assignment matches.
     */
    const bool HasSamePartialAssignment(const TopTreeNode &other) {
      if (assigned_lits_.size() != other.assigned_lits_.size())
        return false;

      for (VariableIndex i = 0; i < assigned_lits_.size(); i++)
        if (assigned_lits_[i] != other.assigned_lits_[i])
          return false;
      return true;
    }
    /**
     * Used to create an ordering of SampleSetTreeNode() objects.  It orders by the length of
     * the assigned literals first.  If they have equal lengths, then it sorts by the assigned
     * literal stack contents themselves (starting from index 0).
     *
     * @param first A SampleSetTreeNode() object.
     * @param second Another SampleSetTreeNode() object.
     * @return True if \p first should be before \p second.
     */
    static const bool compare_partial_assn(const TopTreeNode &first, const TopTreeNode &second) {
      if (first.assigned_lits_.size() != second.assigned_lits_.size())
        return first.assigned_lits_.size() < second.assigned_lits_.size();
      for (VariableIndex i = 0; i < first.assigned_lits_.size(); i++)
        if (first.assigned_lits_[i] != second.assigned_lits_[i])
          return first.assigned_lits_[i] < second.assigned_lits_[i];
      return false;
    }

   private:
    /**
     * Verify that there are no duplicate literals in the assigned literal set.
     *
     * @return false if the node has duplicate assigned literals.
     */
    const bool verifyNoDuplicatesAssignedLits() const {
      std::vector<VariableIndex> temp_assn_vars;
      temp_assn_vars.reserve(assigned_lits_.size());
      for (auto &lit : assigned_lits_)
        temp_assn_vars.emplace_back(std::abs(lit));

      std::sort(temp_assn_vars.begin(), temp_assn_vars.end());
      for (VariableIndex i = 1; i < temp_assn_vars.size(); i++)
        if (temp_assn_vars[i] == temp_assn_vars[i-1])
          return false;
      return true;
    }
  };

  class BaseAssignments {
    friend class TopTreeSampler;
   public:
    /**
     * Gets the next unfinished formula that exists (if any).
     *
     * @param num_samples Number of samples associated with the partial formula.
     * @param part_assn
     * @return True if a formula was gotten.
     */
    const bool GetNext(SampleSize &num_samples, PartialAssignment &part_assn,
                       std::vector<VariableIndex> &emancipated_vars) {
      if (!is_input_open()) {
        assert(!input_filename_.empty());
        fin_.open(input_filename_);
        has_added_any_ = false;
      }

      std::string assigned_lits_str, emancipated_vars_str;
      if (!(fin_ >> num_samples) || !(fin_ >> assigned_lits_str)
          || !(fin_ >> emancipated_vars_str)) {
        fin_.close();
        return false;
      }

      // Build the partial assignment
      std::vector<TopTreeLiteral> assigned_lits;
      TokenizeLiteralStack(assigned_lits_str, assigned_lits);
      part_assn.clear();

      emancipated_vars.clear();
      TokenizeEmancipatedVars(emancipated_vars_str, emancipated_vars);
      BuildPartialAssignment(assigned_lits, part_assn);
      return true;
    }
    /**
     * Adds a new unfinished formula to the output file.
     *
     * @param sample Top tree sample with partial assignment and literals.
     */
    void addNew(const TopTreeNode &sample) {
      if (!is_output_open())
        fout_.open(output_filename_);
      else
        fout_ << "\n";
      WriteOutputAssignment(fout_, sample);

      has_added_any_ = true;
    }
    /**
     * Checks whether any formulas were added since the last time this object was opened as an
     * input.
     *
     * @return True if this object has any assignments/formulas associated with it.
     */
    const bool has_added_any() const { return has_added_any_; }
    /**
     * Clean up the input and output files.
     */
    virtual void cleanup() {
      assert(FileExists(input_filename_));
      std::remove(input_filename_.c_str());

      if (input_filename_ != output_filename_) {
        assert(FileExists(output_filename_));
        std::remove(output_filename_.c_str());
      }
    }

   private:
    /**
     * Sets the input filename for these formulas.
     *
     * @param input_filename New name for reach the input file names.
     */
    void set_input_filename(const std::string &input_filename) {
      input_filename_ = input_filename;
    }
    /**
     * Sets the filename that formulas are written to.
     *
     * @param output_filename New output filename.
     */
    void set_output_filename(const std::string &output_filename) {
      output_filename_ = output_filename;

      // Delete any previous version of this file.
      if (FileExists(output_filename_))
        std::remove(output_filename_.c_str());
    }
    /**
     * Close the formula output file.
     */
    void close_output() {
      assert(fout_.is_open());
      fout_ << std::flush;
      fout_.close();
    }
    /**
     * Checks whether the file to write output of the formulas is already open.
     *
     * @return True if the output file is already open.
     */
    const bool is_output_open() const { return fout_.is_open(); }
    /**
     * Checks whether the file to where outputted formulas are read back is already open.
     *
     * @return True if the input file is already open.
     */
    const bool is_input_open() const { return fin_.is_open(); }
    /**
     * Stores whether any incomplete formulas still exist.
     */
    std::string output_filename_ = "";
    /**
     * Filename for the input set of formulas.
     */
    std::string input_filename_ = "";
    /**
     * Open the output formula file stream.
     */
    std::ofstream fout_;
    /**
     * Open the input formula file stream.
     */
    std::ifstream fin_;
    /**
     * Stores whether this has any added any formulas since the last time the output was opened.
     */
    bool has_added_any_ = false;
  };

  class RemainingFormulas : public BaseAssignments {
   public:
    /**
     * Initialize data structures associated with the formulas that still require additional top
     * tree.
     */
    RemainingFormulas() {
      set_output_filename(buildOutputFilename(num_rounds_));
    }
    /**
     * Updates the unfinished formula file for another round of analysis.
     */
    void RoundCompleted() {
      if (!config_->quiet) {
        std::stringstream ss;
        ss << "Round #" << num_rounds_ << " completed.";
        PrintInColor(ss, COLOR_CYAN);
      }
      num_rounds_++;
      if (!has_added_any())
        return;

      assert(is_output_open() && !is_input_open());

      fout_ << std::flush;
      fout_.close();
      set_input_filename(output_filename_);

      // New inputs are read from the last output.
      set_output_filename(buildOutputFilename(num_rounds_));

      // Restore the configuration and statistics to the root solver's
      config_ = root_config_;
      statistics_ = root_statistics_;
    }
    /**
     * Delete any files created as part of the remaining formula storage.
     */
    void cleanup() override {
      for (auto file_cnt = FIRST_FILE_NUMBER; file_cnt < num_rounds_; file_cnt++)
        std::remove((buildOutputFilename(file_cnt)).c_str());
    }

   private:
    /**
     * Stores the number for the first remaining formulas file created.
     */
    const uint32_t FIRST_FILE_NUMBER = 1;
    /**
     * Number of rounds of incomplete top tree performed.
     */
    uint32_t num_rounds_ = FIRST_FILE_NUMBER;
    /**
     * Maximum zero padding width when building the output filename.
     */
    const int WIDTH_OF_FILE_NUMBER = 5;
    /**
     * Builds the output filename for writing the remaining formulas.
     *
     * @param file_number File number to be given to the output file.  It must be greater than 0.
     * @return Filename.
     */
    inline std::string buildOutputFilename(uint32_t file_number) {
      assert(file_number >= FIRST_FILE_NUMBER);
      std::stringstream ss;
      ss << "./__incomplete_top_tree_"
         << std::setfill('0') << std::setw(WIDTH_OF_FILE_NUMBER) << num_rounds_ << ".txt";
      return ss.str();
    }
  };

 public:
  /**
   * Unfinished formulas object used for managing and controlling write files.
   */
  static RemainingFormulas remaining_formulas;
  /**
   * Set of partial assignments output from the top tree solver.
   */
  static BaseAssignments output_partials;
  /**
   * Initializes all data structures needed by the top tree sampler.
   *
   * @param config Master solver configuration
   *
   * @param statistics Master solver statistics object.
   */
  static void initialize(SolverConfiguration &config, DataAndStatistics &statistics) {
    all_top_tree_filename_ = "." FILE_PATH_SEPARATOR "__all_top_tree_samples.txt";
    if (FileExists(all_top_tree_filename_))
      std::remove(all_top_tree_filename_.c_str());

    total_models_ = 0;
    tot_num_rough_samples_ = 0;

    root_config_ = &config;
    root_statistics_ = &statistics;

    // Initialize the unfinished formulas data structures
//    remaining_formulas = RemainingFormulas();

    // Configure where the complete assignments will go.
//    output_partials = BaseAssignments();
    output_partials.set_input_filename(root_config_->top_tree_samples_output_file_);
    output_partials.set_output_filename(root_config_->top_tree_samples_output_file_);
  }
  /**
   * Called before beginning top tree sampling.  It resets any residual data structures that
   * may have been set in a previous run.
   */
  static void start(SolverConfiguration &config, DataAndStatistics &statistics) {
    config_ = &config;
    assert(config_->perform_top_tree_sampling);
    statistics_ = &statistics;

    statistics_->num_top_tree_nodes_ = 0;

    total_models_ = 0;
    tot_num_rough_samples_ = 0;

    if (FileExists(all_top_tree_filename_))
      std::remove(all_top_tree_filename_.c_str());

    assert(!all_top_tree_fout_.is_open());
  }
  /**
   *
   * Each output sample is on its own line.  The basic format is as follows:
   *
   * <SampleModelCount> <LiteralStack1,LiteralStack2,...>
   *
   * @param literal_stack Variables in the literals stack.
   * @param Set of emancipated variables.
   * @param model_count Number of models for this node in the tree.
   */
  //@{
  /**
   * @brief Uses the count multiplier with the unscaled model count.
   *
   * @param unscaled_model_count Unscaled model count that does not account for variables
   * emancipated by earlier variable setting.
   *
   * @param count_multiplier Multiplier for the unscaled model count.
   */
  static void StoreSample(const mpz_class &unscaled_model_count, const mpz_class &count_multiplier,
                          const std::vector<LiteralID> &literal_stack,
                          const std::vector<VariableIndex> &emancipated_vars,
                          TopTreeNodeType node_type) {
    mpz_class total_model_count = unscaled_model_count;
    if (count_multiplier != 1)
      total_model_count *= count_multiplier;
    StoreSample(total_model_count, literal_stack, emancipated_vars, node_type);
  }
  /**
   * @brief Function with the node's total model count supplied.
   */
  static void StoreSample(const mpz_class &model_count,
                          const std::vector<LiteralID> &literal_stack,
                          const std::vector<VariableIndex> &emancipated_vars,
                          TopTreeNodeType node_type) {
    assert(config_ && statistics_);
    assert(model_count > 0 && !literal_stack.empty());

    // Update the stats
    TopTreeNodeType stats_node_type;
    if (component_split_depth() > 0) {
      stats_node_type = TopTreeNodeType::COMPONENT_SPLIT;
    } else {
      if (node_type == TopTreeNodeType::MAX_DEPTH)
        stats_node_type = TopTreeNodeType::MAX_DEPTH;
      else
        stats_node_type = TopTreeNodeType::CYLINDER;
    }
    statistics_->UpdateNodeTypeStatistics(stats_node_type, model_count);

    StartTopTreeFoutLine();
    if (node_type == TopTreeNodeType::MAX_DEPTH) {
      all_top_tree_fout_ << TOP_TREE_MAX_DEPTH_KEY;
    } else if (node_type == TopTreeNodeType ::CYLINDER) {
      all_top_tree_fout_ << TOP_TREE_CYLINDER_KEY;
    }
    all_top_tree_fout_ << TOP_TREE_DATA_TYPE_SEP << model_count << TOP_TREE_DATA_TYPE_SEP;

    WriteLiteralStackOut(all_top_tree_fout_, literal_stack);
    WriteEmancipatedVars(all_top_tree_fout_, emancipated_vars);
    all_top_tree_fout_ << std::flush;

    // Variables kept for simplified bookkeeping
    if (component_split_depth() == 0)
      total_models_ += model_count;
    statistics_->num_top_tree_nodes_++;
  }
  //@}
  /**
   * Opens a new component split in the results file.  The format of the output is:
   *
   * } <ComponentSplitDepth> <LiteralStack> <EmancipatedVariables>
   *
   * @param depth Component split depth.
   */
  static void StartComponentSplit(const std::vector<LiteralID> &literal_stack,
                                  const std::vector<VariableIndex> &emancipated_vars) {
    StartTopTreeFoutLine();

    all_top_tree_fout_ << TOP_TREE_COMPONENT_SPLIT_START << TOP_TREE_DATA_TYPE_SEP
                       << (component_split_depth() + 1)  // Add 1 since not updated depth yet
                       << TOP_TREE_DATA_TYPE_SEP;
    WriteLiteralStackOut(all_top_tree_fout_, literal_stack);
    WriteEmancipatedVars(all_top_tree_fout_, emancipated_vars);
    all_top_tree_fout_ << std::flush;

    VariableIndex ofs = (!literal_stack_start_loc_.empty()) ? literal_stack_start_loc_.back() : 0;
    debug_lit_stack_.insert(std::end(debug_lit_stack_), std::begin(literal_stack) + ofs,
                            std::end(literal_stack));
    assert(debug_lit_stack_ == literal_stack);
    // Update component split depth.
    literal_stack_start_loc_.emplace_back(literal_stack.size());
  }
  /**
   * Mark the end of a component split in the top tree. It scales the specified model
   * count by the \p actual_multiplier.
   *
   * @param depth Component  split depth.
   * @param model_count Model count of the component split unscaled for any preceding
   * emancipated variables.
   * @param report_multiplier Multiplier reported in the file for the count but not what is
   * multiplied.
   * @param count_multiplier Multiplier for the component split actually used to calculated the
   * total models of this component split.
   */
  static void CloseComponentSplit(VariableIndex depth, const mpz_class &model_count,
                                  const mpz_class &report_multiplier,
                                  const mpz_class &actual_multiplier = 1) {
    assert(component_split_depth() == depth);
    mpz_class total_models = model_count;
    if (actual_multiplier != 1)
      total_models *= actual_multiplier;
    literal_stack_start_loc_.pop_back();
    if (component_split_depth() == 0)
      total_models_ += total_models;

    StartTopTreeFoutLine();
    all_top_tree_fout_ << TOP_TREE_COMPONENT_SPLIT_END << TOP_TREE_DATA_TYPE_SEP
                       << depth << TOP_TREE_DATA_TYPE_SEP
                       << report_multiplier << TOP_TREE_DATA_TYPE_SEP
                       << total_models << std::flush;

    // DEBUG ONLY
    if (!literal_stack_start_loc_.empty())
      debug_lit_stack_.resize(literal_stack_start_loc_.back());
    else
      debug_lit_stack_.clear();
  }
  /**
   * Mark the beginning of a component.
   */
  static void StartComponent() {
    StartTopTreeFoutLine();
    all_top_tree_fout_ << TOP_TREE_COMPONENT_START << TOP_TREE_DATA_TYPE_SEP << std::flush;
  }
  /**
   * For tracking purposes. It marks the end of a component split.
   *
   * @param unscaled_model_count Model count for the component without any preceding variable
   * multipliers.
   */
  static void CloseComponent(const mpz_class &unscaled_model_count) {
    StartTopTreeFoutLine();
    all_top_tree_fout_ << TOP_TREE_COMPONENT_END << TOP_TREE_DATA_TYPE_SEP
                       << unscaled_model_count << std::flush;
  }
  /**
   * Finalizes the top tree sampler.  It processes the complete list of all tree nodes
   * and then selects the specified number of samples.  These are then written to a file.
   */
  static void FinalizeTopTreeSamples() {
    all_top_tree_fout_ << std::flush;
    all_top_tree_fout_.close();

    if (!config_->quiet)
      printTreeSolverStatistics(std::cout);
    if (statistics_->final_solution_count_ != total_models_)
      ExitWithError("Invalid top tree solution count.");

    TopTreeNode tree_sample_set(TopTreeSampleNodeType::ROOT);
    std::ifstream all_top_tree_fin(all_top_tree_filename_);
    BuildRoughSampleSet(all_top_tree_fin, tree_sample_set, 0);
    all_top_tree_fin.close();
    assert(tree_sample_set.model_count == statistics_->final_solution_count_);

    tree_sample_set.num_samples = config_->num_samples_;
    DownsampleRoughSet(tree_sample_set);
    assert(TopTreeSampler::verifySampleSizes(tree_sample_set));

    TopTreeNode final_tree_set(TopTreeSampleNodeType::ROOT);
    BuildFinalTopTreeSet(tree_sample_set, final_tree_set, TOP_COMPONENT_SPLIT_DEPTH_);
    // For simplified analysis of elements needing re-analysis, sort the list of nodes.
    final_tree_set.SortBySampleCount();
    assert(verifySampleSizes(final_tree_set));

    WriteFinalTopTreeSampleSet(final_tree_set);
  }
  /**
   * Prints statistics about the number and type of top tree nodes found during execution.
   *
   * @param out Output stream to write to (e.g., console, file, etc.)
   */
  static void printTreeSolverStatistics(std::ostream &out) {
    bool delete_stats = false;
    if (!statistics_) {
      statistics_ = new DataAndStatistics();
      delete_stats = true;
    }
    out << "num_top_tree_nodes," << statistics_->num_top_tree_nodes_ << "\n"
        << "num_cylinder_nodes,"
        << statistics_->num_tree_nodes(TopTreeNodeType::CYLINDER) << "\n"
//        << "num_max_depth_nodes,"
//        << statistics_->num_tree_nodes(TopTreeNodeType::MAX_DEPTH) << "\n"
//        << "num_component_split_nodes,"
//        << statistics_->num_tree_nodes(TopTreeNodeType::COMPONENT_SPLIT) << "\n"
        << "num_cylinder_models,"
        << statistics_->num_tree_node_models(TopTreeNodeType::CYLINDER) << "\n"
//        << "num_max_depth_models,"
//        << statistics_->num_tree_node_models(TopTreeNodeType::MAX_DEPTH) << "\n"
//        << "num_component_split_models,"
//        << statistics_->num_tree_node_models(TopTreeNodeType::COMPONENT_SPLIT) << "\n"
//        << "percent_cylinder_models,"
//        << statistics_->percent_tree_node_models(TopTreeNodeType::CYLINDER) << "\n"
//        << "percent_max_depth_models,"
//        << statistics_->percent_tree_node_models(TopTreeNodeType::MAX_DEPTH) << "\n"
//        << "percent_component_split_models,"
//        << statistics_->percent_tree_node_models(TopTreeNodeType::COMPONENT_SPLIT) << "\n"
        << std::flush;
    if (delete_stats) {
      delete statistics_;
      statistics_ = nullptr;
    }
  }
  /**
   * Converts the literal stack string in a top tree sample file into a partial assignment
   * that can be passed to initialize a solver.
   *
   * @param literal_string Assigned literals representation string in the top tree sampler file.
   *
   * @param lits Vector of the literals as they would be in a DIMACS file.
   */
  inline static void TokenizeLiteralStack(const std::string &literal_string,
                                          std::vector<TopTreeLiteral> &lits) {
    assert(!literal_string.empty());
    char * lit_str;
    auto rest = const_cast<char*>(literal_string.c_str());
    lits.clear();
    // When only a trivial formula left at max depth, can have any empty literal string
    if (literal_string == TOP_TREE_EMPTY_LIT_STACK)
      return;
    while ((lit_str = strtok_r(rest, TOP_TREE_FILE_ITEM_SEP, &rest))) {
      long lit_val = strtol(lit_str, nullptr, STR_DECIMAL_BASE);
      lits.emplace_back(static_cast<TopTreeLiteral>(lit_val));
    }
  }
  /**
   * Emancipated variable comma-separated list into a
   *
   * @param literal_string Assigned literals representation string in the top tree sampler file.
   *
   * @param emancipated_vars Vector of the identification numbers of the variables in the
   * emancipated variable list.
   */
  inline static void TokenizeEmancipatedVars(const std::string &literal_string,
                                             std::vector<VariableIndex> &emancipated_vars) {
    assert(!literal_string.empty());
    char * var_str;
    auto rest = const_cast<char*>(literal_string.c_str());
    emancipated_vars.clear();
    // When only a trivial formula left at max depth, can have any empty literal string
    if (literal_string == TOP_TREE_EMPTY_LIT_STACK)
      return;

    while ((var_str = strtok_r(rest, TOP_TREE_FILE_ITEM_SEP, &rest))) {
      long var_id = strtol(var_str, nullptr, STR_DECIMAL_BASE);
      assert(var_id > 0);
      emancipated_vars.emplace_back(static_cast<VariableIndex>(var_id));
    }
  }
  /**
   * Delete any files created during top tree sampling.
   */
  inline static void cleanup() {
    output_partials.cleanup();
    remaining_formulas.cleanup();
    std::remove(all_top_tree_filename_.c_str());
  }

 private:
  /**
   * Disallows creating an instance.  All methods are static.
   */
  TopTreeSampler() = default;
  /**
   * Writes a subset of the literal stack to the specified stream.  The portion of the literal stack
   * written is that which is new after the last component split.  If no component split has
   * occurred, it writes the complete literal stack.
   *
   * The string of literal stack is formatted as a comma-separated list of the variables in
   * the literal stack.  The variables are negated if set to "false" and positive it set to "true."
   *
   * @param out Output stream to which to write the literal stack.
   * @param literal_stack Current literal stack.
   */
  inline static void WriteLiteralStackOut(std::ostream &out,
                                          const std::vector<LiteralID> &literal_stack) {
    VariableIndex start_loc;
    start_loc = (literal_stack_start_loc_.empty()) ? 0 : literal_stack_start_loc_.back();

    if (start_loc == literal_stack.size()) {
      all_top_tree_fout_ << TOP_TREE_EMPTY_LIT_STACK;
    } else {
      for (VariableIndex i = start_loc; i < literal_stack.size(); i++) {
        if (i != start_loc)
          all_top_tree_fout_ << TOP_TREE_FILE_ITEM_SEP;
        if (!literal_stack[i].sign())
          all_top_tree_fout_ << "-";
        all_top_tree_fout_ << literal_stack[i].var();
      }
    }
  }
  /**
   * Writes the complete emancipated variable list ot the function.  For a leaf (i.e., cylinder,
   * max depth node), the node's complete emancipated variable split will be broken up at the
   * component split points.
   *
   * @param out Output stream to which to write the literal stack.
   * @param emancipated_vars A set of emancipated variables.
   */
  inline static void WriteEmancipatedVars(std::ostream &out,
                                          const std::vector<VariableIndex> &emancipated_vars) {
    if (emancipated_vars.empty()) {
      all_top_tree_fout_ << TOP_TREE_EMPTY_LIT_STACK;
    } else {
      for (auto free_var : emancipated_vars) {
        if (free_var != emancipated_vars.front())
          all_top_tree_fout_ << TOP_TREE_FILE_ITEM_SEP;
        all_top_tree_fout_ << free_var;
      }
    }
  }
  /**
   * Constructs a partial assignment from the tokenized literal information in the top-tree
   * file.
   *
   * @param lits Literals in the top tree file.
   *
   * @return Partial assignment associated with the top of the tree samples.
   */
  inline static void BuildPartialAssignment(const std::vector<TopTreeLiteral> &lits,
                                            PartialAssignment &partial_assn) {
    if (partial_assn.empty())
      partial_assn.resize(statistics_->num_variables_ + FIRST_VAR, ASSN_U);
    assert(partial_assn.size() == statistics_->num_variables_ + FIRST_VAR);
    for (auto lit : lits) {
      assert(partial_assn[abs(lit)] == ASSN_U);
      partial_assn[abs(lit)] = (lit < 0) ? ASSN_F : ASSN_T;
    }
  }
  /**
   * If the file is open, add a new line to end the last entry otherwise open the file.
   */
  static void StartTopTreeFoutLine() {
    if (all_top_tree_fout_.is_open())
      all_top_tree_fout_ << "\n";
    else
      all_top_tree_fout_.open(all_top_tree_filename_);
    // Create number of tabs equal to component split depth.
    if (component_split_depth() > 0)
      all_top_tree_fout_ << std::string(component_split_depth(), '\t');
  }
  /**
   * Parses the complete top tree file and builds from it the top tree as a data structure.  This
   * function handles component splits recursively.
   *
   * Sample counts will eventually be added to the nodes of this tree.
   *
   * @param rough_sample_set Rough sample set to be created from the top tree samples.
   */
  static void BuildRoughSampleSet(std::ifstream &all_top_tree_fin,
                                  TopTreeNode &rough_sample_set,
                                  const VariableIndex &component_split_depth) {
    std::vector<TopTreeLiteral> sample_stack_lits;
    std::vector<VariableIndex> emancipated_vars;
    mpz_class node_model_count;
    std::string entry_header;
    mpz_class total_model_count = 0;
    while (all_top_tree_fin >> entry_header) {
      // At the end of the component, return and it will be processed by the component split.
      if (entry_header == TOP_TREE_COMPONENT_END) {
        assert(component_split_depth > 0);
        rough_sample_set.model_count = total_model_count;
        return;
      // Process a leaf node and store its partial formula and model count
      } else if (entry_header == TOP_TREE_CYLINDER_KEY || entry_header == TOP_TREE_MAX_DEPTH_KEY) {
        std::string lit_stack_list, emancipated_vars_list;
        all_top_tree_fin >> node_model_count >> lit_stack_list >> emancipated_vars_list;
        assert(node_model_count > 0);

        TopTreeSampleNodeType type;
        type = (entry_header == TOP_TREE_CYLINDER_KEY) ? TopTreeSampleNodeType::CYLINDER_NODE
                                                       : TopTreeSampleNodeType::MAX_DEPTH_NODE;
        TokenizeLiteralStack(lit_stack_list, sample_stack_lits);
        TokenizeEmancipatedVars(emancipated_vars_list, emancipated_vars);
        rough_sample_set.descendants.emplace_back(TopTreeNode(type, sample_stack_lits,
                                                              emancipated_vars, node_model_count));
      // Process a component split.
      } else if (entry_header == TOP_TREE_COMPONENT_SPLIT_START) {
        ProcessRoughSetCompSplit(all_top_tree_fin, rough_sample_set,
                                 component_split_depth + 1);
      }
      total_model_count += rough_sample_set.descendants.back().model_count;

      // A parent is only fully explored if its descendants are as well.
      rough_sample_set.fully_explored_ &= rough_sample_set.descendants.back().fully_explored_;
    }
    rough_sample_set.model_count = total_model_count;
  }
  /**
   * Recursively processes component splits in the all samples file.
   *
   * @param all_top_tree_fin All samples input_file stream.
   * @param rough_sample_set Rough set of samples created by parsing the file so far.
   * @param component_split_depth Current component split depth in the file parsing.
   */
  inline static void ProcessRoughSetCompSplit(std::ifstream &all_top_tree_fin,
                                              TopTreeNode &rough_sample_set,
                                              const VariableIndex &component_split_depth) {
    // Extract the remaining parameters for the component split from the file
    std::string lit_stack_list, emancipated_vars_list;
    VariableIndex new_comp_depth;
    all_top_tree_fin >> new_comp_depth >> lit_stack_list >> emancipated_vars_list;
    assert(new_comp_depth == component_split_depth);

    std::vector<TopTreeLiteral> sample_stack_lits;
    TokenizeLiteralStack(lit_stack_list, sample_stack_lits);
    std::vector<VariableIndex> emancipated_vars;
    TokenizeEmancipatedVars(emancipated_vars_list, emancipated_vars);
    rough_sample_set.descendants.emplace_back(TopTreeNode(TopTreeSampleNodeType::COMPONENT_SPLIT,
                                                          sample_stack_lits, emancipated_vars));

    TopTreeNode &comp_split = rough_sample_set.descendants.back();
    mpz_class component_split_model_count = 1;
    std::string entry_header;
    mpz_class node_model_count;
    while ((all_top_tree_fin >> entry_header) && entry_header == TOP_TREE_COMPONENT_START) {
      comp_split.descendants.emplace_back(TopTreeNode(TopTreeSampleNodeType::COMPONENT));
      BuildRoughSampleSet(all_top_tree_fin, comp_split.descendants.back(),
                          component_split_depth);
      all_top_tree_fin >> node_model_count;
      assert(comp_split.descendants.back().model_count == node_model_count);
      component_split_model_count *= node_model_count;

      // A component split is only fully explored if all of its descendants are as well
      comp_split.fully_explored_ &=  comp_split.descendants.back().fully_explored_;
    }

    mpz_class comp_split_multiplier;
    all_top_tree_fin >> new_comp_depth >> comp_split_multiplier >> node_model_count;
    assert(new_comp_depth == component_split_depth);

    comp_split.component_split_multiplier = comp_split_multiplier;
    component_split_model_count *= comp_split.component_split_multiplier;
    comp_split.model_count = node_model_count;
    assert(component_split_model_count == comp_split.model_count);
  }
  /**
   * Constructs the final set of partial assignments.
   *
   * @param input_set Input top tree which may be the complete top sample tree or only a subtree
   * from that tree.
   * @param output_set Output set of partial assignments.
   * @param component_split_depth Current component split depth.
   */
  static void BuildFinalTopTreeSet(const TopTreeNode &input_set,
                                   TopTreeNode &output_set,
                                   VariableIndex component_split_depth) {
    SampleSize num_samples = 0;
    output_set.descendants.clear();
    for (const auto &input : input_set.descendants) {
      if (input.num_samples == 0)
        continue;
      // No special processing required for a standard leaf.
      if (input.node_type_ == TopTreeSampleNodeType::MAX_DEPTH_NODE
          || input.node_type_ == TopTreeSampleNodeType::CYLINDER_NODE) {
        output_set.descendants.emplace_back(input);
      // Recursively handle component splits.
      } else if (input.node_type_ == TopTreeSampleNodeType::COMPONENT_SPLIT) {
        TopTreeNode comp_split(TopTreeSampleNodeType::COMPONENT_SPLIT);
        BuildFinalComponentSplitPartialAssignment(input, comp_split,
                                                  component_split_depth + 1);
        output_set.addDescendants(comp_split);
      } else {
        ExitWithError("Unknown top tree node type.");
      }
      num_samples += input.num_samples;
    }
    output_set.num_samples = num_samples;
    assert(verifySampleSizes(output_set));
    assert(component_split_depth != TOP_COMPONENT_SPLIT_DEPTH_
           || num_samples == config_->num_samples_);
  }
  /**
   * Recursive function used for forming the component split partial assignments.
   *
   * @param input_node Component split node in the partial assignment tree.
   * @param output_set Set of partial assignments for this top tree node.
   * @param component_split_depth Current component split depth in the top tree.
   */
  static void BuildFinalComponentSplitPartialAssignment(const TopTreeNode &input_node,
                                                        TopTreeNode &output_set,
                                                        VariableIndex component_split_depth) {
    assert(input_node.node_type_ == TopTreeSampleNodeType::COMPONENT_SPLIT);

    // Create a new component split object and give each the literal stack of the component split.
    output_set.num_samples = input_node.num_samples;
    output_set.createEmptyDescendants(output_set.num_samples);
    for (auto &desc : output_set.descendants)
      desc.addAssignedLitsAndEmancipatedVars(input_node);

    for (const auto &comp : input_node.descendants) {
      TopTreeNode comp_out(TopTreeSampleNodeType::COMPONENT);
      BuildFinalTopTreeSet(comp, comp_out, component_split_depth);
      output_set.randomComponentMerge(comp_out);
    }

    if (component_split_depth == TOP_COMPONENT_SPLIT_DEPTH_ + 1) {
      CollapseComponentSplitPartialAssignment(output_set);
      assert(verifySampleSizes(output_set));
    }
  }
  /**
   * Converts the partial assignment multiset into a partial assignment with an integer.
   *
   * @param output_set Output set of partial assignments.  It is modified in place if applicable.
   */
  static void CollapseComponentSplitPartialAssignment(TopTreeNode &output_set) {
    // ToDo Written for quick development and debug. Can be made more efficient.
    std::sort(output_set.descendants.begin(), output_set.descendants.end(),
              TopTreeNode::compare_partial_assn);
    SampleSize i = 0;
    std::vector<TopTreeNode> &desc = output_set.descendants;
    while (i < desc.size() - 1) {
      if (desc[i].HasSamePartialAssignment(desc[i + 1])) {
        assert(desc[i].node_type_ == desc[i + 1].node_type_);

        desc[i].num_samples += desc[i + 1].num_samples;
        desc.erase(desc.begin() + i + 1);  // ToDo change to a faster way than delete.
      } else {
        i++;
      }
    }
  }
  /**
   * Verifies the sample size of the top tree node matches that of all descendant nodes.
   *
   * @param tree_node SampleTopTree node which may be the root of the tree or some inner or leaf
   * vertex.
   * @return True if the sample sizes of the pass node match those of descendants.
   */
  inline static bool verifySampleSizes(const TopTreeNode &tree_node) {
    SampleSize num_samples = 0;
    for (const auto & descendant : tree_node.descendants) {
      num_samples += descendant.num_samples;
      if (descendant.node_type_ == TopTreeSampleNodeType::COMPONENT_SPLIT) {
        for (const auto & component : descendant.descendants) {
          if (descendant.num_samples != component.num_samples || !verifySampleSizes(component))
            return false;
        }
      }
    }
    return num_samples == tree_node.num_samples;
  }
  /**
   * Helper function used to write the final set of top tree samples to a file.
   *
   * Each line in the file represents a single top tree formula.  The format of each line is:
   *
   * <NumbSamplesFromFormula> <CommaSeparatedSignedSetLiterals>
   *
   * The set literals in this context represent the constrained (i.e., set) variables.  If variable
   * #150 is set to "True", it will be written as "150" in the comma separated list.  If it is set
   * to "False", it will be written as "-150".  This is to mimic the DIMACS file format.
   *
   * @param final_sample_sets Set of partial formulas and the associated final samples associated
   * with them.
   */
  static void WriteFinalTopTreeSampleSet(TopTreeNode &final_samples) {
    SampleSize num_completed_samples = 0, num_unfinished_samples = 0;
    for (const auto &sample : final_samples.descendants) {
      TopTreeSampleNodeType node_type = sample.node_type_;
      assert(sample.num_samples > 0 && (node_type == TopTreeSampleNodeType::CYLINDER_NODE
                                        || node_type == TopTreeSampleNodeType::MAX_DEPTH_NODE));

      if (sample.num_samples > config_->max_top_tree_leaf_sample_count
          && node_type != TopTreeSampleNodeType::CYLINDER_NODE) {
        num_unfinished_samples += sample.num_samples;
        remaining_formulas.addNew(sample);
      } else {
        num_completed_samples += sample.num_samples;
        output_partials.addNew(sample);
      }
    }
    assert(num_completed_samples + num_unfinished_samples == config_->num_samples_);
  }
  /**
   * Writes the specified output assignment to a file for processing.
   *
   * @param fout Output file to write.
   * @param sample Output sample including number of samples and assigned literals.
   */
  static void WriteOutputAssignment(std::ofstream &fout, const TopTreeNode &sample) {
    assert(sample.num_samples > 0 && sample.descendants.empty());

    // Literal stack - always at least one.
    fout << sample.num_samples << TOP_TREE_DATA_TYPE_SEP;
    for (auto set_lit : sample.assigned_lits_) {
      if (set_lit != sample.assigned_lits_[0])
        fout << TOP_TREE_FILE_ITEM_SEP;
      fout << set_lit;
    }
    fout << TOP_TREE_DATA_TYPE_SEP;
    // There may be no emancipated variables.
    if (!sample.emancipated_vars_.empty()) {
      fout << TOP_TREE_EMPTY_LIT_STACK;
    } else {
      for (auto var : sample.emancipated_vars_) {
        if (var != sample.emancipated_vars_.front())
          fout << TOP_TREE_FILE_ITEM_SEP;
        fout << var;
      }
    }
    fout << std::flush;
  }
  /**
   * Given a top tree set of size larger than the desired number of samples, this function reduces
   * the set uniformly at random to exactly the desired number of samples.
   *
   * @param sample_set Rough set of top tree sample nodes.  This is modified to be the
   * final set of top tree sample nodes.
   */
  static void DownsampleRoughSet(TopTreeNode &sample_set) {
    if (sample_set.num_samples == 0)
      return;

    SetOversamplerScalar(sample_set.num_samples);

    // Reset the rough counts
    tot_num_rough_samples_ = 0;
    for (auto &tree_node : sample_set.descendants) {
      if (tree_node.model_count > 0) {
        tree_node.num_samples = GetRoughSampleCount(sample_set.model_count, tree_node.model_count,
                                                    sample_set.num_samples);
        tot_num_rough_samples_ += tree_node.num_samples;
      }
    }
    // If insufficient samples were selected, reset and start over
    if (tot_num_rough_samples_ < sample_set.num_samples)
      return DownsampleRoughSet(sample_set);

    std::vector<TreeNodeIndex> final_sample_node_ids;
    final_sample_node_ids.reserve(tot_num_rough_samples_);

    // Build the array to be used for the modified Fisher-Yates shuffle to build the fine
    // sample set.
    for (TreeNodeIndex node_index = 0; node_index < sample_set.descendants.size(); node_index++) {
      for (TreeNodeIndex j = 0; j < sample_set.descendants[node_index].num_samples; j++)
        final_sample_node_ids.emplace_back(node_index);
      // Reset and will be updated after downsampling.
      sample_set.descendants[node_index].num_samples = 0;
    }
    assert(tot_num_rough_samples_ == final_sample_node_ids.size()
           && tot_num_rough_samples_ > sample_set.num_samples);

    Random::DownsampleList<TreeNodeIndex, TreeNodeIndex>(sample_set.num_samples,
                                                         final_sample_node_ids);
    assert(final_sample_node_ids.size() == sample_set.num_samples);

    // Update the sample count for each of the tree sample sets.
    for (TreeNodeIndex sample_num = 0; sample_num < sample_set.num_samples; sample_num++) {
      TreeNodeIndex node_idx = final_sample_node_ids[sample_num];
      sample_set.descendants[node_idx].num_samples++;
    }
    // Recursively push the correct sample counts down the sample stack.
    for (auto &subsample : sample_set.descendants) {
      if (subsample.node_type_ != TopTreeSampleNodeType::COMPONENT_SPLIT)
        continue;
      for (auto &subcomponent : subsample.descendants) {
        subcomponent.num_samples = subsample.num_samples;
        DownsampleRoughSet(subcomponent);
      }
    }
  }
  /**
   * Calculates the over-specification scalar for building the initial set of samples.
   */
  static void SetOversamplerScalar(SampleSize num_samples) {
    if (num_samples <= 5) {
      sample_count_scalar_ = 5.0;
    } else if (num_samples <= 20) {
      sample_count_scalar_ = 3.0;
    } else if (num_samples <= 125) {
      sample_count_scalar_ = 2.0;
    } else if (num_samples <= 1100) {
      sample_count_scalar_ = 1.5;
    } else {
      sample_count_scalar_ = 1.25;
    }
  }
  /**
   * For a given top tree node with A models and a targeted number of samples |S|, this function
   * provides a rough sample count for that node that will be down-sampled later to get exactly
   * |S| models.
   *
   * @param total_model_count Total number of models.
   *
   * @param tree_node_model_count Number of models in the top tree node.
   *
   * @param target_sample_count Target for the number of samples to select
   *
   * @return Binomially distributed random variable with expected value E[scalar * |S| * a / T]
   */
  static SampleSize GetRoughSampleCount(mpz_class &total_model_count,
                                        mpz_class &tree_node_model_count,
                                        SampleSize target_sample_count) {
    assert(total_model_count >= tree_node_model_count);
    return Random::Mpz::binom(static_cast<SampleSize>(sample_count_scalar_ * target_sample_count),
                              total_model_count, tree_node_model_count);
  }
  /**
   * Current depth of component splits.
   *
   * @return Component split depth.  If there has been no component splits in the current literal
   * stack, it returns 0.
   */
  static VariableIndex component_split_depth() {
    return literal_stack_start_loc_.size();
  }
  /**
   * File used to write the top tree samples during the normal execution.
   * This is NOT the final set of saved samples.
   */
  static std::ofstream all_top_tree_fout_;
  /**
   * Filename where the top tree samples are written.
   */
  static std::string all_top_tree_filename_;
  /**
   * Active solver's configuration.
   */
  static SolverConfiguration * config_;
  /**
   * Active solver's statistics.
   */
  static DataAndStatistics * statistics_;
  /**
   * Root solver's configuration.
   */
  static SolverConfiguration * root_config_;
  /**
   * Root solver's statistics.
   */
  static DataAndStatistics * root_statistics_;
  /**
   * Total number of models observed across all top tree nodes.
   */
  static mpz_class total_models_;
  /**
   * Total number of samples after running the rough sample count.
   */
  static SampleSize tot_num_rough_samples_;
  /**
   * Debug only data structure used for tracing the literal stack building.
   *
   * ToDo Remove the debug literal stack.
   */
  static std::vector<LiteralID> debug_lit_stack_;
  /**
   * Starting point in the literal stack corresponding to the component split.
   */
  static std::vector<VariableIndex> literal_stack_start_loc_;
  /**
   * This setting specifies how much to over sample the initial set of top tree samples.
   */
  static float sample_count_scalar_;
  /**
   * Denotes what is used to represent no component split.
   */
  static const VariableIndex TOP_COMPONENT_SPLIT_DEPTH_;
};

#endif // TOP_TREE_SAMPLER_H
