/**
 * instance.cpp
 *
 * Purpose: Defines the methods for the Instance() class.
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


#include <sys/stat.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "instance.h"
#include "model_sampler.h"

void Instance::cleanClause(ClauseOfs cl_ofs) {
  bool satisfied = false;
  for (auto it = beginOf(cl_ofs); *it != SENTINEL_LIT; it++)
    if (isSatisfied(*it)) {
      satisfied = true;
      break;
    }
  // mark the clause as empty if satisfied
  if (satisfied) {
    *beginOf(cl_ofs) = SENTINEL_LIT;
    return;
  }
  auto jt = beginOf(cl_ofs);
  auto it = beginOf(cl_ofs);
  // from now, all inactive literals are resolved
  for (; *it != SENTINEL_LIT; it++, jt++) {
    while (*jt != SENTINEL_LIT && !isActive(*jt))
      jt++;
    *it = *jt;
    if (*jt == SENTINEL_LIT)
      break;
  }
  unsigned length = it - beginOf(cl_ofs);
  // if it has become a unit clause, it should have already been asserted
  if (length == 1) {
    *beginOf(cl_ofs) = SENTINEL_LIT;
    // if it has become binary, transform it to binary and delete it
  } else if (length == 2) {
    addBinaryClause(*beginOf(cl_ofs), *(beginOf(cl_ofs) + 1));
    *beginOf(cl_ofs) = SENTINEL_LIT;
  }
}

void Instance::compactClauses() {
  std::vector<ClauseOfs> clause_ofs;
  clause_ofs.reserve(statistics_.num_long_clauses_);

  // clear watch links and occurrence lists
  for (auto it_lit = literal_pool_.begin(); it_lit != literal_pool_.end();
      it_lit++) {
    if (*it_lit == SENTINEL_LIT) {
      if (it_lit + 1 == literal_pool_.end())
        break;
      it_lit += ClauseHeader::overheadInLits();
      clause_ofs.push_back(1 + it_lit - literal_pool_.begin());
    }
  }

  for (auto ofs : clause_ofs)
    cleanClause(ofs);

  for (auto &l : literals_)
    l.resetWatchList();

  occurrence_lists_.clear();
  occurrence_lists_.resize(variables_.size());

  std::vector<LiteralID> tmp_pool = literal_pool_;
  literal_pool_.clear();
  literal_pool_.push_back(SENTINEL_LIT);
  ClauseOfs new_ofs;
  unsigned num_clauses = 0;
  for (auto ofs : clause_ofs) {
    auto it = (tmp_pool.begin() + ofs);
    if (*it != SENTINEL_LIT) {
      for (unsigned i = 0; i < ClauseHeader::overheadInLits(); i++)
        literal_pool_.emplace_back(0);
      new_ofs = literal_pool_.size();
      literal(*it).addWatchLinkTo(new_ofs);
      literal(*(it + 1)).addWatchLinkTo(new_ofs);
      num_clauses++;
      for (; *it != SENTINEL_LIT; it++) {
        literal_pool_.push_back(*it);
        occurrence_lists_[*it].push_back(new_ofs);
      }
      literal_pool_.push_back(SENTINEL_LIT);
    }
  }

  std::vector<LiteralID> tmp_bin;
  unsigned bin_links = 0;
  for (auto &l : literals_) {
    tmp_bin.clear();
    for (auto it = l.binary_links_.begin(); *it != SENTINEL_LIT; it++)
      if (isActive(*it))
        tmp_bin.push_back(*it);
    bin_links += tmp_bin.size();
    tmp_bin.push_back(SENTINEL_LIT);
    l.binary_links_ = tmp_bin;
  }
  statistics_.num_long_clauses_ = num_clauses;
  statistics_.num_binary_clauses_ = bin_links >> 1;
}

void Instance::compactVariables() {
  std::vector<unsigned> var_map(variables_.size(), 0);
  unsigned last_ofs = 0;
  unsigned num_isolated = 0;  // Emancipated/freed variables no longer in any clause
  LiteralIndexedVector<std::vector<LiteralID> > _tmp_bin_links(1);
  LiteralIndexedVector<TriValue> _tmp_values = literal_values_;

  for (auto l : literals_)
    _tmp_bin_links.push_back(l.binary_links_);

  assert(_tmp_bin_links.size() == literals_.size());
  for (unsigned v = 1; v < variables_.size(); v++)
    if (isActive(v)) {
      if (isolated(v)) {
        std::cout << "Variable " << v << " is isolated." << std::endl;
        num_isolated++;
        continue;
      }
      last_ofs++;
      var_map[v] = last_ofs;
    }

  variables_.clear();
  variables_.resize(last_ofs + 1);
  occurrence_lists_.clear();
  occurrence_lists_.resize(variables_.size());
  literals_.clear();
  literals_.resize(variables_.size());
  literal_values_.clear();
  literal_values_.resize(variables_.size(), X_TRI);

  unsigned bin_links = 0;
  LiteralID newlit;
  for (auto l = LiteralID(0, false); l != _tmp_bin_links.end_lit(); l.inc()) {
    if (var_map[l.var()] != 0) {
      newlit = LiteralID(var_map[l.var()], l.sign());
      for (auto it = _tmp_bin_links[l].begin(); *it != SENTINEL_LIT; it++) {
        assert(var_map[it->var()] != 0);
        literals_[newlit].addBinLinkTo(
            LiteralID(var_map[it->var()], it->sign()));
      }
      bin_links += literals_[newlit].binary_links_.size() - 1;
    }
  }

  std::vector<ClauseOfs> clause_ofs;
  clause_ofs.reserve(statistics_.num_long_clauses_);
  // clear watch links and occurrence lists
  for (auto it_lit = literal_pool_.begin(); it_lit != literal_pool_.end();
      it_lit++) {
    if (*it_lit == SENTINEL_LIT) {
      if (it_lit + 1 == literal_pool_.end())
        break;
      it_lit += ClauseHeader::overheadInLits();
      clause_ofs.push_back(1 + it_lit - literal_pool_.begin());
    }
  }

  for (auto ofs : clause_ofs) {
    literal(LiteralID(var_map[beginOf(ofs)->var()], beginOf(ofs)->sign())).addWatchLinkTo(
        ofs);
    literal(LiteralID(var_map[(beginOf(ofs) + 1)->var()],
            (beginOf(ofs) + 1)->sign())).addWatchLinkTo(ofs);
    for (auto it_lit = beginOf(ofs); *it_lit != SENTINEL_LIT; it_lit++) {
      *it_lit = LiteralID(var_map[it_lit->var()], it_lit->sign());
      occurrence_lists_[*it_lit].push_back(ofs);
    }
  }

  literal_values_.clear();
  literal_values_.resize(variables_.size(), X_TRI);
  unit_clauses_.clear();

  statistics_.num_variables_ = variables_.size() - 1 + num_isolated;

  statistics_.num_used_variables_ = num_variables();
  statistics_.num_free_variables_ = num_isolated;
}

void Instance::compactConflictLiteralPool() {
  auto write_pos = conflict_clauses_begin();
  std::vector<ClauseOfs> tmp_conflict_clauses = conflict_clauses_;
  conflict_clauses_.clear();
  for (auto clause_ofs : tmp_conflict_clauses) {
    auto read_pos = beginOf(clause_ofs) - ClauseHeader::overheadInLits();
    for (unsigned i = 0; i < ClauseHeader::overheadInLits(); i++)
      *(write_pos++) = *(read_pos++);
    ClauseOfs new_ofs =  write_pos - literal_pool_.begin();
    conflict_clauses_.push_back(new_ofs);
    // first substitute antecedent if clause_ofs implied something
    if (isAntecedentOf(clause_ofs, *beginOf(clause_ofs)))
      var(*beginOf(clause_ofs)).ante = Antecedent(new_ofs);

    // now redo the watches
    literal(*beginOf(clause_ofs)).replaceWatchLinkTo(clause_ofs, new_ofs);
    literal(*(beginOf(clause_ofs)+1)).replaceWatchLinkTo(clause_ofs, new_ofs);
    // next, copy clause data
    assert(read_pos == beginOf(clause_ofs));
    while (*read_pos != SENTINEL_LIT)
      *(write_pos++) = *(read_pos++);
    *(write_pos++) = SENTINEL_LIT;
  }
  literal_pool_.erase(write_pos, literal_pool_.end());
}


//bool Instance::deleteConflictClauses() {
//  statistics_.times_conflict_clauses_cleaned_++;
//  std::vector<ClauseOfs> tmp_conflict_clauses = conflict_clauses_;
//  conflict_clauses_.clear();
//  std::vector<double> tmp_ratios;
//  double score, lifetime;
//  for (auto clause_ofs: tmp_conflict_clauses) {
//    score = getHeaderOf(clause_ofs).score();
//    lifetime = statistics_.num_conflicts_ - getHeaderOf(clause_ofs).creation_time();
//    tmp_ratios.push_back(score/lifetime/(getHeaderOf(clause_ofs).length()));
//  }
//  std::vector<double> tmp_ratiosB = tmp_ratios;
//
//  sort(tmp_ratiosB.begin(), tmp_ratiosB.end());
//
//  double cutoff = tmp_ratiosB[tmp_ratiosB.size()/2];
//
//  for (unsigned i = 0; i < tmp_conflict_clauses.size(); i++) {
//    if (tmp_ratios[i] < cutoff) {
//      if (!markClauseDeleted(tmp_conflict_clauses[i]))
//        conflict_clauses_.push_back(tmp_conflict_clauses[i]);
//    } else
//      conflict_clauses_.push_back(tmp_conflict_clauses[i]);
//  }
//  return true;
//}

bool Instance::deleteConflictClauses(bool delete_all) {
  statistics_.times_conflict_clauses_cleaned_++;
  std::vector<ClauseOfs> tmp_conflict_clauses = conflict_clauses_;
  conflict_clauses_.clear();
  std::vector<double> tmp_ratios;
  for (auto clause_ofs : tmp_conflict_clauses) {
    double score = getHeaderOf(clause_ofs).score();
//    lifetime = statistics_.num_conflicts_ - getHeaderOf(clause_ofs).creation_time();
//    tmp_ratios.push_back(score/lifetime);
    tmp_ratios.push_back(score);
  }
  std::vector<double> tmp_ratiosB = tmp_ratios;

  sort(tmp_ratiosB.begin(), tmp_ratiosB.end());

  double cutoff = -1;
  if (!tmp_ratiosB.empty())
    cutoff = tmp_ratiosB[tmp_ratiosB.size()/2];

  for (unsigned i = 0; i < tmp_conflict_clauses.size(); i++) {
    if (delete_all) {
      markClauseDeleted(tmp_conflict_clauses[i]);
      continue;
    }
    if (tmp_ratios[i] < cutoff) {
      if (!markClauseDeleted(tmp_conflict_clauses[i]))
        conflict_clauses_.push_back(tmp_conflict_clauses[i]);
    } else {
      conflict_clauses_.push_back(tmp_conflict_clauses[i]);
    }
  }
  return true;
}


bool Instance::markClauseDeleted(ClauseOfs cl_ofs) {
  // only first literal may possibly have cl_ofs as antecedent
  if (isAntecedentOf(cl_ofs, *beginOf(cl_ofs)))
    return false;

  literal(*beginOf(cl_ofs)).removeWatchLinkTo(cl_ofs);
  literal(*(beginOf(cl_ofs)+1)).removeWatchLinkTo(cl_ofs);
  return true;
}


bool Instance::createfromFile(const std::string &file_name) {
  if (!file_name.empty())
    statistics_.input_file_ = file_name;
  unsigned max_ignore = 1000000;  // Max number of characters on line to ignore.

  // Initialize the literal pool
  literal_pool_.clear();
  literal_pool_.push_back(SENTINEL_LIT);

  // Builds the list
  variables_.clear();
  variables_.emplace_back();  //initializing the Sentinel
  literal_values_.clear();
  unit_clauses_.clear();
  unused_vars_.clear();

  // BEGIN File input
  std::ifstream input_file(statistics_.input_file_);
  if (!input_file)
    ExitWithError("Cannot open file: " + statistics_.input_file_, EX_IOERR);

  char c;
  while ((input_file >> c) && c != 'p')
    ParseCnfCommentForSupport(input_file);
  std::string idstring;
  long nVars = -1, nCls = -1;
  if (!((input_file >> idstring) && idstring == "cnf" && (input_file >> nVars)
      && (input_file >> nCls)))
    ExitWithError("Invalid CNF file", EX_PROTOCOL);
  else if (nVars <= 0)
    ExitWithError("Invalid variable count.  At least one variable is required.", EX_PROTOCOL);
  else if (nCls < 0)
    ExitWithError("Invalid clause count.  A negative clause count is invalid.", EX_PROTOCOL);

  variables_.resize(nVars + FIRST_VAR);
  literal_values_.resize(nVars + FIRST_VAR, X_TRI);
  struct stat file_status;
  stat(statistics_.input_file_ .c_str(), &file_status);
  literal_pool_.reserve(file_status.st_size);
  conflict_clauses_.reserve(2 * nCls);
  occurrence_lists_.clear();
  occurrence_lists_.resize(nVars + FIRST_VAR);

  std::vector<LiteralID> literals;
  literals_.clear();
  literals_.resize(nVars + FIRST_VAR);
  std::vector<bool> var_used(nVars + FIRST_VAR, false);

  // Analyze each clause and add literals/variables as appropriate.
  unsigned clauses_added = 0;
  while ((input_file >> c) && clauses_added < nCls) {
    input_file.unget();  // extracted a nonspace character to determine if
                         // we have a clause, so put it back

    if ((c == '-') || isdigit(c)) {
      literals.clear();  // Empty the literal set in the clause
      bool skip_clause = false;
      long lit;
      while ((input_file >> lit) && lit != 0) {
        var_used[abs(lit)] = true;  // Mark the variable as used.
        bool duplicate_literal = false;
        for (auto i : literals) {
          // Checks if same literal already in the clause
          if (i.toInt() == lit) {
            duplicate_literal = true;
            break;
          }
          // Checks if a clause has a literal and its complement in the same clause.
          if (i.toInt() == -lit) {
            skip_clause = true;
            break;
          }
        }
        if (!duplicate_literal)
          literals.emplace_back(lit);
      }
      if (!skip_clause) {
        assert(!literals.empty());  // May report a fail in an UNSAT file.
        clauses_added++;
        statistics_.incorporateClauseData(literals);
        long cl_ofs = addClause(literals);
        if (cl_ofs == CLAUSE_ADDING_ERROR)
          return false;
        // If the clause is non-binary and non-unit, then
        // Add to the occurrence list for the literal.
        if (literals.size() >= 3)
          for (auto l : literals)
            occurrence_lists_[l].push_back(cl_ofs);
      }
    } else if (c == 'c') {
      input_file >> c;
      ParseCnfCommentForSupport(input_file);
      continue;
    }
    input_file.ignore(max_ignore, '\n');
  }
  // END NEW
  input_file.close();
  //  /// END FILE input
  statistics_.num_variables_ = statistics_.num_original_variables_ = (VariableIndex)nVars;
  statistics_.num_used_variables_ = num_variables();
  statistics_.num_indep_support_variables_ = independent_support.size();
//  statistics_.num_free_variables_ = nVars - num_variables();  // This line seems not to work right
  if (!independent_support.empty()) {
    has_independent_support_ = true;
    std::sort(independent_support.begin(), independent_support.end());
  }

  // Get a list of unused variables.
  statistics_.num_free_variables_ = 0;
  for (VariableIndex i = 1; i < var_used.size(); i++) {
    if (!var_used[i]) {
      statistics_.num_free_variables_++;
//      unused_vars_.push_back(i);  // Done instead at the top of var stack
    }
  }
  statistics_.num_original_clauses_ = static_cast<ClauseIndex>(nCls);

  statistics_.num_original_binary_clauses_ = statistics_.num_binary_clauses_;
  statistics_.num_original_unit_clauses_ = statistics_.num_unit_clauses_ = unit_clauses_.size();

  original_lit_pool_size_ = literal_pool_.size();
  return true;
}


void Instance::ParseCnfCommentForSupport(std::ifstream &input_file) {
  std::string ind_str, line;
  std::istringstream iss;
  // Check if the comment line contains an independent support
  std::getline(input_file, line);
  iss.str(line);
  iss.clear();
  iss >> ind_str;
  if (ind_str != "ind")
    return;
  VariableIndex support_var;
  while ((iss >> support_var) && support_var != 0)
    independent_support.emplace_back(support_var);
}

