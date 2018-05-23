/**
 * instance.h
 *
 * Purpose: Defines the "Instance()" class that is the superclass of the main "Solver()" class.
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

#ifndef INSTANCE_H_
#define INSTANCE_H_

#include <sysexits.h>
#include <limits.h>

#include <cassert>
#include <string>
#include <vector>

#include "statistics.h"
#include "structures.h"
#include "containers.h"


#define CLAUSE_ADDING_ERROR INT_MIN

class Instance {
 protected:
  void unSet(LiteralID lit) {
    var(lit).ante = Antecedent(NOT_A_CLAUSE);
    var(lit).decision_level = INVALID_DL;
    literal_values_[lit] = X_TRI;
    literal_values_[lit.neg()] = X_TRI;
  }

  Antecedent & getAntecedent(LiteralID lit) {
    return variables_[lit.var()].ante;
  }
  /**
   * Accessor the antecedent object.  It is a constant alternative
   * to the method getAntecedent(LiteralID).
   *
   * @param lit Literal identification number
   *
   * @return Variable corresponding to the literal.
   */
  const Antecedent & antecedent(LiteralID lit) const {
    return variables_[lit.var()].ante;
  }

  /**
   * Checks if the literal has an antecedent clause.
   *
   * @param lit Literal identification object
   *
   * @return True if the literal has an antecedent.
   */
  bool hasAntecedent(LiteralID lit) const {
    return variables_[lit.var()].ante.isAnt();
  }

  bool isAntecedentOf(ClauseOfs ante_cl, LiteralID lit) {
    return var(lit).ante.isAClause() && (var(lit).ante.asCl() == ante_cl);
  }

  bool isolated(VariableIndex v) {
    LiteralID lit(v, false);
    return (literal(lit).binary_links_.size() <= 1)
        & occurrence_lists_[lit].empty()
        & (literal(lit.neg()).binary_links_.size() <= 1)
        & occurrence_lists_[lit.neg()].empty();
  }

  bool free(VariableIndex v) {
    return isolated(v) & isActive(v);
  }

  bool deleteConflictClauses(bool delete_all = false);
  bool markClauseDeleted(ClauseOfs cl_ofs);

  // Compact the literal pool erasing all the clause
  // information from deleted clauses
  void compactConflictLiteralPool();

  // we assert that the formula is consistent
  // and has not been found UNSAT yet
  // hard wires all assertions in the literal stack into the formula
  // removes all set variables and essentially reinitiallizes all
  // further data
  void compactClauses();
  void compactVariables();
  void cleanClause(ClauseOfs cl_ofs);

  /////////////////////////////////////////////////////////
  // END access to variables and literals
  /////////////////////////////////////////////////////////


  unsigned long num_conflict_clauses() const {
    return conflict_clauses_.size();
  }
  /**
   * Number of variables in the formula.
   *
   * @return Size of the variables array
   */
  const VariableIndex num_variables() const {
    return (VariableIndex)variables_.size() - 1;
  }
  /**
   * Reads the CNF formula from a file.  It does minimal analysis of the
   * clauses.  It adds the clauses to the basic clause data structures
   * including the occurence lists.
   *
   * If there is a self contradictory clause, it may catch that
   * failure but it is not guaranted.
   *
   * @param file_name Location of the CNF formula in a file
   *
   * @return "true" if the file was parsed and the formula properly built.
   */
  bool createfromFile(const std::string &file_name);

  DataAndStatistics statistics_;

  /** literal_pool_: the literals of all clauses are stored here
   *   INVARIANT: first and last entries of literal_pool_ are a SENTINEL_LIT
   *
   *   Clauses begin with a ClauseHeader structure followed by the literals
   *   terminated by SENTINEL_LIT
   */
  std::vector<LiteralID> literal_pool_;

  /**
   * this is to determine the starting offset of
   * conflict clauses - MT
   *
   * This contains the number of unique literals that appear in the original
   * formula.  It will be less than or equal to 2 * |variables_|
   */
  unsigned long original_lit_pool_size_;

  LiteralIndexedVector<Literal> literals_;

  /**
   * Set of long clauses that the specified literal appears in.
   */
  LiteralIndexedVector<std::vector<ClauseOfs>> occurrence_lists_;

  std::vector<ClauseOfs> conflict_clauses_;

  /**
   * Set of all unit clauses (i.e., list of unit literals)
   */
  std::vector<LiteralID> unit_clauses_;

  std::vector<Variable> variables_;
  LiteralIndexedVector<TriValue> literal_values_;


  /*-------------------------------------------
  --         Begin Sampler Objects           --
  -------------------------------------------*/
  /**
   * Stores the variable identification numbers of
   * any variables that do not appear in a
   * formula at all.
   */
  std::vector<VariableIndex> unused_vars_;
  /*-------------------------------------------
  --         End Sampler Objects           --
  -------------------------------------------*/

  /**
   * Decay the activity scores of the literal based on a regular decision
   * frequency (128 decisions in the original code).
   */
  void decayActivities() {
    for (auto l_it = literals_.begin(); l_it != literals_.end(); l_it++)
      l_it->activity_score_ *= 0.5;

    for (auto clause_ofs : conflict_clauses_)
        getHeaderOf(clause_ofs).decayScore();
  }

  void updateActivities(ClauseOfs clause_ofs) {
    getHeaderOf(clause_ofs).increaseScore();
    for (auto it = beginOf(clause_ofs); *it != SENTINEL_LIT; it++) {
      literal(*it).increaseActivity();
    }
  }

  /**
   * Checks whether the specified literal is unit clause.
   *
   * @param lit Literal identification object.
   *
   * @return True if the literal is a unit clause.
   */
  bool isUnitClause(const LiteralID lit) {
    for (auto l : unit_clauses_)
      if (l == lit)
        return true;
    return false;
  }
  /**
   * Checks whether a unit clause exists for the specified variable.
   *
   * @param v Variable identification number
   *
   * @return true if the specified variable is in a unit clause.
   */
  bool existsUnitClauseOf(VariableIndex v) {
    for (auto l : unit_clauses_)
      // Extract the unit clause's variable number
      if (l.var() == v)
        return true;
    return false;
  }
  /**
   * Parses a comment line in a DIMACs file to extract the independent support if it is
   * present.
   *
   * @param input_file Input file stream of a DIMACs file.
   */
  void ParseCnfCommentForSupport(std::ifstream &input_file);

  // addUnitClause checks whether lit or lit.neg() is already a
  // unit clause
  // a negative return value implied that the Instance is UNSAT
  /**
   * Attempts to add a unit clause literal to the set of unit clauses.
   *
   * If the unit clause (or its complement) already exists, the function
   * merely returns and does nothing.
   *
   * @param lit Unit clause literal
   *
   * @return True if the literal was added or already exists in the unit clause list
   * A return of false implies the formula is UNSAT.
   */
  bool addUnitClause(const LiteralID lit) {
    for (auto l : unit_clauses_) {
      if (l == lit)
        return true;
      if (l == lit.neg())
        return false;
    }
    unit_clauses_.push_back(lit);
    return true;
  }
  /**
   * Adds a clause (be it unit, binary, or long) to the set of clauses
   * for the formula.
   *
   * @param literals Set of literals in the clause
   * @return Identification number of the clause.  If there is an error adding the
   * clause, it returns "INT_MIN" (i.e., the minimum integer value).
   */
  inline long addClause(std::vector<LiteralID> &literals);

  // adds a UIP Conflict Clause
  // and returns it as an Antecedent to the first
  // literal stored in literals
  inline Antecedent addUIPConflictClause(std::vector<LiteralID> &literals);

  inline bool addBinaryClause(LiteralID litA, LiteralID litB);

  /////////////////////////////////////////////////////////
  // BEGIN access to variables, literals, clauses
  /////////////////////////////////////////////////////////

  /**
   * Accessor to get a variable object from a literal.
   *
   * @param lit Literal identification object.
   *
   * @return Variable object associated with the literal
   */
  inline Variable& var(const LiteralID lit) {
    //return variables_[lit.var()];
    return const_cast<Variable&>(var_const(lit));
  }
  /**
   * Const Version of the Variable Function.
   *
   * @param lit Literal identification object.
   *
   * @return Variable object associated with the literal
   */
  inline const Variable& var_const(const LiteralID lit) const {
    return variables_[lit.var()];
  }
  /**
   * Accessor to extract the specified literal.
   *
   * @param lit Identification number for a positive or negative literal.
   *
   * @return Literal with the specified ID
   */
  Literal & literal(LiteralID lit) {
    return literals_[lit];
  }

  inline bool isSatisfied(const LiteralID &lit) const {
    return literal_values_[lit] == T_TRI;
  }

  bool isResolved(LiteralID lit) {
    return literal_values_[lit] == F_TRI;
  }
  /**
   *
   * @param lit Literal ID number
   * @return true if the variable is unset.
   */
  bool isActive(LiteralID lit) const {
    return literal_values_[lit] == X_TRI;
  }

  std::vector<LiteralID>::const_iterator beginOf(ClauseOfs cl_ofs) const {
    return literal_pool_.begin() + cl_ofs;
  }
  std::vector<LiteralID>::iterator beginOf(ClauseOfs cl_ofs) {
    return literal_pool_.begin() + cl_ofs;
  }

  decltype(literal_pool_.begin()) conflict_clauses_begin() {
    return literal_pool_.begin() + original_lit_pool_size_;
  }

  ClauseHeader &getHeaderOf(ClauseOfs cl_ofs) {
    return *reinterpret_cast<ClauseHeader *>(&literal_pool_[cl_ofs
                                                            - ClauseHeader::overheadInLits()]);
  }

  bool isSatisfied(ClauseOfs cl_ofs) {
    for (auto lt = beginOf(cl_ofs); *lt != SENTINEL_LIT; lt++)
      if (isSatisfied(*lt))
        return true;
    return false;
  }
  /**
   * Stores whether the solver has an independent support.  This is based on the contents of the
   * CNF file.
   */
  bool has_independent_support_ = false;
  /**
   * Variables that represent the independent support.
   */
  std::vector<VariableIndex> independent_support;
};

long Instance::addClause(std::vector<LiteralID> &literals) {
  if (literals.size() == 1) {
    // MT_TODO Deal properly with the situation that opposing unit clauses are learned
//    assert(!isUnitClause(literals[0].neg()));
    if (isUnitClause(literals[0].neg()))
      return CLAUSE_ADDING_ERROR;
    unit_clauses_.push_back(literals[0]);
    return 0;
  }
  if (literals.size() == 2) {
    addBinaryClause(literals[0], literals[1]);
    return 0;
  }
  for (unsigned i = 0; i < ClauseHeader::overheadInLits(); i++)
    literal_pool_.emplace_back(0);
  ClauseOfs cl_ofs = literal_pool_.size();

  for (auto l : literals) {
    literal_pool_.push_back(l);
    literal(l).increaseActivity(1);
  }
  // make an end: SENTINEL_LIT
  literal_pool_.push_back(SENTINEL_LIT);
  literal(literals[0]).addWatchLinkTo(cl_ofs);
  literal(literals[1]).addWatchLinkTo(cl_ofs);
  getHeaderOf(cl_ofs).set_creation_time(statistics_.num_conflicts_);
  return cl_ofs;
}


Antecedent Instance::addUIPConflictClause(std::vector<LiteralID> &literals) {
    Antecedent ante(NOT_A_CLAUSE);
    statistics_.num_clauses_learned_++;
    ClauseOfs cl_ofs = addClause(literals);
    if (cl_ofs != 0) {
      conflict_clauses_.push_back(cl_ofs);
      getHeaderOf(cl_ofs).set_length(literals.size());
      ante = Antecedent(cl_ofs);
    } else if (literals.size() == 2) {
      ante = Antecedent(literals.back());
      statistics_.num_binary_conflict_clauses_++;
    } else if (literals.size() == 1) {
      statistics_.num_unit_clauses_++;
    }
    return ante;
  }

bool Instance::addBinaryClause(LiteralID litA, LiteralID litB) {
  if (literal(litA).hasBinaryLinkTo(litB))
    return false;
  literal(litA).addBinLinkTo(litB);
  literal(litB).addBinLinkTo(litA);
  literal(litA).increaseActivity();
  literal(litB).increaseActivity();
  return true;
}


#endif /* INSTANCE_H_ */
