/**
 * structures.h
 *
 * Purpose: Defines a set of classes originally created by Marc Thurley include LiteralID, Literal,
 * Antecedent, Variable, and ClauseHeader.
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

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <assert.h>
#include <vector>
#include <iostream>

#include "primitive_types.h"

#define INVALID_DL -1

typedef unsigned char TriValue;
/// Macro value for an unsatisfied LITERAL
#define   F_TRI  0
/// Macro value for a satisfied LITERAL
#define   T_TRI  1
/// Macro value for an unassigned/free LITERAL
#define   X_TRI  2

class LiteralID {
 public:
  /**
   * Initializes the literal to negated and ID 0.
   */
  LiteralID() {
    value_ = 0;
  }
  /**
   * Converts a literal number into a literal ID.  The absolute value of the passed literal ID
   * is the variable number.  The sign represents whether the literal is negated or non-negated.
   *
   * @param lit Literal ID with a negative number representing a negated literal
   */
  // This code does screwy casting.  Cannot use explicit
  LiteralID(int lit) {  // NOLINT (runtime/explicit)
    value_ = (abs(lit) << 1) + static_cast<unsigned>(lit > 0);
  }
  /**
   * Build a literal from the variable ID and a sign term.
   *
   * @param var Variable identification number.
   * @param sign true for a positive literal and false for a negative literal
   */
  LiteralID(VariableIndex var, bool sign) {
    value_ = (var << 1) + static_cast<unsigned>(sign);
  }

  /**
   * Gets the variable identification number associated with the literal.  It ignores the
   * literal's sign (i.e., negation/positive)
   *
   * @return Variable ID number
   */
  VariableIndex var() const {
    return (value_ >> 1);
  }
  /**
   * Converts the literal into the form (sign)VariableNumber
   * where sign is +/- depending whether the literal
   * is negated.
   *
   * @return Signed version of the variable number.
   */
  int toInt() const {
    return (static_cast<int>(value_) >> 1) * ((sign()) ? 1 : -1);
  }
  /**
   * Increments the literal ID number.  This may cause a negated
   * literal to become positive or cause a positive literal's
   * ID to become the next negated literal.
   */
  void inc() {++value_;}
  /**
   * Sets the literal ID to the processed value.  It does
   * not process its integrity or correctness.
   *
   * @param v New literal value
   */
  void copyRaw(unsigned int v) {
    value_ = v;
  }
  /**
   * Returns whether the literal is positive or negated.
   *
   * @return True if the literal is positive (i.e., unnegated)
   */
  const bool sign() const {
    return static_cast<bool>(value_ & 0x01);
  }
  const bool operator!=(const LiteralID &rL2) const {
    return value_ != rL2.value_;
  }
  const bool operator==(const LiteralID &rL2) const {
    return value_ == rL2.value_;
  }
  const bool operator<(const LiteralID &rL2) const {
    return value_ < rL2.value_;
  }
  /**
   * Creates a literal that has the opposite sign as the
   * current literal.  For example, if the implicit literal is
   * negated, this returns a positive/unnegated literal.  If it is positive,
   * this function returns a negated literal.
   *
   * @return LiteralID of this literal just negated.
   */
  const LiteralID neg() const {
    return LiteralID(var(), !sign());
  }

  void print() const {
    std::cout << (sign() ? " " : "-") << var() << " ";
  }
  /**
   * Accesses the raw value of how the literal is stored in the data
   * structure.
   *
   * @return Unmodified verison of the data structure value
   */
  inline unsigned raw() const { return value_;}

 private:
  /**
   * Encodes the literal value.  Bit 0 is the sign with 0b1 being
   * a positive literal while 0b0 is a negative literal.
   *
   * The literal ID number is shifted by 1 bit to the left.
   */
  unsigned value_;

  template <class _T> friend class LiteralIndexedVector;
};

static const LiteralID NOT_A_LIT(0, false);
#define SENTINEL_LIT NOT_A_LIT

/**
 * Reference class storing information on each literal (i.e., both positive
 * and negative0.  The primary information this tracks is:
 * <ul>
 *  <li>binary_links_ - Represents variables linked by binary clauses</li>
 *  <li>Watch List - Clauses associated with the two watch literals principle</li>
 * </ul>
 */
class Literal {
 public:
  std::vector<LiteralID> binary_links_ = std::vector<LiteralID>(1, SENTINEL_LIT);
  /// List of watch clauses for the two watch literals paradigm.
  std::vector<ClauseOfs> watch_list_ = std::vector<ClauseOfs>(1, SENTINEL_CL);  // Watch clauses
  float activity_score_ = 0.0f;

  /**
   * Increase the literal's supplied literals activity score by the specified amount.
   *
   * @param u Amount to increase the literals activity score by.  This should generally be
   * a positive number.
   */
  void increaseActivity(unsigned u = 1) {
    activity_score_+= u;
  }
  /**
   * Resets the literals activity score which represents how frequently it appears in
   * a clause.
   */
  void resetActivity() { activity_score_ = 0;}
  /**
   *
   * @param clause_ofs
   */
  void removeWatchLinkTo(ClauseOfs clause_ofs) {
    for (auto it = watch_list_.begin(); it != watch_list_.end(); it++)
      if (*it == clause_ofs) {
        *it = watch_list_.back();  // Get last element in vector and replace item to be removed
        watch_list_.pop_back();  // Delete the last element
        return;
      }
  }

  void replaceWatchLinkTo(ClauseOfs clause_ofs, ClauseOfs replace_ofs) {
    for (auto it = watch_list_.begin(); it != watch_list_.end(); it++)
      if (*it == clause_ofs) {
        *it = replace_ofs;
        return;
      }
  }
  /**
   * Associates the clause to be monitored with this literal.
   *
   * @param clause_ofs Clause to be monitored
   */
  void addWatchLinkTo(ClauseIndex clause_ofs) {
    watch_list_.push_back(clause_ofs);
  }
  /**
   * Creates a link between the implicit literal and the explicitly supplied literal.
   * This implies the two literals are in the same binary clause.
   *
   * @param lit Other literal to be linked to.
   */
  void addBinLinkTo(LiteralID lit) {
    binary_links_.back() = lit;
    binary_links_.push_back(SENTINEL_LIT);
  }
  /**
   * Resets the binary watch list for a literal.  It removes all binary links from
   * the literal to all binary clauses.
   *
   * It does put the sentinel clause onto the stack.
   */
  void resetWatchList() {
    watch_list_.clear();
    watch_list_.push_back(SENTINEL_CL);
  }
  /**
   * Checks if the implicit literal has a link to the explicitly supplied literal.
   *
   * @param lit Other literal to check for a binary clause linkage to.
   *
   * @return true if the two literals appear in the same binary clause.
   */
  bool hasBinaryLinkTo(LiteralID lit) {
    for (auto l : binary_links_) {
      if (l == lit)
        return true;
    }
    return false;
  }
  /**
   * Checks if the literal is associated with any binary links.
   *
   * @return True if there is at least one binary link.
   */
  bool hasBinaryLinks() {
    return !binary_links_.empty();
  }
};

/**
 * Antecedent can be either a clause or a literal.
 */
class Antecedent {
  unsigned int val_;

 public:
  /**
   * Creates an empty antecedent clause with ID zero.
   */
  Antecedent() {
    val_ = 1;
  }
  /**
   * Shifts the clause offset by 1 to the left and then adds 1.
   *
   * @param cl_ofs Clause offset.
   */
  explicit Antecedent(const ClauseOfs cl_ofs) {
    val_ = (cl_ofs << 1) | 1;
  }
  explicit Antecedent(const LiteralID idLit) {
    val_ = (idLit.raw() << 1);
  }
  /**
   * Checks if the antecendent is a clause.  This is done by checking if
   * the least significant bit is a 0b0 or a 0b1.
   *
   * @return True if the antecedent is a clause.  False if it is a literal.
   */
  bool isAClause() const {
    return val_ & 0x01;
  }
  /**
   * Returns the antecedent information formatted as a clause
   *
   * @return Antecedent clause
   */
  ClauseOfs asCl() const {
    assert(isAClause());  // ZH - Verify the antecedent actually is a CLAUSE.
    return val_ >> 1;
  }
  /**
   * Returns the antecedent information formatted as a literal.
   *
   * @return Antecedent literal
   */
  LiteralID asLit() {
    assert(!isAClause());  // ZH - Verify the antecedent actually is a LITERAL.
    LiteralID idLit;
    idLit.copyRaw(val_ >> 1);
    return idLit;
  }
  /**
   * MT - A NON-Antecedent will only be A NOT_A_CLAUSE Clause Id
   *
   * Checks whether this is an actual antecedent.
   *
   * @return true if the this is not a true antecendent (i.e., the value equals 1).
   */
  bool isAnt() const {
    return val_ != 1;  // i.e. NOT a NOT_A_CLAUSE;
  }
};

/**
 * Simple variable structure used for storing variable related information
 * including its antecedent and the decision level in the stack.
 */
struct Variable {
  Antecedent ante;
  long decision_level = INVALID_DL;
};


class ClauseHeader {
  unsigned creation_time_;  // number of conflicts seen at creation time
  unsigned score_;
  unsigned length_;

 public:
  void increaseScore() {
    score_++;
  }
  void decayScore() {
    score_ >>= 1;
  }
  unsigned score() const {
    return score_;
  }

  unsigned creation_time() {
    return creation_time_;
  }
  unsigned length() { return length_; }
  void set_length(unsigned length) { length_ = length; }

  void set_creation_time(unsigned time) {
    creation_time_ = time;
  }
  static unsigned overheadInLits() { return sizeof(ClauseHeader)/sizeof(LiteralID); }
};

#endif /* STRUCTURES_H_ */
