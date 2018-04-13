/**
 * base_packed_component.h
 *
 * Purpose: Defines the classes BitStuffer() and BasePackedComponent().
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

#ifndef BASE_PACKED_COMPONENT_H_
#define BASE_PACKED_COMPONENT_H_

#include <assert.h>
#include <gmpxx.h>
#include <iostream>

//using namespace std;

/**
 * Manages the data in a byte array.  It shift the bits into the array
 * of objects of type T.
 *
 * @tparam T Generally the "unsigned" type.
 */
template <class T> class BitStuffer {
 public:
  explicit BitStuffer(T *data) : data_start_(data), p(data) {
    *p = 0;
  }

  /**
   * Packs the specified signal "val" into the memory.
   *
   * @param val  Binary value to pack into memory.
   * @param num_bits_val The expected bit length of val.
   */
  void stuff(const unsigned val, const unsigned num_bits_val) {
    assert(num_bits_val > 0);  // Verify a number of bits is specified
    assert((val >> num_bits_val) == 0);  // Verify no hanging over bits in the passed in value

    // Clear the memory at the specified location to make sure it is empty
    if (end_of_bits_ == 0)
      *p = 0;

    // Ensure there are no 1 bits after the end_of_bits_ bit location.
    assert((*p >> end_of_bits_) == 0);
    // Insert val into the memory location pointed to by p.
    *p |= val << end_of_bits_;
    // Shift the bit location by the number of bits in the vlue
    end_of_bits_ += num_bits_val;
    if (end_of_bits_ > _bits_per_block) {
      //assert(*p);
      // Roll over the bit count
      end_of_bits_ -= _bits_per_block;
      *(++p) = val >> (num_bits_val - end_of_bits_);
      assert(end_of_bits_ != 0 || (*p == 0));
    } else if (end_of_bits_ == _bits_per_block) {
      end_of_bits_ -= _bits_per_block;
      p++;
    }
  }
   /**
    * Verify the specified size matches the actual contents of the BitStuff
    * object.
    *
    * @param size number of bytes between the data_start_ pointer and the
    * end of the BitStuff managed area.
    */
  void assert_size(unsigned size) {
    if (end_of_bits_ == 0)
       p--;
    assert(p - data_start_ == size - 1);
  }

 private:
  /**
   * Pointer to the start of a type <T> array used to store the cache data
   */
  T *data_start_ = nullptr;
  /**
   * Pointer to the end of the data
   */
  T *p = nullptr;
  /**
   * In the current block the bit position just after the last bit written
   */
  unsigned end_of_bits_ = 0;
  /**
   * Equals 32 (= 4 * 8 bits)
   *
   * sizeof(unsigned) = 4
   * Three left bit shifts equals 8
   */
  static const unsigned _bits_per_block = (sizeof(T) << 3);
};


/**
 * Stores the specifications for package cache components including
 * the number of bits in a variable, clause, and other information.
 */
class BasePackedComponent {
 public:
  /**
   * Accessor to get n := \ceil{lg |var(F)|}, i.e., the number of bits
   * that represent a VARIABLE.
   *
   * @return Number of bits to store a VARIABLE
   */
  static unsigned bits_per_variable() {
    return _bits_per_variable;
  }
  static unsigned variable_mask() {
      return _variable_mask;
  }
  /**
   * Accessor to get m := \ceil{lg |cl(F)|}, i.e., the number of bits
   * that represent a CLAUSE.
   *
   * @return Number of bits to store a CLAUSE.
   */
  static unsigned bits_per_clause() {
    return _bits_per_clause;
  }

  static unsigned bits_per_block() {
    return _bits_per_block;
  }

  static unsigned bits_of_data_size() {
    return _bits_of_data_size;
  }

  /**
   * Configures that variables that define the packing size
   * variables.
   *
   * @param maxVarId Maximum ID number for a variable
   * @param maxClId Maximum ID number for a clause
   */
  static void adjustPackSize(unsigned int maxVarId, unsigned int maxClId);

  BasePackedComponent() = default;

  /**
   * Creates a packed component and records with it the creation time.
   * @param creation_time
   */
  explicit BasePackedComponent(unsigned creation_time) : creation_time_(creation_time) {}

  ~BasePackedComponent() {
//    if (data_)  // Deleting nullptr has no effect so no need for if
//      delete data_;
    delete data_;
  }
  /**
   * Prints a bit string to the console.  The MSB is printed first and
   * the LSB last.
   *
   * @param v Value to be printed by bit to the console.
   */
  static void outbit(unsigned v) {
    for (auto i=0; i < 32; i++) {
      std::cout << ((v&2147483648)?"1":"0");  // 2147483648 == 2^31
      v&=2147483648-1;
      v <<= 1;
    }
  }


  static unsigned log2(unsigned v) {
    // taken from
    // http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup
    static const char LogTable256[256] = {
      #define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
      -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
      LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
      LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
    };

    unsigned r;     // r will be lg(v)
    unsigned int t, tt;  // temporaries

    if ((tt = (v >> 16))) {
      r = (t = (tt >> 8)) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
    } else {
      r = (t = (v >> 8)) ? 8 + LogTable256[t] : LogTable256[v];
    }
    return r;
  }

  unsigned creation_time() {
    return creation_time_;
  }

  const mpz_class &model_count() const {
    return model_count_;
  }

  unsigned alloc_of_model_count() const {
    return sizeof(mpz_class) + model_count_.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
  }

  void set_creation_time(unsigned time) {
    creation_time_ = time;
  }

  void set_model_count(const mpz_class &rn, unsigned time) {
    model_count_ = rn;
    length_solution_period_and_flags_ = (time - creation_time_)
                                        | (length_solution_period_and_flags_ & 1);
  }

  unsigned hashkey() const  {
    return hashkey_;
  }

  bool modelCountFound() {
    return (length_solution_period_and_flags_ >> 1);
  }

  /**
   * Determines if a cache entry is deletable.  This entails
   * that the entry is not connected to an active component
   * in the component stack.
   *
   * @return true if the cache entry is disconnected from an
   * active component and can be deleted.
   */
  bool isDeletable() const {
    return length_solution_period_and_flags_ & 1;
  }
  /**
   * Mark the packed cache entry as deletable.
   */
  void set_deletable() {
    length_solution_period_and_flags_ |= 1;
  }

  /**
   * Free the memory associated with the packed cache entry's data.
   */
  void clear() {
    // before deleting the contents of this component,
    // we should make sure that this component is not present in the component stack anymore!
    assert(isDeletable());
//    if (data_)  // If statement unnecessary as deleting nullptr as no effect
//      delete data_;
    delete data_;
    data_ = nullptr;
  }

  static unsigned _debug_static_val;

 protected:
  /**
   * Packed data array to be stored in memory,
   *
   * data_ contains in packed form the variable indices
   * and clause indices of the component ordered
   * structure is
   *
   * var var ... clause clause ...
   *
   * clauses begin at clauses_ofs_
   */
  unsigned* data_ = nullptr;

  unsigned hashkey_ = 0;

  mpz_class model_count_;

  unsigned creation_time_ = 1;


  // this is:  length_solution_period = length_solution_period_and_flags_ >> 1
  // length_solution_period == 0 means unsolved
  // and the first bit is "delete_permitted"
  unsigned length_solution_period_and_flags_ = 0;

  // deletion is permitted only after
  // the copy of this component in the stack
  // does not exist anymore

 protected:
  /**
   * In the sharpSAT paper, bits per clause (m) is defined as:
   *
   * m := \ceil{lg |cl(F)|}
   *
   * It is the minimum number of bits per clause when packing is enabled.
   *
   * The implementation is slightly different in the code.  It is defined
   * as:
   *
   * m := \floor{lg |cl(F)| + 1}
   */
  static unsigned _bits_per_clause;
  /**
   * In the sharpSAT paper, bits per variable (n) is defined as:
   *
   * n := \ceil{lg |var(F)|}
   *
   * It is the minimum number of bits per clause when packing is enabled.
   *
   * The implementation is slightly different in the code.  It is defined
   * as:
   *
   * n := \floor{lg |var(F)| + 1}
   */
  static unsigned _bits_per_variable;  // bitsperentry
  static unsigned _bits_of_data_size;  // number of bits needed to store the data size.
  static unsigned _data_size_mask;
  /**
   * Bit mask in the form of bits of 0b1 of length "n"
   *
   * @see BasePackedComponent#_bits_per_variable.
   */
  static unsigned _variable_mask;
  /**
   * Bit mask in the form of bits of 0b1 of length "m"
   *
   * @see BasePackedComponent#_bits_per_variable.
   */
  static unsigned _clause_mask;
  /**
   * Equals 32 (= 4 * 8 bits)
   *
   * sizeof(unsigned) = 4
   * Three left bit shifts equals 8
   */
  static const unsigned _bits_per_block = (sizeof(unsigned) << 3);
};

#endif /* BASE_PACKED_COMPONENT_H_ */
