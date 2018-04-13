/**
 * containers.h
 *
 * Purpose: Defines the LiteralIndexVector() class.  This stores items and access them via
 * literal identification numbers.
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

#ifndef CONTAINERS_H_
#define CONTAINERS_H_

#include <vector>

#include "structures.h"

template<class _T>
/***
 * @tparam _T Type of element to store in the vector.
 */
class LiteralIndexedVector: protected std::vector<_T> {
 public:
  /**
   * Creates a vector twice the specified size.
   *
   * @param size Expected number of variables
   */
  explicit LiteralIndexedVector(VariableIndex size = 0) :
      std::vector<_T>(size * 2) {
  }
  /**
   * Create a vector twice the specified size initialized to the specified
   * value.
   *
   * @param size Expected number of variables.
   *
   * @param __value Value to write into all locations in the vector
   */
  LiteralIndexedVector(VariableIndex size, const typename std::vector<_T>::value_type& __value)
      : std::vector<_T>(size * 2, __value) {
  }
  inline _T &operator[](const LiteralID lit) {
    return *(std::vector<_T>::begin() + lit.raw());
  }

  inline const _T &operator[](const LiteralID &lit) const {
    return *(std::vector<_T>::begin() + lit.raw());
  }
  /**
   * Creates an iterator that points to the first element in the vector.
   *
   * @return Iterator to first element
   */
  inline typename std::vector<_T>::iterator begin() {
    return std::vector<_T>::begin() + 2;
  }
  /**
   *
   *
   * @param _size Expected number of variables.
   */
  void resize(VariableIndex _size) {
    std::vector<_T>::resize(_size * 2);
  }
  /**
   * Resizes the vector to twice the specified size.  If a value is specified,
   * it sets all the elements to the specified value.
   *
   * @param _size
   * @param _value
   */
  void resize(VariableIndex _size, const typename std::vector<_T>::value_type& _value) {
    std::vector<_T>::resize(_size * 2, _value);
  }

  void reserve(VariableIndex _size) {
    std::vector<_T>::reserve(_size * 2);
  }
  /**
   * Calculated version of the last literal ID number.  It is based off the size of the literal array
   * and may have unexpected behavior if the numbers are not incrementally increasing.
   *
   * @return Negated version of the largest possible literal
   */
  LiteralID end_lit() {
    return LiteralID(size() / 2, false);
  }

  // Methods reused from the imported class
  using std::vector<_T>::end;
  using std::vector<_T>::size;
  using std::vector<_T>::clear;
  using std::vector<_T>::push_back;
};

#endif /* CONTAINERS_H_ */
