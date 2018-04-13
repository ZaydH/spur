/**
 * base_packed_component.cpp
 *
 * Purpose: Defines the static variables and methods for the BasePackedComponent() class.
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

#include "base_packed_component.h"
#include <math.h>  // Needed since taking log based 2.
#include <iostream>

unsigned BasePackedComponent::_bits_per_clause = 0;
unsigned BasePackedComponent::_bits_per_variable = 0;  // bitsperentry
unsigned BasePackedComponent::_variable_mask = 0;
unsigned BasePackedComponent::_clause_mask = 0;  // bitsperentry
unsigned BasePackedComponent::_debug_static_val = 0;
unsigned BasePackedComponent::_bits_of_data_size = 0;
unsigned BasePackedComponent::_data_size_mask = 0;


void BasePackedComponent::adjustPackSize(unsigned int maxVarId, unsigned int maxClId) {
  _bits_per_variable = log2(maxVarId) + 1;
  _bits_per_clause   = log2(maxClId) + 1;

  _bits_of_data_size = log2(maxVarId + maxClId) + 1;

  _variable_mask = _clause_mask = _data_size_mask = 0;
  // Iteratively bit shift to create the mask variable
  for (unsigned int i = 0; i < _bits_per_variable; i++)
    _variable_mask = (_variable_mask << 1) + 1;
  for (unsigned int i = 0; i < _bits_per_clause; i++)
    _clause_mask = (_clause_mask << 1) + 1;
  for (unsigned int i = 0; i < _bits_of_data_size; i++)
    _data_size_mask = (_data_size_mask << 1) + 1;
}




