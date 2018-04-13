/**
 * component_archetype.h
 *
 * Purpose: Defines the static variables for the ComponentArchetype() class.
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

#include "component_archetype.h"

CA_SearchState *ComponentArchetype::seen_ = nullptr;
unsigned ComponentArchetype::seen_byte_size_ = 0;
