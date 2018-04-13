/**
 * primitive_type.h
 *
 * Purpose: Base class and type definitions used throughout this project.
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

#include <vector>

#ifndef PRIMITIVE_TYPES_H_
#define PRIMITIVE_TYPES_H_

#define varsSENTINEL  0
#define clsSENTINEL   NOT_A_CLAUSE

/**
 * Identification numbers for variables.
 */
typedef uint64_t VariableIndex;
/**
 * Identification number for clauses.
 */
typedef VariableIndex ClauseIndex;
/// Unsigned Int
typedef VariableIndex ClauseOfs;
typedef VariableIndex ComponentVarAndCls;
/// Unsigned Int
typedef unsigned CacheEntryID;

/// May be negative due to Thurley's implementation.
typedef int64_t DecisionLevel;

static const ClauseIndex NOT_A_CLAUSE(0);
#define SENTINEL_CL NOT_A_CLAUSE

#define BITS_PER_BYTE 8
#define FIRST_VAR 1

/**
 Enumerated class the stores the result of a solver execution.  It is
 returned by the Solver::countSAT method.  It is also stored as
 part of the run statistics in @see DataAndStatistics#
 */
enum class SolverExitState {
  NO_STATE, SUCCESS, TIMEOUT//, ABORTED
};
/**
 * Associated with backtracking and whether a conflict is resolved or backtracking is required.
 *
 * EXIT means the program is completed and should return success.
 */
enum class SolverNextAction {
  EXIT, RESOLVED, PROCESS_COMPONENT, BACKTRACK
};


#ifdef DEBUG
  #define toDEBUGOUT(X) std::cout << X;
#else
  #define toDEBUGOUT(X)
#endif

enum AssignmentEncoding{
  ASSN_F = 0x0, ASSN_T = 0x1, ASSN_U = 0x3,
};
typedef std::vector<AssignmentEncoding> PartialAssignment;
/**
 * Defines the number of samples that can be requested.
 */
typedef uint32_t SampleSize;


/**
 * Type used for encoding literals in the top tree analyzer.
 */
typedef int32_t TopTreeLiteral;
/**
 * Index used to represent top tree nodes.
 */
typedef uint64_t TreeNodeIndex;
/**
 * Enumerated type to represent each of the top tree node types.
 */
enum class TopTreeNodeType {
  MAX_DEPTH,
  CYLINDER,
  COMPONENT_SPLIT,
  NUM_TREE_NODE_TYPES // Always be last
};

#endif /* PRIMITIVE_TYPES_H_ */
