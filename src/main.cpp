/**
 * main.cpp
 *
 * Purpose: Runs the uniform sampler.
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

#include "sysexits.h"

#include <vector>
#include <iostream>
#include "solver.h"

/**
 * Prints a description of the input parameters to the console.
 */
int printInputArgumentDescription();
/**
 * Main sharpSAT function.
 *
 * @param argc Number of input arguments.
 * @param argv Input argument list
 * @return 0 showing successful execution.  Any other value reports that the tool was not executed.
 */
int main(int argc, char *argv[]) {
  if (argc == 1 || (argc == 2 && strcmp(argv[0], "-h") == 0))
    printInputArgumentDescription();

  Solver theSolver(argc, argv);

  if (theSolver.config().perform_random_sampling_)
    theSolver.sample_models();
  else
    theSolver.solve();

  return EXIT_SUCCESS;
}
/**
 * Input Argument Descriptor
 *
 * When no input arguments are specified, this function provides a description
 * of the possible input arguments and settings.
 *
 * Exits the entire program after executing.
 */
int printInputArgumentDescription() {
  std::cout << "Usage: spur [settings]\n"
            << "\t -cnf [cnf_file] Path to the CNF file\n"
            << "\t -s [s] \t Number of models \"s\" to sample uniformly at random\n"
            << "\t -tp    \t Forces two-pass sampling. (Only applicable when s=1)\n"
            << "\t -seed  \t If set to -1 (default), random seed is current time.\n"
            << "\t        \t Otherwise, if non-negative, fixes the random seed.\n"
            << "\t -out [out_file] Path to write the specified samples\n"
            << "\t -no-sample-write Disable writing the final samples to a file.\n"
            << "\t -count-only\t Perform only model counting.  Disable sampling.\n"
            << "\t -q     \t Quiet mode\n"
            << "\t -v     \t Verbose and trace mode\n"
            << "\t -d     \t Debug mode\n"
            << "\t -t [s] \t Set time bound to s seconds\n"
            << "\t -cs [n]\t Set max cache size to n MB\n"
            << std::flush;
  exit(EX_USAGE);
}
