/**
 * statistics.cpp
 *
 * Purpose: Defines methods for the DataAndStatistics() class.
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

#include <fstream>
#include <string>
#include "statistics.h"

//void DataAndStatistics::writeToFile(const std::string & file_name) {
//  std::ofstream out(file_name, std::ios_base::app);
//  unsigned long pos = input_file_.find_last_of("/\\");
//  out << "<tr>" << std::endl;
//  out << "<td>" << input_file_.substr(pos + 1) << "</td>" << std::endl;
//  out << "<td>" << num_original_variables_ << "</td>" << std::endl;
//  out << "<td>" << num_original_clauses_ << "</td>" << std::endl;
//  out << "<td>" << num_decisions_ << "</td>" << std::endl;
//  out << "<td>" << sampler_time_elapsed_ << "</td>" << std::endl;
//
//  std::string s = final_solution_count_.get_str();
//  if (final_solution_count_ == 0)
//    s = "UNSAT";
//  out << "<td>" << s << "</td>" << std::endl;
//  out << "</tr>" << std::endl;
//}

void DataAndStatistics::printShort() {
  if (exit_state_ == SolverExitState::TIMEOUT) {
    std::cout << "\n" << " TIMEOUT !" << std::endl;
    return;
  }
  std::cout << "\n\n"
            << "variables (total / active / free)\t" << num_variables_ << "/"
            << num_used_variables_ << "/" << num_variables_ - num_used_variables_
            << "\n"
            << "clauses (removed) \t\t" << num_original_clauses_ << " ("
            << num_original_clauses_ - num_clauses() << ")" << "\n"
            << "decisions \t\t\t\t" << num_decisions_ << "\n"
            << "conflicts \t\t\t\t" << num_conflicts_ << "\n"
            << "conflict clauses (all/bin/unit) \t"
            << num_conflict_clauses()
            << "/" << num_binary_conflict_clauses_ << "/" << num_unit_clauses_
            << "\n"
            << "failed literals found by implicit BCP \t "
            << num_failed_literals_detected_ << "\n";


  std::cout << "implicit BCP miss rate \t" << implicitBCP_miss_rate() * 100 << "%\n"
            << "bytes cache size     \t" << cache_bytes_memory_usage()  << "\t\n";

  std::cout << "bytes cache (overall) \t" << overall_cache_bytes_memory_stored()
            << "\n"
            << "bytes cache (infra / comps) "
            << (cache_infrastructure_bytes_memory_usage_) << "/"
            << sum_bytes_cached_components_  << "\t\n";

  std::cout << "bytes pure comp data (curr)    "
            << sum_bytes_pure_cached_component_data_  << "\n"
            << "bytes pure comp data (overall) "
            << overall_bytes_pure_stored_component_data_ << "\n";

  std::cout << "bytes cache with sysoverh (curr)    "
            << sys_overhead_sum_bytes_cached_components_  << "\n"
            << "bytes cache with sysoverh (overall) "
            << sys_overhead_overall_bytes_components_stored_ << "\n";


  std::cout << "cache (stores / hits) \t\t\t" << num_cached_components_ << "/"
            << num_cache_hits_ << "\n"
            << "cache miss rate \t\t" << cache_miss_rate() * 100 << "%\n"
            << "avg. variable count (stores / hits) \t" << getAvgComponentSize()
            << "/" << getAvgCacheHitSize() << "\n\n"
            << "\n# solutions " << "\n"
            << final_solution_count_.get_str() << "\n"
            << "\n# END\n\n"
            << "time: " << sampler_time_elapsed_ << "s\n\n";
}
