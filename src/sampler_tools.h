/**
 * sampler_tools.h
 *
 * Purpose: Contains generic helper functions used by the solver.
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

#ifndef SAMPLER_TOOLS_
#define SAMPLER_TOOLS_

#define ESCAPE_CHAR "\033["
#define PRINT_BOLD "1"
#define ITEM_SEP ";"
#define END_ESCAPE "m"
#define RESET_CONSOLE_COLOR (ESCAPE_CHAR "0" END_ESCAPE)

#include <sysexits.h>
#include <sys/stat.h>

#include <iostream>
#include <string>
#include <sstream>

#include "solver_config.h"
#include "statistics.h"

#define STR_DECIMAL_BASE 10

#if defined(WIN32) || defined(_WIN32)
  #define FILE_PATH_SEPARATOR "\\"
#else
  #define FILE_PATH_SEPARATOR "/"
#endif

enum PrintColor {
//  INVERT_COLORS  =  7,
  COLOR_BLACK    = 30,
  COLOR_RED      = 31,
  COLOR_GREEN    = 32,
  COLOR_YELLOW   = 33,
  COLOR_BLUE     = 34,
  COLOR_MAGENTA  = 35,
  COLOR_CYAN     = 36,
//  COLOR_WHITE    = 37
};

/**
 * Helper function to write to an output stream in different colors.
 *
 * @param out Output stream (e.g., cerr, cout) to which to write.
 * @param msg Message to write.
 * @param color Color of the message.
 * @param bold True if the text should be bold.
 */
//@{
/**
 * @brief String error message
 */
inline void PrintInColor(std::ostream &out, const std::string &msg,
                         const PrintColor color = COLOR_BLACK, bool bold = true)  {
  out << ESCAPE_CHAR;
  if (bold)
    out << PRINT_BOLD ITEM_SEP;
  out << color;
  out << END_ESCAPE << msg << RESET_CONSOLE_COLOR << std::endl;
}
/**
 * @brief String error message with default of printing to std::cout.
 */
inline void PrintInColor(const std::string &msg,
                         PrintColor color = COLOR_BLACK, bool bold = true)  {
  PrintInColor(std::cout, msg, color, bold);
}
/**
 * @brief Stringstream error message
 */
inline void PrintInColor(std::ostream &out, const std::stringstream &msg,
                         const PrintColor color = COLOR_BLACK, bool bold = true) {
  PrintInColor(out, msg.str(), color, bold);
}
/**
 * @brief Stringstream error message with default of printing to std::cout.
 */
inline void PrintInColor(const std::stringstream &msg,
                         PrintColor color = COLOR_BLACK, bool bold = true)  {
  PrintInColor(std::cout, msg.str(), color, bold);
}
//@}
/**
 * Generic handler for printing warning messages.
 *
 * @param msg Warning message
 */
inline void PrintWarning(const std::string &msg) {
  PrintInColor(std::cout, "WARNING: " + msg, COLOR_YELLOW);
}
/**
 * Generic handler for printing error messages.
 *
 * @param msg Error message
 */
//@{
/**
 * @brief String version of the error printer
 */
inline void PrintError(const std::string &msg) {
  PrintInColor(std::cerr, "ERROR: " + msg, COLOR_RED);
}
/**
 * @brief Stringstream version of the error function.
 */
inline void PrintError(const std::stringstream &msg) { PrintError(msg.str()); }
//@}
/**
 * Exits the entire program after printing the associated message
 *
 * @param msg Error message to be printed
 */
inline void ExitWithError(const std::string &msg, const int err_code=EXIT_FAILURE) {
  PrintError(msg);
  exit(err_code);
}
/**
 * Exits the program due to an invalid input argument/
 *
 * @param msg Error message to be printed
 */
inline void ExitInvalidParam(const std::string &msg) {
  ExitWithError(msg, EX_DATAERR);
}
/**
 * Checks if the specified file exists on disk.
 *
 * @param file_path Path to the file of interest.
 * @return True if the file exists.
 */
inline const bool FileExists(const std::string &file_path) {
  struct stat buffer;
  return (stat(file_path.c_str(), &buffer) == 0);
}

#endif // SAMPLER_TOOLS_
