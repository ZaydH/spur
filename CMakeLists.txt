cmake_minimum_required (VERSION 2.8)
project (SPUR)

# MakeFile Option Explanation
# -Wall                     Enable all compiler warnings
# -O3                       Maximum compiler optimization
# -DNDEBUG                  Disable assert checking
# -g                        Produce debugging information
# -fno-omit-frame-pointer   Preserve the frame pointer for use with perf
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -Wall")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -O3 -DNDEBUG -Wall")

set(CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_PROFILING} -std=c++11 -O3 -g
    -DNDEBUG -Wall -fno-omit-frame-pointer")

IF(UNIX)
  IF(APPLE)
    include_directories(/opt/local/include)
    include_directories(/opt/homebrew/include)
  ELSE(APPLE)


  ENDIF(APPLE)
ELSE(UNIX)

ENDIF(UNIX)

# GNU Multiprecision library
find_library(GMP_LIB gmp)
find_library(GMPXX_LIB gmpxx)

add_subdirectory(src)

add_executable (spur src/main.cpp)
target_link_libraries (spur libSpur ${GMP_LIB} ${GMPXX_LIB})
