cmake_minimum_required(VERSION 3.20)
project(aerial_autonomy_lab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(AAL_BUILD_TESTS "Build the C++ test suite" ON)

# Floating point determinism.
#
# Golden-trajectory tests compare numbers across machines and across compilers,
# so the compiler must not be allowed to reassociate arithmetic or contract
# a*b+c into a single fused multiply-add at its discretion. -ffp-contract=off is
# the important one: GCC defaults to `fast` even without -ffast-math, and an FMA
# changes the last bits of every dot product and distance computation. MSVC's
# /fp:precise is the equivalent default, so no flag is needed there.
if(MSVC)
  add_compile_options(/fp:precise /W4)
else()
  add_compile_options(-ffp-contract=off -Wall -Wextra -Wpedantic)
endif()

add_library(aal_core
  src/collision/esdf.cpp
  src/planner/rrt_star.cpp
)
target_include_directories(aal_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

if(AAL_BUILD_TESTS)
  enable_testing()
  find_package(GTest REQUIRED)
  add_executable(aal_tests
    tests/cpp/test_rng.cpp
    tests/cpp/test_esdf.cpp
    tests/cpp/test_rrt_star.cpp
  )
  target_link_libraries(aal_tests PRIVATE aal_core GTest::gtest GTest::gtest_main)
  include(GoogleTest)
  gtest_discover_tests(aal_tests)
endif()

add_executable(bench_planners bench/bench_planners.cpp)
target_link_libraries(bench_planners PRIVATE aal_core)
