cmake_minimum_required(VERSION 3.20)
project(aerial_autonomy_lab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(AAL_BUILD_TESTS "Build the C++ test suite" ON)
option(AAL_DISABLE_CUDA "Force the CPU backend even if a CUDA compiler exists" OFF)

# CUDA is optional and auto-detected. Every kernel has a CPU fallback selected at
# runtime, and the full test suite must pass with AAL_DISABLE_CUDA=ON, because
# not every machine that reads this repository has a GPU.
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER AND NOT AAL_DISABLE_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)

  # Real target hardware for this project, so build a fatbin rather than pinning
  # one architecture.
  #   89  Ada Lovelace   RTX 4090
  #   90  Hopper         H100, the training box
  #   120 Blackwell      RTX 5090, consumer Blackwell and NOT the same as 100,
  #                      which is datacenter Blackwell
  # sm_120 needs CUDA 12.8 or newer, so drop it on older toolkits rather than
  # failing the configure step.
  set(AAL_CUDA_ARCHS 89 90)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 12.8)
    list(APPEND AAL_CUDA_ARCHS 120)
  else()
    message(STATUS "CUDA ${CMAKE_CUDA_COMPILER_VERSION} predates sm_120, RTX 50 series excluded")
  endif()
  set(CMAKE_CUDA_ARCHITECTURES ${AAL_CUDA_ARCHS})
  message(STATUS "CUDA enabled, architectures: ${CMAKE_CUDA_ARCHITECTURES}")

  # Match the host FP policy. --fmad=false is the CUDA equivalent of
  # -ffp-contract=off, and without it the GPU and CPU paths disagree in the low
  # bits, which makes the CPU fallback useless as a reference.
  add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>)
  add_compile_definitions(AAL_WITH_CUDA=1)
else()
  message(STATUS "CUDA disabled, building the CPU backend only")
endif()

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
