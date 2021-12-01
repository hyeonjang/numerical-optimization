cmake_minimum_required(VERSION 3.14)
project(numerical-optimization)

set(CMAKE_CXX_STANDARD 17) 
set(CMAKE_CXX_STANDARD_REQUIRED ON) 
set(CMAKE_CXX_EXTENSIONS OFF)

# =======================================
# Google test
# ========================================
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# =========================================
# Google Benchmark
# =========================================
FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG origin/main
    GIT_SHALLOW TRUE
)
FetchContent_GetProperties(benchmark)
if(NOT benchmark_POPULATED)
    FetchContent_Populate(benchmark)
    set(BENCHMARK_ENABLE_TESTING OFF CACHE INTERNAL "" FORCE)
    add_subdirectory(${benchmark_SOURCE_DIR} ${benchmark_BINARY_DIR})
endif()

# =========================================
# Real code
# =========================================
add_subdirectory(code)

# =========================================
# Testing
# =========================================
option(BUILD_PLOTTING "Build with plotting examples" ON)
if(BUILD_PLOTTING)
# =========================================
# Json parsing lib
# ========================================
  FetchContent_Declare(
    json 
    GIT_REPOSITORY https://github.com/nlohmann/json.git 
    GIT_TAG v3.7.3
  )
  FetchContent_GetProperties(json)
  if(NOT json_POPULATED)
    FetchContent_Populate(json)
    add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
endif(BUILD_PLOTTING)
add_subdirectory(test)