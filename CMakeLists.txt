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
add_subdirectory(ext)

# =========================================
# Real code
# =========================================
add_subdirectory(code)

# =========================================
# Testing
# =========================================
option(BUILD_PLOTTING "Build with plotting examples" ON)
add_subdirectory(test)