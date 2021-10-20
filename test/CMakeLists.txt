enable_testing()

# =========================================================
# test set
# =========================================================
add_executable(hw1_test hw1/test.cpp)
add_executable(hw2_test hw2/test.cpp)
add_executable(hw3_test_gradient hw3/test/gradient.cpp)
add_executable(hw3_test_termination hw3/test/termination.cpp)
add_executable(hw4_test hw4/test.cpp)

set(TEST_LIST 
        hw1_test 
        hw2_test 
        hw3_test_gradient hw3_test_termination
        hw4_test    
    )
include_directories("../code/")

include(GoogleTest)

foreach(test ${TEST_LIST})
    target_link_libraries(${test} numerical_optimization)
    target_link_libraries(${test} gtest_main)
    gtest_discover_tests(${test})
endforeach()

# =========================================================
# benchmark set
# =========================================================
add_executable(hw1_benchmark hw1/benchmark.cpp)
add_executable(hw2_benchmark hw2/benchmark.cpp)
add_executable(hw3_benchmark hw3/benchmark.cpp)
add_executable(hw4_benchmark hw4/benchmark.cpp)

set(BENCHMARK_LIST hw1_benchmark hw2_benchmark hw3_benchmark hw4_benchmark)
foreach(benchmark ${BENCHMARK_LIST})
    target_link_libraries(${benchmark} numerical_optimization)
    target_link_libraries(${benchmark} benchmark::benchmark)
endforeach()

# =========================================================
# plot set
# =========================================================
if(BUILD_PLOTTING)
add_compile_definitions(BUILD_WITH_PLOTTING=True)
add_executable(hw3_plot hw3/plot.cpp)
add_executable(hw4_plot hw4/plot.cpp)

set(PLOT_LIST hw3_plot hw4_plot)
find_package(Boost REQUIRED filesystem iostreams)
    if(Boost_FOUND)
        foreach(plot ${PLOT_LIST})
            target_link_libraries(${plot} numerical_optimization)
            target_link_libraries(${plot} Boost::filesystem Boost::iostreams)
        endforeach()
    endif(Boost_FOUND)
endif(BUILD_PLOTTING)
    