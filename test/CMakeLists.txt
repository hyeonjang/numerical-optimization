enable_testing()

# test set
add_executable(hw1_test hw1/test.cpp)
add_executable(hw2_test hw2/test.cpp)
add_executable(hw3_test hw3/test.cpp)

set(TEST_LIST hw1_test hw2_test hw3_test)
include_directories("../code/")

include(GoogleTest)

foreach(test ${TEST_LIST})
    target_link_libraries(${test} numerical_optimization)
    target_link_libraries(${test} gtest_main)
    gtest_discover_tests(${test})
endforeach()

# benchmark set
add_executable(hw1_benchmark hw1/benchmark.cpp)
add_executable(hw2_benchmark hw2/benchmark.cpp)
add_executable(hw3_benchmark hw3/benchmark.cpp)

set(BENCHMARK_LIST hw1_benchmark hw2_benchmark hw3_benchmark)
foreach(benchmark ${BENCHMARK_LIST})
    target_link_libraries(${benchmark} numerical_optimization)
    target_link_libraries(${benchmark} benchmark::benchmark)
endforeach()

# plot set
add_executable(hw3_plot hw3/plot.cpp)
target_link_libraries(hw3_plot numerical_optimization)