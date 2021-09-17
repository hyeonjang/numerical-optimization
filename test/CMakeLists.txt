enable_testing()

# test set
add_executable(assign1_test assign1_test.cpp)
add_executable(assign2_test assign2_test.cpp)

set(TEST_LIST assign1_test assign2_test)
include_directories("../code/")

include(GoogleTest)

foreach(test ${TEST_LIST})
    target_link_libraries(${test} numerical_optimization)
    target_link_libraries(${test} gtest_main)
    gtest_discover_tests(${test})
endforeach()

# benchmark set
add_executable(assign1_benchmark assign1_benchmark.cpp)
add_executable(assign2_benchmark assign2_benchmark.cpp)

set(BENCHMARK_LIST assign1_benchmark assign2_benchmark)
foreach(benchmark ${BENCHMARK_LIST})
    target_link_libraries(${benchmark} numerical_optimization)
    target_link_libraries(${benchmark} benchmark::benchmark)
endforeach()


