enable_testing()

# test set
add_executable(assign1_test assign1_test.cpp)

include_directories("../code/")
target_link_libraries(assign1_test numerical_optimization)

target_link_libraries(assign1_test gtest_main)

include(GoogleTest)
gtest_discover_tests(assign1_test)

# benchmark set
add_executable(assign1_benchmark assign1_benchmark.cpp)

target_link_libraries(assign1_benchmark numerical_optimization)

target_link_libraries(assign1_benchmark benchmark::benchmark)
