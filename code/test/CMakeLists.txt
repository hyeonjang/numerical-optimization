
enable_testing()

add_executable(
    test
    test_bisection.cpp
)

target_link_libraries(
    test
    gtest_main
)

include(GoogleTest)
gtest_discover_tests(test)