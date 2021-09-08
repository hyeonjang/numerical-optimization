enable_testing()

add_executable(
    assignment1
    assignment1.cpp
)

target_link_libraries(
    assignment1
    gtest_main
)

include_directories("../code/")
target_link_libraries(assignment1 numerical_optimization)


include(GoogleTest)
gtest_discover_tests(assignment1)