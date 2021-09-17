cmake_minimum_required(VERSION 3.10)
project(numerical_optimization)

# create main
add_library(numerical_optimization STATIC method.cpp)

# add boost lib
find_package(Boost REQUIRED)

if (Boost_FOUND)
    target_link_libraries(numerical_optimization ${Boost_LIBRARIES})
endif()

