cmake_minimum_required(VERSION 3.10)
project(numerical_optimization)

# eigen dependencies
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# create main
add_library(numerical_optimization STATIC 
    fwd.h
    method.h
    univariate.h
    multivariate.h

    method.cpp
    univariate.cpp
    multivariate.cpp
    )

# linking    
target_link_libraries(numerical_optimization Eigen3::Eigen)