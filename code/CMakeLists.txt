cmake_minimum_required(VERSION 3.10)
project(numerical_optimization)

# create main
add_library(numerical_optimization STATIC 
    fwd.h
    univariate.h
    multivariate.h

    univariate.cpp
    multivariate.cpp
    )