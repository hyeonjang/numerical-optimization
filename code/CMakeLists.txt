cmake_minimum_required(VERSION 3.10)
project(numerical_optimization)

# =====================================================
# create main
# =====================================================
# build option
option(BUILD_WITH_PLOTTING "build with plotting vectors" ON)

# build lib
add_library(numerical_optimization STATIC 
    fwd.h
    method.h
    univariate.h
    multivariate.h

    method.cpp
    univariate.cpp

    multivariate.cpp
    multi/nelder_mead.hpp
    multi/powells.hpp
    multi/cauchys.hpp
    multi/newtons.hpp
    multi/quasi_newtons.hpp
    )
    
# eigen dependencies
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# linking    
target_link_libraries(numerical_optimization Eigen3::Eigen)