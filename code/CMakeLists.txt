cmake_minimum_required(VERSION 3.10)
project(numerical_optimization)

set(CMAKE_CXX_STANDARD 17) 
set(CMAKE_CXX_STANDARD_REQUIRED ON) 
set(CMAKE_CXX_EXTENSIONS OFF)
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
    univariate.cpp

    multivariate.h
    multivariate.cpp
    termination.h

    multi/nelder_mead.hpp
    multi/powells.hpp
    multi/cauchys.hpp
    multi/newtons.hpp
    multi/quasi_newtons.hpp
    
    cg/linearcg.hpp
    cg/nonlinearcg.hpp

    lsm.h
    lsm/gauss_newtons.hpp
    lsm/LM.hpp

    global/genetic_algorithm.hpp
    )

# eigen dependencies
find_package(Eigen3 3.3 REQUIRED NO_MODULE)


# linking    
target_link_libraries(numerical_optimization Eigen3::Eigen)