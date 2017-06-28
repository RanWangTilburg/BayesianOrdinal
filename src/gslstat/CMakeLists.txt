CMAKE_MINIMUM_REQUIRED(VERSION 3.6)
project(gslstat C CXX)

SET(CMAKE_C_COMPILER icc)
SET(CMAKE_CXX_COMPILER icpc)
set(CMAKE_Fortran_COMPILER ifort)

if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0   -Wall -std=c++11 -qopenmp")
    set(CMAKE_Fortran_FLAGS "-O3 -qopenmp  ")
ELSEIF (CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3  -Wall -std=c++11 -qopenmp")
    set(CMAKE_Fortran_FLAGS "-O3 -qopenmp  ")
ENDIF ()

include_directories(../../lib/Eigen3)

#############################Require GSL to be installed##################################################
FIND_PACKAGE(GSL REQUIRED)
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
#LINK_DIRECTORIES(${GSL_LIBRARIES})

#############################Require MKL to be installed##################################################
set(CMAKE_MODULE_PATH
        ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
FIND_PACKAGE(MKL REQUIRED)
INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIRS})
#LINK_DIRECTORIES(${MKL_LIBRARIES})
#############################Find TBB#####################################################################
FIND_PACKAGE(TBB REQUIRED)

ADD_LIBRARY(gslstat SHARED
        src/linalg/matrix.h
        src/linalg/vector.h
        src/random/random_v2
        src/optimization/optimizer.h
        )
SET_TARGET_PROPERTIES(gslstat PROPERTIES LINKER_LANGUAGE CXX)
TARGET_LINK_LIBRARIES(gslstat ${GSL_LIBRARIES} ${MKL_LIBRARIES})

############################Trucated Normal###############################################################
ADD_LIBRARY(tnormal
        tnormal/rtnorm.cpp
        tnormal/rtnorm.hpp)
TARGET_LINK_LIBRARIES(tnormal ${GSL_LIBRARIES} ${MKL_LIBRARIES})
############################TEST##########################################################################
ADD_EXECUTABLE(test_eigen_interop test/eigen_interop.cpp)
TARGET_LINK_LIBRARIES(test_eigen_interop gslstat gtest gtest_main)

ADD_EXECUTABLE(test_random test/test_random.cpp)
TARGET_LINK_LIBRARIES(test_random gslstat gtest gtest_main tnormal)

ADD_EXECUTABLE(tbb_example tbb_cilk/tbb.cpp)
TARGET_LINK_LIBRARIES(tbb_example tbb gtest gtest_main)

ADD_EXECUTABLE(cilk_example tbb_cilk/cilk.cpp)
TARGET_LINK_LIBRARIES(cilk_example tbb gtest gtest_main)


ADD_EXECUTABLE(misc_test test/misc_test.cpp)
TARGET_LINK_LIBRARIES(misc_test gtest gtest_main)
########################################Bayesian Ordinal Regression#######################################
ADD_LIBRARY(bOrdinal
        bOrdinal/bOrdinal.cpp
        bOrdinal/bOrdinal.h)
SET_TARGET_PROPERTIES(bOrdinal PROPERTIES LINKER_LANGUAGE CXX)
TARGET_LINK_LIBRARIES(bOrdinal tnormal ${GSL_LIBRARIES} ${MKL_LIBRARIES})

ADD_EXECUTABLE(test_b_ord test/test_b_ord.cpp)
TARGET_LINK_LIBRARIES(test_b_ord bOrdinal gtest gtest_main)