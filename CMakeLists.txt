cmake_minimum_required(VERSION 3.9)
project(MasterProjectCpp)

set(CMAKE_CXX_STANDARD 11)


############################# SUNDIALS #############################
# Specify path to SUNDIALS header files
SET(SUNDIALS_INC_DIR
        /Users/Sindre/Sundials/instdir/include
        CACHE STRING
        "Location of SUNDIALS header files")

# Add path to SUNDIALS header files
INCLUDE_DIRECTORIES(${SUNDIALS_INC_DIR})

# Set search path for SUNDIALS libraries
SET(SUNDIALS_LIB_DIR /Users/Sindre/Sundials/instdir/lib)

# Find the SUNDIALS solver's library
FIND_LIBRARY(SUNDIALS_SOLVER_LIB
        sundials_cvode ${SUNDIALS_LIB_DIR}
        DOC "CVODE library")

# Find the NVECTOR library
FIND_LIBRARY(SUNDIALS_NVEC_LIB
        sundials_nvecserial ${SUNDIALS_LIB_DIR}
        DOC "NVECTOR library")

# Set an extra link directory if necessary (false if empty)
IF(EXISTS )
    LINK_DIRECTORIES()
ENDIF()

# Set additional libraries
SET(SUNDIALS_EXTRA_LIB  -lm CACHE STRING "Additional libraries")

# List of Sundials libraries shared across all examples
SET(SUNDIALS_LIBS ${SUNDIALS_SOLVER_LIB} ${SUNDIALS_NVEC_LIB} ${SUNDIALS_EXTRA_LIB})
#####################################################################################


############################# GNU Scientific Library #############################
# Specify path to GSL header files
SET(GSL_INC_DIR
        /Users/sindre/GSL/instdir/include
        CACHE STRING
        "Location of GSL header files")
# Add path to GSL header files
INCLUDE_DIRECTORIES(${GSL_INC_DIR})
# Set search path for GSL libraries
SET(GSL_LIB_DIR /Users/sindre/GSL/instdir/lib)
# Find library
FIND_LIBRARY(GSL_LIBRARY
        gsl ${GSL_LIB_DIR}
        DOC "GSL library")
FIND_LIBRARY(GSL_CBLAS
        gslcblas ${GSL_LIB_DIR}
        DOC "GSL BLAS Library")
# List of GSL libraries
SET(GSL_LIBS ${GSL_LIBRARY} ${GSL_CBLAS})
##################################################################################
add_executable(MasterProjectCpp
    cmake-build-debug/crudeB.csv
    main.cpp
        Fluid.cpp Fluid.h
        SystemProperties.cpp SystemProperties.h
        utilities.cpp utilities.h
        Grid.cpp Grid.h
        Kernels.cpp Kernels.h
        Constants.cpp Constants.h)

TARGET_LINK_LIBRARIES(MasterProjectCpp ${SUNDIALS_LIBS} ${GSL_LIBS})