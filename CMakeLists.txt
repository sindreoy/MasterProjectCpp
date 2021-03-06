cmake_minimum_required(VERSION 3.9)
project(MasterProjectCpp)

set(CMAKE_CXX_STANDARD 11)

set(programs main testSolution testLogNormalInitialCondition checkSensitivity chartMeanSSE
        bruteForceParamEstimation)
set(programs_dependencies
        cmake-build-debug/crudeB.csv
        results/raw_logNormal3.txt
        Fluid.cpp Fluid.h
        SystemProperties.cpp SystemProperties.h
        Grid.cpp Grid.h
        Kernels.cpp Kernels.h
        PBModel.cpp PBModel.h
        )

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
FOREACH(program ${programs})

        # example source files
        ADD_EXECUTABLE(${program} ${program}.cpp ${programs_dependencies})

        # libraries to link against
        TARGET_LINK_LIBRARIES(${program} ${SUNDIALS_LIBS} ${GSL_LIBS})

ENDFOREACH(program ${programs})
