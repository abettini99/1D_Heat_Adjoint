cmake_minimum_required(VERSION 3.22.1)    # CMake version check
set(CMAKE_CXX_STANDARD 23)                # Enable C++23 standard
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN hidden)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

## ============== ##
## Define Project ##
## ============== ##
set(PROJECT "HeatAdjoint")
project(PROJECT)

## ================= ##
## Create Executable ##
## ================= ##
add_executable(${PROJECT} main.cpp)

target_sources(${PROJECT}
    # no need to add headers here, only sources are required
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src/main/core/logger.cpp
        ${PROJECT_SOURCE_DIR}/src/main/preconditioner/Jacobi.cpp
        ${PROJECT_SOURCE_DIR}/src/main/preconditioner/incompleteCholesky.cpp
        ${PROJECT_SOURCE_DIR}/src/main/solver/CG.cpp
        ${PROJECT_SOURCE_DIR}/src/main/mesh/mesh.cpp
        ${PROJECT_SOURCE_DIR}/src/main/mesh/polynomials.cpp
)

target_include_directories(${PROJECT} 
    PRIVATE
        # where the project itself will look for internal headers
        ${PROJECT_SOURCE_DIR}/src/main/
        ${PROJECT_SOURCE_DIR}/src/main/core/
        ${PROJECT_SOURCE_DIR}/src/main/mesh/
        ${PROJECT_SOURCE_DIR}/src/main/preconditioner/
        ${PROJECT_SOURCE_DIR}/src/main/solver/
    PUBLIC
        # where the project will look for public headers
        ${PROJECT_SOURCE_DIR}/external/eigen/
)


## ================= ##
## Rerout Executable ##
## ================= ##
set_target_properties(${PROJECT}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
)
