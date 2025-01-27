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
set(BUILD_SHARED_LIBS ON)

## ================= ##
## Create Executable ##
## ================= ##
add_executable(${PROJECT} main.cpp)
add_subdirectory(${PROJECT_SOURCE_DIR}/external/hypre/src)

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_definitions(${PROJECT} PRIVATE RELEASE=1)
else()
    target_compile_definitions(${PROJECT} PRIVATE RELEASE=0)
    # target_compile_options(${PROJECT} PRIVATE Wall Wextra)
endif()

target_sources(${PROJECT}
    # no need to add headers here, only sources are required
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src/main/core/logger.cpp
        ${PROJECT_SOURCE_DIR}/src/main/mesh/mesh.cpp
        ${PROJECT_SOURCE_DIR}/src/main/mesh/polynomials.cpp
        ${PROJECT_SOURCE_DIR}/src/main/physics/integrator_dOmega.cpp
        ${PROJECT_SOURCE_DIR}/src/main/physics/integrator_Omega.cpp
        ${PROJECT_SOURCE_DIR}/src/main/physics/integrator_assembly.cpp
)

target_include_directories(${PROJECT} 
    PRIVATE
        # where the project itself will look for internal headers
        ${PROJECT_SOURCE_DIR}/src/main/
        ${PROJECT_SOURCE_DIR}/src/main/core/
        ${PROJECT_SOURCE_DIR}/src/main/mesh/
        ${PROJECT_SOURCE_DIR}/src/main/physics/
    PUBLIC
        # where the project will look for public headers
        ${PROJECT_SOURCE_DIR}/external/eigen/
        ${PROJECT_SOURCE_DIR}/external/petsc/
)

## ================= ##
## Rerout Executable ##
## ================= ##
target_link_libraries(${PROJECT} HYPRE)
set_target_properties(${PROJECT}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
)
