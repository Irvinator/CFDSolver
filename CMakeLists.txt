cmake_minimum_required(VERSION 3.20)

# ============================================================================
# vcpkg
# IMPORTANT: Must be set BEFORE project()
# ============================================================================

if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE
        "${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake"
        CACHE STRING "vcpkg toolchain file"
    )
endif()


# ============================================================================
# Project
# ============================================================================

project(
    CFDSolver
    VERSION 0.1.0
    LANGUAGES CXX
)


# ============================================================================
# Options
# ============================================================================

option(BUILD_UI "Build UI application" ON)


# ============================================================================
# C++ settings
# ============================================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# ============================================================================
# Compiler options
# ============================================================================

if(MSVC)

    add_compile_options(
        /W4
        /EHsc
    )

else()

    add_compile_options(
        -Wall
        -Wextra
        -O2
    )

endif()


# ============================================================================
# Core library
# ============================================================================

add_library(CFDCore INTERFACE)

target_include_directories(
    CFDCore
    INTERFACE
        ${CMAKE_SOURCE_DIR}/CFDSolver/src
)


# ============================================================================
# Testing
# ============================================================================

enable_testing()


# ============================================================================
# Vector test
# ============================================================================

add_executable(
    test_vector
    CFDSolver/tests/test_vector.cpp
)

target_link_libraries(
    test_vector
    PRIVATE
        CFDCore
)

add_test(
    NAME VectorTest
    COMMAND test_vector
)


# ============================================================================
# Mesh test
# ============================================================================

add_executable(
    test_mesh2D
    CFDSolver/src/mesh/mesh2D.cpp
    CFDSolver/tests/test_mesh2D.cpp
)

target_link_libraries(
    test_mesh2D
    PRIVATE
        CFDCore
)

add_test(
    NAME Mesh2DTest
    COMMAND test_mesh2D
)


# ============================================================================
# Heat diffusion test
# ============================================================================

add_executable(
    test_heat_diffusion2D
    CFDSolver/tests/test_heat_diffusion2D.cpp
    CFDSolver/src/solvers/HeatSolver2D.cpp
    CFDSolver/src/mesh/mesh2D.cpp
)

target_link_libraries(
    test_heat_diffusion2D
    PRIVATE
        CFDCore
)

add_test(
    NAME HeatDiffusion2DTest
    COMMAND test_heat_diffusion2D
)


# ============================================================================
# Heat diffusion 2D example
# ============================================================================

add_executable(
    heat_diffusion_2D
    CFDSolver/examples/heat_diffusion_2D.cpp
    CFDSolver/src/solvers/HeatSolver2D.cpp
    CFDSolver/src/mesh/mesh2D.cpp
)

target_link_libraries(
    heat_diffusion_2D
    PRIVATE
        CFDCore
)


# ============================================================================
# Heat diffusion 1D example
# ============================================================================

add_executable(
    heat_diffusion_1D
    CFDSolver/examples/heat_diffusion_1D.cpp
)

target_link_libraries(
    heat_diffusion_1D
    PRIVATE
        CFDCore
)


# ============================================================================
# Fields test
# ============================================================================

add_executable(
    test_Fields
    CFDSolver/tests/test_fields.cpp
    CFDSolver/src/mesh/mesh2D.cpp
)

target_link_libraries(
    test_Fields
    PRIVATE
        CFDCore
)

add_test(
    NAME FieldsTest
    COMMAND test_Fields
)


# ============================================================================
# Boundary condition test
# ============================================================================

add_executable(
    tests_BC
    CFDSolver/tests/test_BC.cpp
    CFDSolver/src/BCs/BC.cpp
)

target_include_directories(
    tests_BC
    PRIVATE
        ${CMAKE_SOURCE_DIR}/CFDSolver/src
)


# ============================================================================
# BiCGSTAB test
# ============================================================================

add_executable(
    test_BiCGSTAB
    CFDSolver/tests/test_BiCGSTAB.cpp
    CFDSolver/src/linearAlgebra/Vector.hpp
    CFDSolver/src/linearAlgebra/SparseMatrix.hpp
    CFDSolver/src/linearAlgebra/BiCGSTAB.h
)

target_link_libraries(
    test_BiCGSTAB
    PRIVATE
        CFDCore
)

add_test(
    NAME BiCGSTABTest
    COMMAND test_BiCGSTAB
)


# ============================================================================
# Original SIMPLE test
# ============================================================================

add_executable(
    test_SIMPLE
    CFDSolver/tests/test_SIMPLE.cpp
    CFDSolver/src/solvers/SIMPLE.cpp
    CFDSolver/src/mesh/mesh2D.cpp
    CFDSolver/src/BCs/BC.cpp
)

target_link_libraries(
    test_SIMPLE
    PRIVATE
        CFDCore
)

add_test(
    NAME SIMPLETest
    COMMAND test_SIMPLE
)


# ============================================================================
# Staggered SIMPLE test
# ============================================================================

add_executable(
    test_staggered_SIMPLE
    CFDSolver/tests/test_staggered_SIMPLE.cpp
    CFDSolver/src/solvers/StaggeredSIMPLE.cpp
    CFDSolver/src/mesh/mesh2D.cpp
    CFDSolver/src/BCs/BC.cpp
)

target_link_libraries(
    test_staggered_SIMPLE
    PRIVATE
        CFDCore
)

add_test(
    NAME StaggeredSIMPLETest
    COMMAND test_staggered_SIMPLE
)


# ============================================================================
# CFD GUI APPLICATION
# ============================================================================

if(BUILD_UI)

    message(STATUS "Configuring CFD GUI application...")


    # ------------------------------------------------------------------------
    # UI dependencies
    # ------------------------------------------------------------------------

    find_package(nfd CONFIG REQUIRED)
    find_package(glfw3 CONFIG REQUIRED)
    find_package(glad CONFIG REQUIRED)
    find_package(glm CONFIG REQUIRED)
    find_package(imgui CONFIG REQUIRED)


    # ------------------------------------------------------------------------
    # GUI executable
    # ------------------------------------------------------------------------

    add_executable(
        cfd_app

        # Application
        CFDSolver/examples/cfd_app.cpp

        # Renderer
        CFDSolver/src/renderer/MeshEditor2D.cpp

        # Mesh
        CFDSolver/src/mesh/mesh2D.cpp

        # Heat solver
        CFDSolver/src/solvers/HeatSolver2D.cpp

        # Navier-Stokes solver
        CFDSolver/src/solvers/StaggeredSIMPLE.cpp

        # Boundary conditions
        CFDSolver/src/BCs/BC.cpp

        # File I/O
        CFDSolver/src/IO/MeshReader.cpp
    )


    # ------------------------------------------------------------------------
    # Include directories
    # ------------------------------------------------------------------------

    target_include_directories(
        cfd_app
        PRIVATE
            ${CMAKE_SOURCE_DIR}/CFDSolver/src
            ${CMAKE_SOURCE_DIR}/CFDSolver/examples
    )


    # ------------------------------------------------------------------------
    # Libraries
    # ------------------------------------------------------------------------

    target_link_libraries(
        cfd_app
        PRIVATE

            # CFD core
            CFDCore

            # Native File Dialog
            nfd::nfd

            # GLFW
            glfw

            # OpenGL loader
            glad::glad

            # GLM
            glm::glm

            # ImGui
            imgui::imgui
    )


    # ------------------------------------------------------------------------
    # Make cfd_app the Visual Studio startup project
    # ------------------------------------------------------------------------

    if(MSVC)

        set_property(
            DIRECTORY
            ${CMAKE_CURRENT_SOURCE_DIR}
            PROPERTY
                VS_STARTUP_PROJECT cfd_app
        )

    endif()


endif()


# ============================================================================
# Configuration information
# ============================================================================

message(STATUS "")
message(STATUS "==============================================")
message(STATUS "CFDSolver configuration")
message(STATUS "==============================================")
message(STATUS "C++ standard : C++17")
message(STATUS "BUILD_UI     : ${BUILD_UI}")

if(BUILD_UI)
    message(STATUS "GUI target   : cfd_app")
endif()

message(STATUS "==============================================")
message(STATUS "")