
cmake_minimum_required(VERSION 3.20)
project(CFDSolver VERSION 0.1.0 LANGUAGES CXX)

option(BUILD_UI "Build UI application" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/W4 /EHsc)
else()
    add_compile_options(-Wall -Wextra -O2)
endif()

add_library(CFDCore INTERFACE)
target_include_directories(CFDCore
    INTERFACE ${CMAKE_SOURCE_DIR}/CFDSolver/src)

enable_testing()

add_executable(test_vector
    CFDSolver/tests/test_vector.cpp)
target_link_libraries(test_vector PRIVATE CFDCore)
add_test(NAME VectorTest COMMAND test_vector)

add_executable(test_mesh2D
    CFDSolver/src/mesh/mesh2D.cpp
    CFDSolver/tests/test_mesh2D.cpp)
target_link_libraries(test_mesh2D PRIVATE CFDCore)
add_test(NAME Mesh2DTest COMMAND test_mesh2D)

add_executable(heat_diffusion_1D
    CFDSolver/examples/heat_diffusion_1D.cpp)
target_link_libraries(heat_diffusion_1D PRIVATE CFDCore)

add_executable(test_Fields
    CFDSolver/tests/test_fields.cpp
    CFDSolver/src/mesh/mesh2D.cpp
)

target_link_libraries(test_Fields PRIVATE CFDCore)

add_test(NAME FieldsTest COMMAND test_Fields)

add_executable(tests_BC
    CFDSolver/tests/test_BC.cpp
    CFDSolver/src/BCs/BC.cpp
)

target_include_directories(tests_BC
    PRIVATE
        ${CMAKE_SOURCE_DIR}/CFDSolver/src
)

add_executable(test_BiCGSTAB
    CFDSolver/tests/test_BiCGSTAB.cpp
    CFDSolver/src/linearAlgebra/Vector.hpp
    CFDSolver/src/linearAlgebra/SparseMatrix.hpp
    CFDSolver/src/linearAlgebra/BiCGSTAB.h
)

target_link_libraries(test_BiCGSTAB
    PRIVATE CFDCore
)

add_test(NAME BiCGSTABTest COMMAND test_BiCGSTAB)

if(BUILD_UI)
    find_package(glfw3 CONFIG REQUIRED)
    find_package(glad  CONFIG REQUIRED)
    find_package(glm   CONFIG REQUIRED)
    find_package(imgui CONFIG REQUIRED)

    add_executable(cfd_app
        CFDSolver/examples/cfd_app.cpp)
    target_link_libraries(cfd_app
        PRIVATE
        CFDCore
        glfw
        glad::glad
        glm::glm
        imgui::imgui)
endif()

message(STATUS "CFDSolver configured")