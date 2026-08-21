cmake_minimum_required(VERSION 3.20)
project(CFDSolver VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/W4 /EHsc)
else()
    add_compile_options(-Wall -Wextra -O2)
endif()

include_directories(src)

add_library(CFDCore INTERFACE)
target_include_directories(CFDCore INTERFACE src)

enable_testing()
add_executable(test_vector tests/test_vector.cpp)
target_link_libraries(test_vector PRIVATE CFDCore)
add_test(NAME VectorTest COMMAND test_vector)

add_executable(heat_diffusion_1D examples/heat_diffusion_1D.cpp)
target_link_libraries(heat_diffusion_1D PRIVATE CFDCore)

message(STATUS "CFDSolver v${PROJECT_VERSION} configured")