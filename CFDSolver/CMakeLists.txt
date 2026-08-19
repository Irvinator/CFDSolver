cmake_minimum_required(VERSION 3.20)
project(CFDSolver VERSION 0.1.0 LANGUAGES CXX)

# ── C++17 Standard ───────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Compiler flags ───────────────────────────────────────────────
if(MSVC)
    add_compile_options(/W4 /O2 /EHsc)
else()
    add_compile_options(-Wall -Wextra -O2)
endif()

# ── Include directories ──────────────────────────────────────────
include_directories(src)

# ── Core library (header-only) ───────────────────────────────────
add_library(CFDCore INTERFACE)
target_include_directories(CFDCore INTERFACE src)

# ── Tests ────────────────────────────────────────────────────────
enable_testing()
add_executable(test_vector tests/test_vector.cpp)
target_link_libraries(test_vector PRIVATE CFDCore)
add_test(NAME VectorTest COMMAND test_vector)

# ── Example executables (uncomment when ready) ───────────────────
# add_executable(heat_diffusion_1D examples/heat_diffusion_1D.cpp)
# target_link_libraries(heat_diffusion_1D PRIVATE CFDCore)

# ── Output message ───────────────────────────────────────────────
message(STATUS "CFDSolver v${PROJECT_VERSION} configured")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")