# CFDSolver

**Shape the model. Solve the physics.**

CFDSolver is a C++ computational fluid dynamics and heat-transfer project focused on learning-oriented, modular solver development. It currently includes structured 2D mesh support, linear algebra utilities, scalar/vector field types, and finite-volume solvers for diffusion problems, with work in progress toward incompressible flow methods.[doc4]

## Current Features

- Structured 2D Cartesian mesh generation with cell-center connectivity for finite-volume methods.[doc4]
- Core linear algebra components including vectors, sparse matrices, and iterative solvers.[doc4]
- Scalar and vector field containers for solver variables.[doc4]
- 1D heat diffusion example and 2D heat diffusion solver/test workflow.[doc4]
- CMake-based build system with optional UI build support via GLFW, GLAD, GLM, and ImGui when enabled.[doc4]
- Unit-style test executables under the `tests/` directory.[doc4]

## Repository Structure

```text
CFDSolver/
├── docs/
├── examples/
│   ├── cfd_app.cpp
│   ├── heat_diffusion_1D.cpp
│   └── heat_diffusion_2D.cpp
├── src/
│   ├── BCs/
│   ├── fields/
│   ├── geometry/
│   ├── IO/
│   ├── linearAlgebra/
│   ├── mesh/
│   ├── numerics/
│   ├── renderer/
│   ├── solvers/
│   └── turbulence/
├── tests/
│   ├── test_BC.cpp
│   ├── test_fields.cpp
│   ├── test_heat_diffusion2D.cpp
│   ├── test_mesh2D.cpp
│   └── test_vector.cpp
├── CMakeLists.txt
└── CMakePresets.json
```

This layout follows the usual README role of explaining what the project does, why it is useful, and how to get started.[doc4]

## Build Requirements

- C++17 compiler.[doc4]
- CMake 3.20 or newer.[doc4]
- Ninja or another supported CMake generator.[doc4]
- Optional UI dependencies through vcpkg when `BUILD_UI=ON`.[doc4]

## Recommended Build Setup

This project uses CMake. If you are using vcpkg for optional UI dependencies, pass the vcpkg toolchain file during configuration rather than setting it inside `CMakeLists.txt`.[doc4]

### Configure without UI

```bash
cmake -S . -B out/build/x64-debug -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DBUILD_UI=OFF
```

### Configure with UI and vcpkg

```bash
cmake -S . -B out/build/x64-debug -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake \
  -DVCPKG_TARGET_TRIPLET=x64-windows \
  -DBUILD_UI=ON
```

### Build

```bash
cmake --build out/build/x64-debug
```

## Running Tests

CTest can run configured test executables once the project has been built.[doc4]

```bash
ctest --test-dir out/build/x64-debug --output-on-failure
```

You can also run individual test executables directly, for example:

```bash
out/build/x64-debug/test_mesh2D.exe
out/build/x64-debug/test_heat_diffusion2D.exe
```

## Running Examples

### 1D Heat Diffusion

```bash
out/build/x64-debug/heat_diffusion_1D.exe
```

### 2D Heat Diffusion

```bash
out/build/x64-debug/heat_diffusion_2D.exe
```

The 2D example writes output to a folder such as:

```text
output/heat2D_example
```

## Numerical Focus

The project is currently most mature for:

- structured finite-volume discretization,
- diffusion-type equations,
- sparse linear system assembly,
- iterative solution strategies.

Planned and ongoing directions include:

- improved 2D/3D solver structure,
- pressure-based incompressible flow methods,
- staggered-grid support,
- mesh/geometry tooling,
- visualization and UI integration.

## Development Notes

- `mesh2D` is currently a structured Cartesian mesh abstraction intended for orthogonal finite-volume discretizations.
- `HeatSolver2D` is a diffusion solver built on that mesh and the core sparse linear algebra types.
- More advanced incompressible flow work will likely require additional staggered-grid or face-based data structures.

## Contributing

If you contribute new solver modules or tests, keep the project modular:

- mesh logic in `src/mesh/`
- field types in `src/fields/`
- linear algebra in `src/linearAlgebra/`
- solver logic in `src/solvers/`
- runnable demos in `examples/`
- checks and regression-style tests in `tests/`

README files are intended to help users get started quickly, while deeper documentation can live in the `docs/` directory.[doc4]

## Maintainer

Repository: [Irvinator/CFDSolver](https://github.com/Irvinator/CFDSolver)

- references and numerical background.

