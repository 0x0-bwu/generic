# Generic Library

[![License](https://img.shields.io/badge/License-Boost%201.0-blue.svg)](https://www.boost.org/LICENSE_1_0.txt)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/std/the-standard)

A modern, header-only C++ library providing comprehensive tools for geometry processing, spatial data structures, concurrency patterns, and general-purpose utilities.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
  - [Quick Start](#quick-start)
- [Components](#components)
  - [Geometry](#geometry)
  - [Tree Structures](#tree-structures)
  - [Concurrency](#concurrency)
  - [Mathematics](#mathematics)
  - [Utilities](#utilities)
- [Building and Testing](#building-and-testing)
- [Documentation](#documentation)
- [Dependencies](#dependencies)
- [Namespace Reference](#namespace-reference)
- [License](#license)

## Overview

Generic is a modern C++17 header-only library designed for high-performance computational geometry, spatial indexing, and parallel computing. With its templated design, the library is easy to integrate and customize for various applications.

### Key Highlights

- **Header-Only**: Easy integration, no compilation required
- **Modern C++17**: Leverages latest C++ features for performance and expressiveness
- **Template-Based**: Flexible and type-safe design
- **Boost Compatible**: Seamless integration with Boost.Polygon and Boost.Geometry
- **Production-Ready**: Comprehensive test suite included

## Features

### 🔷 Geometry Processing
- **2D/3D Geometry Models**: Points, segments, lines, polygons, boxes, spheres, and more
- **Boolean Operations**: Union, intersection, difference operations on geometries
- **Triangulation**: Delaunay triangulation construction and refinement
- **Polygon Operations**: Merging, offsetting, and connectivity analysis
- **Occupancy Grid Mapping**: Spatial occupancy representation
- **Boost Integration**: Direct compatibility with Boost.Polygon and Boost.Geometry algorithms

### 🌳 Spatial Data Structures
- **KdTree**: k-dimensional tree for efficient spatial queries
- **QuadTree**: Hierarchical 2D spatial partitioning
- **RTree**: Rectangle tree for bounding box queries
- **BVH Tree**: Bounding Volume Hierarchy for collision detection and ray tracing

### ⚡ Concurrency Patterns
- **ThreadPool**: Flexible thread pool for task execution
- **TaskFlow**: Directed acyclic graph (DAG) based task scheduling
- **MapReduce**: Parallel map-reduce pattern implementation
- **Thread-Safe Containers**: Lock-free and thread-safe data structures

### 🔢 Mathematics
- **Linear Algebra**: Vector and matrix operations (Eigen3 integration)
- **Interpolation**: Various interpolation methods
- **Polynomial Fitting**: Curve fitting utilities
- **Fast Math**: Optimized mathematical operations

### 🛠️ Utilities
- **Logging**: Flexible logging system
- **Program Options**: Command-line argument parsing
- **File System**: Cross-platform filesystem operations
- **String Formatting**: String manipulation and formatting
- **Color Management**: Color definitions and conversions
- **QR Code**: QR code generation
- **Graph Theory**: Graph models and algorithms

## Getting Started

### Prerequisites

- **Compiler**: C++17 compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)
- **Required**:
  - [Eigen3](http://eigen.tuxfamily.org/) - Linear algebra library
  - [Boost](https://www.boost.org/) (header-only parts) - Some components require Boost
- **Optional** (for specific features):
  - `libboost_serialization` - For serialization support
  - `libpng` - For PNG image I/O support

### Installation

Since Generic is a header-only library, installation is straightforward:

1. Clone the repository:
```bash
git clone https://github.com/0x0-bwu/generic.git
```

2. Add the library to your include path:
```bash
# Add to your CMakeLists.txt
include_directories(/path/to/generic)

# Or when compiling
g++ -std=c++17 -I/path/to/generic your_code.cpp
```

3. Include the headers you need:
```cpp
#include "geometry/Point.hpp"
#include "tree/KdTree.hpp"
#include "thread/ThreadPool.hpp"
// ... and more
```

### Quick Start

Here's a simple example using the geometry and tree components:

```cpp
#include "geometry/Point.hpp"
#include "tree/KdTree.hpp"
#include <vector>

int main() {
    using namespace generic;
    
    // Create some 2D points
    std::vector<geometry::Point2D<double>> points = {
        {0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}, {3.0, 3.0}
    };
    
    // Build a KdTree for efficient spatial queries
    tree::KdTree<geometry::Point2D<double>> kdtree(points);
    
    // Find nearest neighbor
    auto nearest = kdtree.findNearest({1.5, 1.5});
    
    return 0;
}
```

## Components

### Geometry

The geometry module provides comprehensive 2D and 3D geometric primitives and algorithms:

- **Primitives**: Point, Segment, Line, Polygon, PolygonWithHoles, Box, Triangle, Sphere
- **Operations**: Boolean operations, polygon merging, transformations
- **Algorithms**: Triangulation, tetrahedralization, connectivity extraction
- **Integration**: Boost.Polygon and Boost.Geometry adapters

**Example Use Cases**:
- CAD/CAM applications
- Computational geometry
- Geographic Information Systems (GIS)
- Robot path planning

### Tree Structures

Spatial indexing data structures for efficient geometric queries:

- **KdTree**: Optimal for nearest neighbor searches and range queries
- **QuadTree**: 2D spatial partitioning for collision detection
- **RTree**: Bounding rectangle indexing for spatial databases
- **BVH**: Acceleration structure for ray tracing and collision detection

**Example Use Cases**:
- Collision detection
- Ray tracing
- Spatial database indexing
- Point cloud processing

### Concurrency

Modern concurrency patterns for parallel computing:

- **ThreadPool**: Distribute tasks across worker threads
- **TaskFlow**: Build complex task dependency graphs
- **MapReduce**: Parallel data processing pattern
- **Lock-Free Structures**: High-performance concurrent containers

**Example Use Cases**:
- Parallel mesh processing
- Batch geometric computations
- Simulation systems
- Data processing pipelines

### Mathematics

Mathematical utilities and linear algebra support:

- Linear algebra operations (via Eigen3)
- Interpolation methods (linear, cubic, spline)
- Polynomial fitting and evaluation
- Fast math operations and lookup tables
- Filter implementations

### Utilities

General-purpose utility components:

- **Logging**: Configurable logging system
- **Program Options**: Command-line parsing
- **File System**: Cross-platform file operations
- **String Utilities**: Formatting and parsing
- **Color**: Color space management
- **Hash**: Custom hash functions
- **Circuit**: Circuit simulation (MNA, MOR)
- **Image**: Image processing and QR code generation

## Building and Testing

The library includes a comprehensive test suite using Boost Unit Test Framework.

### Building Tests

```bash
# Set environment variables
export EIGEN_PATH=/path/to/eigen3
export BOOST_PATH=/path/to/boost

# Build tests
cd generic/test
mkdir build && cd build
cmake ..
make

# Run tests
./generic_test
```

### Running Specific Tests

```bash
# Run specific test suites
./generic_test --run_test=geometry_tests
./generic_test --run_test=tree_tests
./generic_test --run_test=thread_tests
```

## Documentation

Generate detailed API documentation using Doxygen:

```bash
doxygen Doxyfile
```

Documentation will be generated in the `docs/` directory. Open `docs/html/index.html` in a browser to view.

## Dependencies

### Required Dependencies

- **Eigen3**: Matrix and vector operations
- **Boost** (header-only): Core utilities and algorithms

### Optional Dependencies

Controlled via macros in `common/Macros.hpp`:

| Macro | Library | Purpose |
|-------|---------|---------|
| `GENERIC_BOOST_SERIALIZATION_SUPPORT` | `libboost_serialization` | Object serialization |
| `GENERIC_BOOST_GIL_IO_PNG_SUPPORT` | `libpng` | PNG image I/O |

### Installing Dependencies

**Ubuntu/Debian**:
```bash
sudo apt-get install libeigen3-dev libboost-all-dev libpng-dev
```

**macOS** (via Homebrew):
```bash
brew install eigen boost libpng
```

**Windows** (via vcpkg):
```bash
vcpkg install eigen3 boost-polygon boost-geometry libpng
```

## Namespace Reference

The library is organized under the `generic` namespace:

```
generic/
├── color                    # Color definitions and conversions
├── common                   # Common definitions and utilities
├── fs                       # Filesystem operations
├── format                   # String formatting
├── geometry                 # Geometric models and algorithms
│   ├── boolean             # Boolean operations
│   └── tri                 # Triangulation
├── log                      # Logging system
├── math                     # Mathematical functions
│   └── la                  # Linear algebra
├── parser                   # Parsing utilities
├── program_options          # Command-line option parsing
├── str                      # String utilities
├── thread                   # Concurrency patterns
│   ├── mapreduce           # MapReduce pattern
│   └── taskflow            # TaskFlow pattern
├── tools                    # General tools and utilities
├── topology                 # Graph theory and topology
├── tree                     # Spatial tree structures
└── unit                     # Unit definitions and conversions
```

### Component Overview

| Component | Description |
|-----------|-------------|
| **boolean** | Boolean expression evaluation and operations |
| **circuit** | Circuit simulation (MNA, MOR, Simulator) |
| **common** | Core macros, traits, and system utilities |
| **geometry** | 2D/3D geometry primitives and algorithms |
| **graph** | Graph models and algorithms |
| **image** | Image processing and QR code generation |
| **math** | Mathematical functions and linear algebra |
| **thread** | Thread pools and concurrency patterns |
| **tools** | Utilities (logging, parsing, file I/O) |
| **topology** | Graph topology and index graphs |
| **tree** | Spatial indexing structures |
| **utils** | Helper utilities (indexing, versioning) |

## License

This library is distributed under the [Boost Software License 1.0](LICENSE).

```
Boost Software License - Version 1.0 - August 17th, 2003

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following conditions...
```

See the [LICENSE](LICENSE) file for the full license text.

---

**Contributions**: Contributions are welcome! Please feel free to submit issues and pull requests.

**Author**: [0x0-bwu](https://github.com/0x0-bwu)

**Repository**: [https://github.com/0x0-bwu/generic](https://github.com/0x0-bwu/generic)
