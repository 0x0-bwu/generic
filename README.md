# Generic Library

## Introduction
A header-only library for geometry handling and others.<br>

Implement with C++17 features, templated design for easy usage and embedding.<br>

It has multiple 2d/3d geometry models, and model adaption with Boost Polygon/Geometry Library, so the algorithms in these two libraries can be applied directly with those models.<br>

It offers lots of basic geometry test utility as well as high level algorithms like polygon merge, multi-layer geometry connectivity extraction, occupancy grid mapping and Delaunay triangulation construction and refinement, etc.<br>

It offers templated tree structures like KdTree, QuadTree, RTree, and BVH Tree for geometry indexing and other algorithms.<br>

It has three concurrency patterns ThreadPool, TaskFlow and MapReduce to support parallel computing in generic library, these patterns can also be applied to parallel programming for others easily.<br>

It also has lots of other header-only utilities like log and program option for convenient usage.<br>

## Dependency
Some components depend on Eigen3 and Boost library(header-only), need include boost library path when compiling.<br>

Some function implementation depends on third part library binaries with the Macro control defined in common/Macros.hpp:<br>       

if defined GENERIC_BOOST_SERIALIZATION_SUPPORT, need link with libboost_serialization.<br>
if defined GENERIC_BOOST_GIL_IO_PNG_SUPPORT, need link with libpng.<br>

## Namespace Description
+ `generic`                         : top namespace<br> 
    + `color`                       : color related definitions, classes and functions<br>  
    + `common`                      : common definitions, classes, and functions<br>  
    + `fs`                          : filesystem related functions<br>  
    + `format`	                    : string format<br>  
    + `geometry`                    : models, boost polygon/geometry adaption, geometry algorithms<br>   
        + `boolean`     	        : geometry Boolean operation APIs<br>    
	    + `tri`			            : triangulation data, construction, and refinement<br>     
    + `log` 				        : log library<br>  
    + `math`				        : math numbers, classes, and functions<br>  
	    + `la`			            : linear algebra<br>  
    + `parser` 			            : parser functions<br>  
    + `program_options`		        : program option library<br>  
    + `str`				            : string related functions<br>  
    + `thread`			            : concurrency patterns, thread pool<br>  
	    + `mapreduce`		        : pattern of map reduce<br>     
	    + `taskflow`		        : pattern of task flow<br>  
    + `tools`				        : tools and basic utilities<br>  
    + `topology`			        : graph related<br>  
    + `tree`				        : tree structures and algorithms<br>  
    + `unit`				        : unit definition and convertion<br>

## Documents
You can generate documents of this library by run `doxygen Doxyfile` in current folder, the documents will be generated in folder `docs/`.<br>  

## API Reference
<!-- AUTO_DOCS_START -->
### Header Files Overview

#### boolean/
- **Expression.hpp**: boolean expresion
- **Operation.hpp**: boolean related operations

#### circuit/
- **MNA.hpp**: fork MNA implementation from jefftrull
- **MOR.hpp**: model order reduction
- **Simulator.hpp**: static/transient circuit simulator

#### common/
- **Archive.hpp**: Serialization header files
- **Exception.hpp**: Exception handle
- **Macros.hpp**: Macro defines
- **System.hpp**: System related functions
- **Traits.hpp**: Common traits

#### geometry/
- **BooleanOperation.hpp**: Boolean operation for geometries
- **BoostGeometryRegister.hpp**: Adaption of the generic geometry models to boost geometry concept
- **BoostPolygonRegister.hpp**: Adaption of the generic geometry models to boost polygon concept
- **Box.hpp**: Model of axis-aligned bounding box2d and box3d concept
- **Clipper.hpp**: Modified version of clipper library, origin: http://www.angusj.com/delphi/clipper.php
- **Common.hpp**: Common geometry definition
- **Connectivity.hpp**: Connectivity extraction algorithm to build connectivity graph on layer based geometries
- **Curves.hpp**: Model of some curve concepts
- **Geometries.hpp**: Header file that include some of the geometris model
- **GeometryIO.hpp**: I/O functions of geometries
- **GeometryTraits.hpp**: Geometry traits
- **HashFunction.hpp**: Hush functions of geometries
- **Line.hpp**: Model of line concept
- **Mesh2D.hpp**: 2d mesh flow
- **OccupancyGridMap.hpp**: Grid map define and generation algorithm
- **Plane.hpp**: Model of plane concept
- **Point.hpp**: Model of point2d and point3d concept
- **Polygon.hpp**: Model of polygon concept
- **PolygonMerge.hpp**: Utility for polygon merge
- **PolygonWithHoles.hpp**: Model of polygon with holes concept
- **Predicates.hpp**: Robust geometric predicates for floating points, modified from predicates by William C. Lenthe.
- **Rasterization.hpp**: Rasterize geometry outline to pixel indices
- **Segment.hpp**: Model of segment concept
- **Serialization.hpp**: Boost serialization functions for geometry classes
- **Sphere.hpp**: Model of circle and sphere concept
- **Tetrahedralization.hpp**: Model of tetrahedralization concept
- **TetrahedralizationIO.hpp**: I/O functions of geometries
- **Topology.hpp**: Geometry topology relationship
- **Transform.hpp**: Model of transform vector, matrix and quaternion concept
- **Trapezoid.hpp**: Model of trapezoid concept
- **Triangle.hpp**: Model of triangle2d and triangle3d concept
- **TriangleEvaluator.hpp**: Utility class for triangulation quality evaluation
- **Triangulation.hpp**: Model of triangulation concept
- **TriangulationOperator.hpp**: Triangulation operator that manipulating the triangulation data
- **TriangulationRefinement.hpp**: Triangulation refinement algorithms
- **Triangulator.hpp**: Triangulator algrithom for triangulation
- **Utility.hpp**: Geometry related utilities
- **Vector.hpp**: Model of vector2d, vector3d and hyper vector concept

#### image/
- **Exception.hpp**: QRcode, modified from https://github.com/soleilpqd/QRMatrix.CPP

#### la/
- **Common.hpp**: common define

#### math/
- **FastMath.hpp**: fast math operations from Paul Mineiro's FastFloat
- **Filter.hpp**: filter related algorithms
- **Interpolation.hpp**: Interpolation method implementation
- **LookupTable.hpp**: Multi-dimensional lookup table with interpolation support
- **MathIO.hpp**: I/O functions of math
- **MathUtility.hpp**: Utility functions for math
- **Numbers.hpp**: Define of const numbers in math
- **PolynomialFit.hpp**: polynomial fit

#### model/
- **Graph.hpp**: basic graph models

#### thread/
- **LockFreeBitSet.hpp**: Lock free bit set, modified from folly concurrent bit set
- **LockFreeHashMap.hpp**: Lock free hash map, modified from folly atomic hash map
- **MapReduce.hpp**: A header only map reduce library on single-machine platform
- **TaskFlow.hpp**: A header only task flow library that run tasks based on dependency parallelly
- **ThreadPool.hpp**: A header only thread pool implementation
- **ThreadSafeContainer.hpp**: Thread-safe container implementations
- **Utility.hpp**: Thread related utilities

#### tools/
- **Color.hpp**: Color related functions
- **FileSystem.hpp**: File system related functions
- **Format.hpp**: String formatting
- **Hash.hpp**: hash related functions
- **ImgIO.hpp**: image i/o
- **Log.hpp**: A header only log library
- **Parser.hpp**: Parser functions
- **ProgramOptions.hpp**: A header only program option library
- **StringHelper.hpp**: String related functions
- **Tools.hpp**: Some tools
- **Units.hpp**: Unit definition and functions

#### topology/
- **Common.hpp**: Common define of topology
- **IndexGraph.hpp**: Model of index graph concept and graph related algorithms

#### traits/
- **BoostGraphTraits.hpp**: boost graph adapter

#### tree/
- **BVH.hpp**: Model of bounding volume hierarchy tree concept and related algorithms
- **BVHUtilityMT.hpp**: BVH tree utility with multi-threads support
- **Builder.hpp**: General tree building helper classes
- **BuilderMT.hpp**: General tree building helper classes in multi-threads
- **IO.hpp**: I/O function for trees
- **KdTree.hpp**: Model of k-dimensional tree and related algorithms
- **KdTreeUtilityMT.hpp**: kd tree utility with multi-threads support
- **QuadTree.hpp**: Model of quad tree concept and related algorithms
- **QuadTreeUtilityMT.hpp**: Quad tree utility with multi-threads support
- **RandomTree.hpp**: Model of random tree
- **RectTree.hpp**: Model of R tree and related algorithms
- **Varification.hpp**: Varification for tree build corectness and quality check

#### utils/
- **Index.hpp**: index utils implement based on https://github.com/verilog-to-routing/tatum
- **LinearMap.hpp**: linear map utils implement based on https://github.com/verilog-to-routing/tatum
- **Version.hpp**: version utility
- **ZipView.hpp**: simple implementation of c++23's std::ranges::zip
<!-- AUTO_DOCS_END -->

## Test
Generic library unit test is written with Boost Unit Test Framework, source code is in folder `test/`.<br>  
build:<br> 
export EIGEN_PATH= ...  <br> 
export BOOST_PATH= ...  <br> 
cd generic/test         <br> 
mkdir build && cd build <br> 
cmake .. && make        <br> 
<br>

## License
See `LICENSE`.<br>
