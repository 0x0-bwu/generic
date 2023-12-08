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

if enable BOOST_SERIALIZATION_SUPPORT, need link with libboost_serialization.<br>
if enable BOOST_GIL_IO_PNG_SUPPORT, need link with libpng.<br>

## Namespace Description
+ `generic`                         : top namespace<br> 
    + `color`                       : color related definitions, classes and functions<br>  
    + `common`                      : common definitions, classes, and functions<br>  
    + `filesystem`                  : filesystem related functions<br>  
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
