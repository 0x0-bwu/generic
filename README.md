# Generic Library

## Introduction
A Header only library for geometry handling and others.

Implement with C++17 features, templated design for easy usage and embedding.

It has multiple 2d/3d geometry models, and model adaption with Boost Polygon/Geometry Library, so the algorithms in these two libraries can be applied directly with those models.

It offers lots of basic geometry test utility as well as high level algorithms like multi-layer geometry connectivity extraction, occupancy grid mapping and Delaunay triangulation construction and refinement, etc.

It offers templated tree structures like KdTree, QuadTree, RTree, and BVH Tree for geometry indexing and other algorithms.

It has three concurrency patterns ThreadPool, TaskFlow and MapReduce to support parallel computing in generic library, these patterns can also be applied to parallel programming for others easily.

It also has lots of other header-only utilities like log and program option for convenient usage.

## Dependency
Some components depend on Boost library(header-only), need include boost library path when compiling.

Some function implementation depends on third part library binaries with the Macro control defined in common/Macros.hpp :
if defined BOOST_SERIALIZATION_SUPPORT, need link with libboost_serialization.
if defined BOOST_FILESYSTEM_SUPPORT, need link with libboost_filesystem.
if defined BOOST_GIL_IO_PNG_SUPPORT, need link with libpng.

## Namespace Description
+ generic                       : top namespace 
    + color                     : color related definitions, classes and functions  
    + common                    : common definitions, classes, and functions  
    + filesystem                : filesystem related functions  
    + format	                : string format  
    + geometry                  : models, boost polygon/geometry adaption, geometry algorithms   
        + boolean     	        : geometry Boolean operation APIs    
	    + tri			        : triangulation data, construction, and refinement     
    + log 				        : log library  
    + math				        : math numbers, classes, and functions  
	    + la			        : linear algebra  
    + parser 			        : parser functions  
    + program_options		    : program option library  
    + str				        : string related functions  
    + thread			        : concurrency patterns, thread pool  
	    + mapreduce		        : pattern of map reduce     
	    + taskflow		        : pattern of task flow  
    + tools				        : tools and basic utilities  
    + topology			        : graph related  
    + tree				        : tree structures and algorithms  
    + unit				        : unit definition and convertion

## Documents
You can generate documents of this library by run `doxygen Doxyfile` in current folder, the documents will be generated in folder `docs/`.  

## Test
Generic library unit test is written with Boost Unit Test Framework, source code is in folder `test/`.  
The test binary can be built with SCons by run `scons` in current folder, the binary will be built in folder `build/`.

## LICENSE
See `LICENSE`
