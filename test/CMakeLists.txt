cmake_minimum_required(VERSION 3.14)
project(unit_test.exe)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
message("-- Build Type: " ${CMAKE_BUILD_TYPE})

add_compile_options(-Wall)
add_compile_options(-Wno-deprecated-declarations)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-fcompare-debug-second) #hide note message
endif()

# Boost library
set(BOOST_INCLUDE_PATH $ENV{BOOST_PATH}/include)
set(BOOST_LIBRARY_PATH $ENV{BOOST_PATH}/lib)

## Ignore BOOST deprecated headers
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS")
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS")

## Link Eigen
if(DEFINED ENV{EIGEN_PATH})
	include_directories($ENV{EIGEN_PATH})
else()
	find_package(Eigen3 REQUIRED)
	include_directories(${EIGEN3_INCLUDE_DIR})
endif()

find_package(PNG REQUIRED)
find_package(Threads REQUIRED)

if(DEFINED ENV{EIGEN_PATH})
	include_directories($ENV{EIGEN_PATH})
endif()

get_filename_component(GENERIC_LIBRARY_DIR ../../ ABSOLUTE)
include_directories(${GENERIC_LIBRARY_DIR} ${PNG_INCLUDE_DIR} ${BOOST_INCLUDE_PATH})

if(EXISTS ${CMAKE_SOURCE_DIR}/test.cpp)
    add_executable(test.exe test.cpp)
    target_link_libraries(test.exe PRIVATE Threads::Threads)
endif()
add_executable(${PROJECT_NAME} UnitTest.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads ${PNG_LIBRARY})