cmake_minimum_required(VERSION 3.18...3.20)
project(generic_test)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
message("Build type: " ${CMAKE_BUILD_TYPE})

find_package(Threads REQUIRED)

get_filename_component(PARENT_DIR ../ ABSOLUTE)
include_directories(${PARENT_DIR}/ ${PARENT_DIR}/../3rdparty/include)
add_executable(${PROJECT_NAME} test/UnitTest.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads)