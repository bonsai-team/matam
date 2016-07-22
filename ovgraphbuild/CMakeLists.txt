cmake_minimum_required (VERSION 2.8)

find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project (ovgraphbuild C CXX)

include_directories ("lib/seqan/include")
set (PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

file(GLOB SOURCE_FILES "${PROJECT_SOURCE_DIR}/*.cpp" "${PROJECT_SOURCE_DIR}/*.h")

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros")

add_executable (ovgraphbuild ${SOURCE_FILES})

find_library (LIBRT_LIB NAMES rt librt)
find_library (PTHREAD_LIB NAMES pthread)

target_link_libraries (ovgraphbuild ${LIBRT_LIB} ${PTHREAD_LIB})

find_package(ZLIB)
if (ZLIB_FOUND)
    include_directories (${ZLIB_INCLUDE_DIRS})
    target_link_libraries (ovgraphbuild ${ZLIB_LIBRARIES})
endif ()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#set_property(TARGET ovgraphbuild PROPERTY CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)
