cmake_minimum_required(VERSION 3.21)

option(CXX "enable C++ compilation" ON)
if (CXX)
  enable_language(CXX)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)

project(StochasticPackageQuery VERSION 1.0.0)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DDEBUG)
endif()

set(EXECUTABLE_NAME main)
set(SRC src)
set(TEST test)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Find Gurobi
find_path(GUROBI_INCLUDE_DIR NAMES gurobi_c++.h PATHS "$ENV{GUROBI_HOME}/include" NO_DEFAULT_PATH)
find_library(GUROBI_LIBRARY NAMES gurobi110 PATHS "$ENV{GUROBI_HOME}/lib" NO_DEFAULT_PATH)
find_library(GUROBI_CXX_LIBRARY NAMES gurobi_c++ PATHS "$ENV{GUROBI_HOME}/lib" NO_DEFAULT_PATH)
include_directories(${GUROBI_INCLUDE_DIR})

find_package(Threads REQUIRED)
find_package(fmt REQUIRED)
find_package(PostgreSQL REQUIRED)
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(CCfits REQUIRED)

if(APPLE)
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)

    if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
find_package(OpenMP REQUIRED)

include_directories(
  ${PROJECT_SOURCE_DIR}/${SRC}
)

set(MAIN_NAME main)
#set(MAIN_NAME test_reproduce)

file(GLOB ALL_SRCS
  "${PROJECT_SOURCE_DIR}/${TEST}/${MAIN_NAME}.cpp"
  "${PROJECT_SOURCE_DIR}/${SRC}/util/*.hpp"
  "${PROJECT_SOURCE_DIR}/${SRC}/util/*.cpp"
  "${PROJECT_SOURCE_DIR}/${SRC}/det/*.hpp"
  "${PROJECT_SOURCE_DIR}/${SRC}/det/*.cpp"
  "${PROJECT_SOURCE_DIR}/${SRC}/core/*.hpp"
  "${PROJECT_SOURCE_DIR}/${SRC}/core/*.cpp"
  "${PROJECT_SOURCE_DIR}/${SRC}/lib/*.hpp"
  "${PROJECT_SOURCE_DIR}/${SRC}/lib/*.cpp"
  "${PROJECT_SOURCE_DIR}/${SRC}/tmp/*.hpp"
  "${PROJECT_SOURCE_DIR}/${SRC}/tmp/*.cpp"
)

add_executable(${EXECUTABLE_NAME} ${ALL_SRCS})
target_link_libraries(${EXECUTABLE_NAME} OpenMP::OpenMP_CXX)
target_link_libraries(${EXECUTABLE_NAME} fmt::fmt)
target_link_libraries(${EXECUTABLE_NAME} Eigen3::Eigen)
target_link_libraries(${EXECUTABLE_NAME} PostgreSQL::PostgreSQL)
target_link_libraries(${EXECUTABLE_NAME} boost::boost)
target_link_libraries(${EXECUTABLE_NAME} ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY})