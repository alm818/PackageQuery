cmake_minimum_required(VERSION 3.21)

option(CXX "enable C++ compilation" ON)
if (CXX)
  enable_language(CXX)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILE_FEATURES cxx_std_17)

project(PackageQuery VERSION 1.0.0)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(EXECUTABLE_NAME "main")

# Matplotlib
# include(GNUInstallDirs)
# set(PACKAGE_NAME matplotlib_cpp)
# set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/${PACKAGE_NAME}/cmake)

# add_library(matplotlib_cpp INTERFACE)
# find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
# target_link_libraries(matplotlib_cpp INTERFACE
#   Python3::Python
#   Python3::Module
# )
# find_package(Python3 COMPONENTS NumPy)
# if(Python3_NumPy_FOUND)
#   target_link_libraries(matplotlib_cpp INTERFACE
#     Python3::NumPy
#   )
# else()
#   target_compile_definitions(matplotlib_cpp INTERFACE WITHOUT_NUMPY)
# endif()
# install(
#   TARGETS matplotlib_cpp
#   EXPORT install_targets
# )

# For window
# set(PYTHONPACKAGE "C:/Users/xuana/AppData/Local/Packages/PythonSoftwareFoundation.Python.3.9_qbz5n2kfra8p0/LocalCache/local-packages/Python39/site-packages")
# set(PYTHONHOME "C:/Program Files/WindowsApps/PythonSoftwareFoundation.Python.3.9_3.9.2032.0_x64__qbz5n2kfra8p0")
# include_directories(${PYTHONHOME}/include)
# include_directories(${PYTHONPACKAGE}/numpy/core/include)
# link_directories(${PYTHONHOME}/libs)
# add_definitions(-DMATPLOTLIBCPP_PYTHON_HEADER=Python.h)

find_package(PythonLibs REQUIRED)
INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS} )
link_directories( ${PYTHON_LIBRARIES})

# MAIN
# set(MAIN_NAME "main_config")
# set(MAIN_NAME "main_dlv")
# set(MAIN_NAME "main_sketchrefine")

# TEST
# set(MAIN_NAME "test_final")
# set(MAIN_NAME "test_det_bound")
set(MAIN_NAME "test_dlv")

message("PROJECT NAME: " ${PROJECT_NAME})
message("SOURCE FOLDER: " ${CMAKE_SOURCE_DIR})
message("BUILD FOLDER: " ${CMAKE_BINARY_DIR})
message("CMAKE MODULE FOLDER: " ${CMAKE_MODULE_PATH})
message("EXECUTABLE NAME: " ${MAIN_NAME})

include_directories("${PROJECT_SOURCE_DIR}")
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)

file(GLOB all_SRCS
  "${CMAKE_SOURCE_DIR}/runnable/${MAIN_NAME}.cpp"
  "${PROJECT_SOURCE_DIR}/include/pb/util/*.h"
  "${PROJECT_SOURCE_DIR}/include/pb/util/*.cpp"
  "${PROJECT_SOURCE_DIR}/include/pb/core/*.h"
  "${PROJECT_SOURCE_DIR}/include/pb/core/*.cpp"
  "${PROJECT_SOURCE_DIR}/include/pb/det/*.h"
  "${PROJECT_SOURCE_DIR}/include/pb/det/*.cpp"
  "${PROJECT_SOURCE_DIR}/include/pb/sto/*.h"
  "${PROJECT_SOURCE_DIR}/include/pb/sto/*.cpp"
  "${PROJECT_SOURCE_DIR}/include/pb/lib/*.h"
  "${PROJECT_SOURCE_DIR}/include/pb/lib/*.hpp"
  "${PROJECT_SOURCE_DIR}/include/pb/lib/*.cpp"
  "${PROJECT_SOURCE_DIR}/include/pb/tmp/*.hpp"
  "${PROJECT_SOURCE_DIR}/include/pb/tmp/*.cpp"
)

#OpenMP for Windows
# find_package(OpenMP)
# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "-g")
endif(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

# Conan
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

# Gurobi
find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIR})

# Cplex
find_package(Cplex)
include_directories(${CPLEX_INCLUDE_DIRS})

add_executable(${EXECUTABLE_NAME} ${all_SRCS})

# OpenMP for Ubuntu
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

if(CXX)
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
endif()

target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${CONAN_LIBS} ${GUROBI_LIBRARY} ${CMAKE_DL_LIBS})
target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${CPLEX_LIBRARIES})

# target_link_libraries(${EXECUTABLE_NAME} PRIVATE matplotlib_cpp)
# set_target_properties(${EXECUTABLE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/resource/matplotlib")

# Matplotlib
# install(FILES
#   "${PROJECT_SOURCE_DIR}/include/matplotlibcpp.h"
#   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# install(EXPORT install_targets
#   FILE
#     ${PACKAGE_NAME}Targets.cmake
#   NAMESPACE
#     ${PACKAGE_NAME}::
#   DESTINATION
#     ${INSTALL_CONFIGDIR}
# )

message("CMAKE CURRENT SOURCE DIR: " ${CMAKE_CURRENT_SOURCE_DIR})
message("CMAKE CURRENT BINARY DIR: " ${CMAKE_CURRENT_BINARY_DIR})
message("INSTALL CONFIGDIR: " ${INSTALL_CONFIGDIR})

# include(CMakePackageConfigHelpers)
# configure_package_config_file(
#   ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PACKAGE_NAME}Config.cmake.in
#   ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Config.cmake
#   INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
# )
# install(FILES
#   ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Config.cmake
#   DESTINATION ${INSTALL_CONFIGDIR}
# )
