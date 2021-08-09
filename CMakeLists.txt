cmake_minimum_required(VERSION 3.17)
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(dpsrvf)
set(CMAKE_CXX_STANDARD 14)
set (SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME})
include(GNUInstallDirs)
include_directories("${SRC_DIR}")

add_library(ldpsrvf ${SRC_DIR}/dpmatch.cpp ${SRC_DIR}/dpmatch.h)

set(PYTHON_MODULE_NAME dpsrvf)
add_subdirectory(lib/pybind11)
pybind11_add_module(${PYTHON_MODULE_NAME} MODULE ${SRC_DIR}/pydpmatch.cpp ${SRC_DIR}/dpmatch.cpp ${SRC_DIR}/dpmatch.h)
set (PYBIND_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/pybind11/include/pybind11)
target_include_directories(${PYTHON_MODULE_NAME} PRIVATE ${PYBIND_INCLUDE_DIR})
